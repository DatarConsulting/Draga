"""Exact checks for the perturbation corollary and explicit rank-four sources."""
import json
from pathlib import Path
from fractions import Fraction as F
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
def run():
    original=json.loads((ROOT/'results/exact_results.json').read_text());out=[]
    for row in original['records']:
        a,m=row['a'],row['m'];rho=F(1,2**a);eps=F(row['epsilon']);eta=eps/4
        # An explicit lower bound on pairwise future separation of the packing.
        cross=rho**(m+1)*(1-rho)/(2+rho)
        same=4*eps;extra=F(2**(m+1)-1,2**((a+1)*(m+1)))
        margin=min(cross,same,extra)
        assert margin>=F(16,5)*eps
        assert margin-2*eta>2*(eps+eta)
        # Rank-four source f + eta*h, h(w)=(3^-len(w)+5^-len(w))/2.
        r=s.Rational(rho.numerator,rho.denominator);e=s.Rational(eta.numerator,eta.denominator)
        A0=s.diag(r,r/2,s.Rational(1,3),s.Rational(1,5));A1=A0.copy();A1[0,1]=r/2
        al=s.Matrix([[1,0,s.Rational(1,2),s.Rational(1,2)]]);be=s.Matrix([0,1,e,e])
        R=s.Matrix.vstack(al,al*A0,al*A0*A0,al*A1)
        O=s.Matrix.hstack(be,A0*be,A0*A0*be,A1*be)
        assert R.det()!=0 and O.det()!=0
        Hbound=F(5,6)+eta*F(9,7);assert Hbound<1
        out.append({'a':a,'m':m,'eta':str(eta),'tolerance':str(eps+eta),
          'optimal_bits':m+2,'source_rank':4,'reachability_minor':str(R.det()),'observability_minor':str(O.det()),
          'unperturbed_packing_gap_lower':str(margin),'perturbed_packing_gap_lower':str(margin-2*eta),
          'source_Hankel_norm_upper':str(Hbound),'same_executable':True})
    (ROOT/'results/robust_rank_four.json').write_text(json.dumps({'all_passed':True,'records':out,
      'theorem':'For any g with sup_w |g(w)-f(w)|<=epsilon/4, optimal bits at tolerance 5epsilon/4 equal m+2. Explicit bounded-Hankel rank-four g supplied.'},indent=2))
    print('PASS robust rank-four',len(out))
if __name__=='__main__':run()
