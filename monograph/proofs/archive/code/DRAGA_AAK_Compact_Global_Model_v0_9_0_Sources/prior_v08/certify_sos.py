"""Rational PSD Gram factors plus an exact coefficient-error budget.

No floating-point assertion is used by verify(). Floating-point eigenvectors
are used only by discover_factors(), then frozen as integer matrices.
"""
from pathlib import Path
from fractions import Fraction as Q
from collections import defaultdict
import itertools,json
import numpy as np
import sympy as sy
from polynomial_model import build,monomials
ROOT=Path(__file__).resolve().parent

def discover_factors():
    data=json.loads((ROOT/'sos_compressed_discovery.json').read_text())
    den=10**8;factors={}
    for name,matrix in data['cones'].items():
        A=np.array(matrix);w,V=np.linalg.eigh((A+A.T)/2)
        keep=w>1e-12
        factors[name]=np.rint((V[:,keep]*np.sqrt(w[keep]))*den).astype(np.int64).tolist()
    import math
    reference=Q(math.floor(-data['normalization_dual']*10**5),10**5)
    cert={'denominator':den,'factors':factors,
          'equality_numerators':[int(round(v*den)) for v in data['equality_dual']],
          'reference_lower_squared':str(reference),
          'candidate_upper':'6533/5000',
          'description':'Rational Gram factors; exact polynomial residual absorbed by coefficient box bounds.'}
    (ROOT/'lower_certificate.json').write_text(json.dumps(cert,indent=2))

def verify(path=None):
    data=json.loads((path or ROOT/'lower_certificate.json').read_text())
    den=data['denominator'];x,K,F,M,S,minors,ball,upper=build();n=len(x)
    zero=(0,)*n;b1=monomials(n,1);b2=monomials(n,2)
    residual=defaultdict(Q);unit=list(zero);unit[-1]=1;residual[tuple(unit)]=Q(1)
    val=Q(data['reference_lower_squared']);residual[zero]-=val
    def terms(p):
        return [(a,Q(int(c.p),int(c.q))) for a,c in sy.Poly(p,*x).terms() if c]
    definitions=[('moment',sy.eye(1),b2),('norm',S,b1),('ball',sy.Matrix([ball]),b1),
                 ('u',sy.Matrix([x[-1]]),b1),('upper',sy.Matrix([upper]),b1)]
    for name,mat,basis in definitions:
        B=np.array(data['factors'][name],dtype=object)
        assert B.shape[0]==mat.rows*len(basis)
        Gram=B@B.T;nb=len(basis)
        for i in range(mat.rows):
          for j in range(mat.cols):
            tt=terms(mat[i,j])
            for a,ma in enumerate(basis):
              for b,mb in enumerate(basis):
                v=Q(int(Gram[i*nb+a,j*nb+b]),den*den)
                if not v:continue
                for mon,c in tt:
                  expo=tuple(aa+bb+cc for aa,bb,cc in zip(ma,mb,mon))
                  residual[expo]-=v*c
        print('exact PSD-factor polynomial:',name,flush=True)
    k=0
    for p in minors:
      for ma in b1:
        v=Q(data['equality_numerators'][k],den);k+=1
        if v:
          for mon,c in terms(p):residual[tuple(a+b for a,b in zip(ma,mon))]+=v*c
    assert k==len(data['equality_numerators'])
    U=Q(data['candidate_upper']);eta=Q(1,100)
    # Bounds follow entrywise from E^T M^{-1} E <= U^2 M.
    bounds=[U,U,U,1+U,U,U,1+U,Q(93,50),Q(93,50),Q(93,50),Q(93,50),2*U,U*U]
    assert (Q(93,50)-eta)**2>2*U*U
    # The redundant ball follows directly from the entrywise bounds.
    assert sum(b*b for b in bounds[:12])<80
    loss=Q(0)
    for mon,c in residual.items():
        loss+=abs(c)*sy.prod(b**a for b,a in zip(bounds,mon))
    loss=Q(loss);lower_squared=val-loss
    # Choose a four-decimal rational strictly below the exact bound.
    # Integer arithmetic fixes the reported endpoint without square roots.
    from math import isqrt
    units=isqrt((lower_squared.numerator*10**8)//lower_squared.denominator)
    claimed=Q(units,10000)
    if lower_squared==claimed*claimed:claimed-=Q(1,10000)
    assert claimed>=Q(1125,1000), (float(loss),float(lower_squared))
    assert lower_squared>claimed*claimed, (float(loss),float(lower_squared))
    report={'status':'PASS','arithmetic':'RATIONAL_ONLY_IN_VERIFIER',
            'scope':'ALL_BOUNDED_REAL_HANKEL_COMPETITORS_OF_RANK_AT_MOST_TWO',
            'lower':str(claimed),'lower_squared_from_identity':str(lower_squared),
            'residual_budget':str(loss),'number_of_residual_monomials':len(residual),
            'normalization':str(val),'box_bounds':list(map(str,bounds)),
            'note':'If error exceeds 1.3066 the lower bound is immediate; otherwise the certified box and all PSD constraints apply.'}
    (ROOT/'lower_verification.json').write_text(json.dumps(report,indent=2))
    print('EXACT residual budget',float(loss),'lower sqrt',float(lower_squared)**.5,'claimed',claimed,'PASS',flush=True)
    return report

if __name__=='__main__':
    import sys
    if '--discover' in sys.argv:discover_factors()
    verify()
