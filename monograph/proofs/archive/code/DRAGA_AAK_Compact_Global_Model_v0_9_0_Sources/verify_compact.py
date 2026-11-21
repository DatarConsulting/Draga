"""Exact rational calibration and structure checks for the compact PMI model.

This verifies a feasible point and the failure of an old local pseudo-candidate.
It does NOT certify global infeasibility at the desired lower threshold.
"""
import itertools,json,math
from pathlib import Path
import sympy as s
from compact_model import build,target,polynomial,monomials,character
ROOT=Path(__file__).resolve().parent

def matrix(x):return s.Matrix(x).applyfunc(s.Rational)
def data(M):return [[str(x) for x in row] for row in M.tolist()]
def pd(M):
    assert M==M.T
    ds=[M[:i,:i].det(method='domain-ge') for i in range(1,M.rows+1)]
    assert all(d>0 for d in ds)
    return list(map(str,ds))

def construct_witness():
    # Floating-free rational approximation of a reverse Cholesky factor.
    # Each square root is bracketed using integer arithmetic; all final
    # properties are checked exactly, so the choice need not be optimal.
    prior=json.loads((ROOT/'prior_upper_exact.json').read_text())
    B=matrix(prior['A_a']);gam=matrix([prior['alpha']]).T;nu=matrix([prior['beta']]).T
    P0=matrix(prior['P']);Q0=matrix(prior['Q']);X=matrix(prior['X'])
    PB=P0[6:,6:]+s.Rational(1,10**8)*X
    def reciprocal_sqrt(a,den=10**12):
        # floor(den / sqrt(a)) / den
        return s.Rational(math.isqrt(int(s.floor(den*den/a))),den)
    s00=reciprocal_sqrt(PB[0,0]-PB[0,1]**2/PB[1,1]);s11=reciprocal_sqrt(PB[1,1])
    s10=-PB[0,1]*s00/PB[1,1]
    # Round s10 to a rational with a small common denominator as well.
    s10=s.Rational(s.floor(s10*10**12),10**12)
    S=s.Matrix([[s00,0],[s10,s11]]);Si=S.inv()
    Bs=Si*B*S;gs=S.T*gam;vs=Si*nu/3
    targetS=s.eye(6);targetS[2,5]=-100;targetS[4,5]=-100;targetS[5,5]=100
    T=s.diag(targetS,S);Ti=T.inv();Pn=T.T*P0*T;Qn=Ti*Q0*Ti.T
    P=Pn.copy();P[6:,6:]=s.eye(2);R=P[:6,6:];W=Qn[:6,6:]/3;Z=Qn[6:,6:]/9
    d=build(False,True);sub={d['u']:s.Rational(prior['upper'])**2}
    for field,M in [('gamma',gs),('v',vs),('Z',Z),('R',R),('W',W)]:
        for expr,val in zip(d[field],M):
            if isinstance(expr,s.Symbol):sub[expr]=val
    # R[0,0] and W[5,:] are eliminated expressions, not extra variables.
    for b in d['B']:
        for expr,val in zip(b,Bs):sub[expr]=val
    assert len(sub)==len(d['variables'])
    assert d['gamma'].subs(sub)==gs and d['R'].subs(sub)==R and d['W'].subs(sub)==W
    assert all(eq.subs(sub)==0 for eq in d['equalities'])
    checks={}
    for name,M,_ in d['matrices']:
        M=M.subs(sub)
        if name=='reachability':assert M==s.zeros(2);checks[name]={'zero_exactly':True}
        else:checks[name]={'positive_leading_minors':pd(M)}
    out={'parameters':{str(x):str(sub[x]) for x in d['variables']},
         'upper':prior['upper'],'B_a':data(Bs),'B_b':data(Bs),
         'gamma':data(gs),'nu':data(3*vs),'checks':checks,
         'scope':'ALL_WORDS','status':'EXACT_PASS'}
    (ROOT/'compact_feasible_exact.json').write_text(json.dumps(out,indent=2))
    return out

def old_pseudocandidate():
    """Its radial two-state recurrence has determinant > 1 after normalization."""
    R=s.Rational
    c0,c1,c2,c3=R(55,100),R(60,100),R(109,200),R(17,40)
    # c_{n+2}=p c_{n+1}-q c_n for S=(Ba+Bb)/2.
    V=s.Matrix([[c1,-c0],[c2,-c1]]);p,q=V.inv()*s.Matrix([c2,c3])
    assert c0*c2-c1*c1!=0
    assert q>R(1,2)
    discr=p*p-4*q
    assert discr<0
    c4=p*c3-q*c2
    assert 4*c4==R(133821,120500)
    H3=s.Matrix([[c0,c1,c2],[c1,c2,c3],[c2,c3,c4]])
    assert H3.det()==0 and H3.rank()==2
    return {'status':'EXACT_PASS','p':str(p),'q':str(q),
            'radial_transition_determinant':str(2*q),'discriminant':str(discr),
            'conclusion':'The old K_star has NO bounded real rank-two Hankel extension.',
            'proof_scope':'Algebraic recurrence checks; the full proof is in the accompanying note.'}

def verify_structure():
    d=build();xs=d['variables'];ch=d['characters'];A,alpha,beta,PA,QA=target()
    for length in range(5):
        for word in itertools.product(range(2),repeat=length):
            M=s.eye(6)
            for a in word:M=M*A[a]
            expected=1 if length==2 and len(set(word))==1 else s.Rational(1,100) if length==3 and len(set(word))==1 else 0
            assert (alpha.T*M*beta)[0]==expected
            if length==4:assert M==s.zeros(6)
    for _,M,states in d['matrices']:
        for i in range(M.rows):
            for j in range(M.cols):
                assert all(character(m,ch)==states[i]^states[j] for m in polynomial(M[i,j],xs))
    for p in d['equalities']:
        cs={character(m,ch) for m in polynomial(p,xs)};assert len(cs)==1
    sizes=[sum(character(m,ch)==c for m in monomials(35,2)) for c in range(4)]
    assert sizes==[195,151,148,172]
    assert 2+1+1+2+s.trace(PA)+s.trace(QA)+d['upper_squared']**2<24
    assert (s.Rational(3)-s.Rational(6533,5000))**2>2*s.Rational(101,100)**2
    pd(2*s.Rational(101,100)**2*PA-PA*QA*PA)
    return {'variables':35,'polynomial_equalities':20,'moment_blocks':sizes,
            'dense_moment_size':666,'invariant_moments':sum(character(m,ch)==0 for m in monomials(35,4)),
            'ball_radius_squared':24,'status':'EXACT_PASS'}

def verify_exact_whitening():
    """Algebraic feasible point for the stronger Stein EQUALITY model.

Only two positive square roots are needed. Calculations take place in the
four-dimensional quotient with a^2=A2 and c^2=C2; positivity is transported
by an invertible real congruence from the rational v0.8 certificate.
"""
    prior=json.loads((ROOT/'prior_upper_exact.json').read_text())
    P=matrix(prior['P']);Q=matrix(prior['Q']);PB=P[6:,6:];QB=Q[6:,6:]
    u=s.Rational(prior['upper'])**2
    base_checks={'P':pd(P),'Q':pd(Q),'norm_slack':pd(u*P-P*Q*P),
                 'Z_upper_congruence':pd(9*PB.inv()-QB)}
    a,c=s.symbols('a c');b=-PB[0,1]/PB[1,1]
    A2=1/(PB[0,0]-PB[0,1]**2/PB[1,1]);C2=1/PB[1,1]
    assert A2>0 and C2>0
    def reduce(p):
        out=0
        for (i,j),v in s.Poly(s.expand(p),a,c).terms():
            out+=v*A2**(i//2)*C2**(j//2)*a**(i%2)*c**(j%2)
        return s.expand(out)
    S=s.Matrix([[a,0],[b*a,c]])
    Si=s.Matrix([[a/A2,0],[-b*c/C2,c/C2]])
    assert (Si*S).applyfunc(reduce)==s.eye(2)
    assert (S.T*PB*S).applyfunc(reduce)==s.eye(2)
    B=matrix(prior['A_a']);gam=matrix([prior['alpha']]).T;nu=matrix([prior['beta']]).T
    Bs=(Si*B*S).applyfunc(reduce);gs=S.T*gam;vs=(Si*nu/3).applyfunc(reduce)
    Z=(Si*QB*Si.T/9).applyfunc(reduce)
    assert (s.eye(2)-gs*gs.T-2*Bs.T*Bs).applyfunc(reduce)==s.zeros(2)
    assert (Z-vs*vs.T-2*Bs*Z*Bs.T).applyfunc(reduce)==s.zeros(2)
    assert gs[1]==0
    d=build(True,True);xs=d['variables'];chars=d['characters']
    assert len(d['equalities'])==26
    for _,M,states in d['matrices']:
        for i in range(M.rows):
            for j in range(M.cols):
                assert all(character(m,chars)==states[i]^states[j] for m in polynomial(M[i,j],xs))
    out={'status':'EXACT_PASS','a_squared':str(A2),'c_squared':str(C2),'b':str(b),
         'S_form':'[[a,0],[b*a,c]], a>0, c>0',
         'upper':prior['upper'],'base_congruence_checks':base_checks,
         'stein_equalities':26,'radial_determinant_cuts':4}
    (ROOT/'stein_feasible_algebraic.json').write_text(json.dumps(out,indent=2))
    return {k:v for k,v in out.items() if k!='base_congruence_checks'}

def main():
    structural=verify_structure();print('STRUCTURE',structural,flush=True)
    pseudo=old_pseudocandidate();print('OLD_PSEUDOCANDIDATE',pseudo,flush=True)
    feasible=construct_witness();print('COMPACT FEASIBLE',feasible['status'],feasible['upper'],flush=True)
    whitening=verify_exact_whitening();print('EXACT WHITENING',whitening['status'],flush=True)
    out={'version':'0.9.0','structure':structural,'old_pseudocandidate':pseudo,
         'compact_feasible_witness':{'status':feasible['status'],'upper':feasible['upper']},
         'stein_equality_witness':whitening,
         'global_lower_verification':'See verify_all.py and verification_report.json',
         'novelty':'NOT_ESTABLISHED','formal_proof_assistant':'NOT_USED'}
    (ROOT/'compact_verification.json').write_text(json.dumps(out,indent=2))

if __name__=='__main__':main()
