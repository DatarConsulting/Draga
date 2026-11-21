"""Replay the lower, upper and relaxation-barrier certificates exactly."""
import json,sys,itertools
from pathlib import Path
import sympy as s
from certify_sos import verify as verify_lower
from polynomial_model import build
ROOT=Path(__file__).resolve().parent

def matrix(a):return s.Matrix([[s.Rational(z) for z in row] for row in a])
def pd(A):
    assert A==A.T
    vals=[A[:j,:j].det() for j in range(1,A.rows+1)]
    assert all(v>0 for v in vals)
    return list(map(str,vals))

def verify_upper():
    data=json.loads((ROOT/'upper_exact.json').read_text())
    B=matrix(data['A_a']);assert B==matrix(data['A_b'])
    X=matrix(data['X']);assert X-2*B.T*X*B==s.eye(2);pd(X)
    A=[]
    for i in range(2):
        a=s.zeros(6);j=1+2*i
        a[0,j]=1;a[j,5]=1;a[j,j+1]=1;a[j+1,5]=s.Rational(1,100);A.append(a)
    D=[s.diag(a,B) for a in A]
    alpha=s.eye(6)[:,0];beta=s.eye(6)[:,5]
    for length in range(5):
        for word in itertools.product(range(2),repeat=length):
            product=s.eye(6)
            for letter in word:product=product*A[letter]
            expected=(s.Integer(1) if length==2 else s.Rational(1,100) if length==3 else 0)
            if not word or len(set(word))!=1:expected=0
            assert (alpha.T*product*beta)[0]==expected
            if length==4:assert product==s.zeros(6)
    c=alpha.col_join(-s.Matrix(list(map(s.Rational,data['alpha']))))
    b=beta.col_join(s.Matrix(list(map(s.Rational,data['beta']))))
    P=matrix(data['P']);Q=matrix(data['Q']);U=s.Rational(data['upper'])
    assert P==c*c.T+sum((d.T*P*d for d in D),s.zeros(8))
    assert Q==b*b.T+sum((d*Q*d.T for d in D),s.zeros(8))
    pd(P);pd(Q);mi=pd(U*U*P-P*Q*P)
    assert mi==data['norm_slack_minors']
    # Negative control: the same candidate does not meet error 1.
    try:pd(P-P*Q*P)
    except AssertionError:pass
    else:raise AssertionError('Negative control was unexpectedly accepted')
    print('EXACT all-word upper',U,'PASS',flush=True)
    return {'status':'PASS','upper':str(U),'states':2,'scope':'ALL_WORDS'}

def verify_barrier():
    x,K,F,M,S,mins,ball,upper=build();R=s.Rational
    G=s.Matrix([[R(55,100),R(60,100),R(60,100),R(109,100)],
                [R(60,100),R(109,200),R(109,200),R(85,100)],
                [R(60,100),R(109,200),R(109,200),R(85,100)],
                [R(109,100),R(85,100),R(85,100),R(133821,120500)]])
    values=[G[0,0],G[0,1],G[0,2],G[1,1],G[1,2],G[2,1],G[2,2],G[1,3],G[2,3],G[3,1],G[3,2],G[3,3]]
    sub=dict(zip(x[:12],values));assert K.subs(sub)==G
    assert G.rank()==2;assert all(p.subs(sub)==0 for p in mins)
    U=R(121,100);Z=U*U*M-(F-G).T*M.inv()*(F-G);mi=pd(Z)
    print('EXACT compressed-relaxation optimum < 1.21: PASS',flush=True)
    return {'status':'PASS','G':[[str(z) for z in row] for row in G.tolist()],
            'upper_on_relaxed_optimum':str(U),'norm_slack_minors':mi,
            'meaning':'This compressed relaxation, even solved exactly, cannot close the global interval to 0.01 with the current upper.'}

def main():
    lower=verify_lower();upper=verify_upper();barrier=verify_barrier()
    perturb_upper=s.Rational(14143,1000000)
    assert perturb_upper**2>s.Rational(2,10000)
    assert s.Rational(lower['lower'])>1+perturb_upper
    width=s.Rational(upper['upper'])-s.Rational(lower['lower'])
    report={'version':'0.8.0','lower':lower,'upper':upper,'relaxation_barrier':barrier,
            'interval_width':str(width),'requested_width':'1/100',
            'requested_closure_achieved':bool(width<=s.Rational(1,100)),
            'novelty':'NOT_ESTABLISHED','proof_assistant':'NOT_USED','independent_review':'NOT_DONE'}
    (ROOT/'verification_report.json').write_text(json.dumps(report,indent=2))
    print('INTERVAL',lower['lower'],upper['upper'],'width',width,'closed_to_0.01',report['requested_closure_achieved'],flush=True)

if __name__=='__main__':main()
