#!/usr/bin/env python3
"""Exact finite certificates for DRAGA v0.7.0 (not a general SOS solver).

Universal convergence and the complete exact search are proved in the note.
This replay checks a genuine polynomial refutation, full-tail Gramian upper
certificates, a multidegree strict-gap witness and the interval contraction.
"""
import argparse
import itertools
import json
import platform
from pathlib import Path
import sympy as S


def require(ok, message):
    if ok != True: raise AssertionError(message)


def matrix_json(M):
    return [[str(M[i,j]) for j in range(M.cols)] for i in range(M.rows)]


def positive_minors(M):
    require(M==M.T,'Symmetry')
    result=[S.factor(M[:i,:i].det()) for i in range(1,M.rows+1)]
    require(all(x>0 for x in result),'Exact Sylvester certificate')
    return [str(x) for x in result]


def target(t=S.Rational(1,100)):
    A=[]
    for i in range(2):
        M=S.zeros(6); j=1+2*i
        M[0,j]=1; M[j,5]=1; M[j,j+1]=1; M[j+1,5]=t
        A.append(M)
    return A,S.eye(6)[:,0],S.eye(6)[:,5]


def evaluate(A,alpha,beta,w):
    row=alpha.T
    for i in w: row=row*A[i]
    return (row*beta)[0]


def polynomial_refutation():
    # Principal block on empty, a, aa of the SAME binary target.
    tau=S.Rational(1,100); ell=S.Rational(49,50); eta=S.Rational(1,100)
    F=S.Matrix([[0,0,1],[0,1,tau],[1,tau,0]])
    z=S.symbols('z0:5'); v=S.Matrix(S.symbols('v0:3'))
    Y=S.Matrix(3,3,S.symbols('y0:9'))
    G=S.Matrix([[z[i+j] for j in range(3)] for i in range(3)])
    E=F-G; q=G*v; p=(v.T*v)[0]-1
    slack=ell**2*S.eye(3)-E.T*E-Y*Y.T
    M=F.T*F-(ell**2+eta)*S.eye(3)
    minors=positive_minors(M)
    L,D=M.LDLdecomposition(hermitian=False)
    require(L*D*L.T==M,'Exact LDL decomposition')
    square_polys=list(L.T*v)+list(Y.T*v)
    weights=[D[i,i]/eta for i in range(3)]+[1/eta]*3
    require(all(w>0 for w in weights),'Positive rational SOS weights')
    sos=sum(w*s*s for w,s in zip(weights,square_polys))
    ideal=p+((v.T*slack*v)[0]-2*(E*v).dot(q)-q.dot(q))/eta
    require(S.expand(sos+ideal+1)==0,'Exact real Nullstellensatz identity')
    require(S.expand(sos-ideal+1)!=0,'Negative control: changed sign rejected')
    # Explicit constraints and multipliers make the certificate replayable.
    equations=[p]+list(q)
    multipliers=[S.Integer(1)]+[-(2*(E*v)[i]+q[i])/eta for i in range(3)]
    for i in range(3):
        for j in range(i,3):
            equations.append(slack[i,j])
            multipliers.append((1 if i==j else 2)*v[i]*v[j]/eta)
    require(S.expand(sos+sum(a*b for a,b in zip(equations,multipliers))+1)==0,
            'Serialized multiplier identity')
    return {'status':'PASS','block_words':['empty','a','aa'],'threshold':str(ell),
            'type':'EXACT_WEIGHTED_SOS_PLUS_IDEAL_EQUALS_MINUS_ONE',
            'variables':[str(x) for x in z]+[str(x) for x in v]+[str(x) for x in Y],
            'squares':[str(S.expand(x)) for x in square_polys],
            'positive_weights':[str(w) for w in weights],
            'equations':[str(S.expand(x)) for x in equations],
            'multipliers':[str(S.expand(x)) for x in multipliers],
            'positive_matrix_minors':minors}


def coherent_system(alphabet_size,N,r,coefficient,threshold):
    """Generate the full quadratic equality system for a finite lower test.

    This generator does not solve the resulting polynomial system.
    """
    words=[w for length in range(N+1)
           for w in itertools.product(range(alphabet_size),repeat=length)]
    coefficient_words=[w for length in range(2*N+1)
                       for w in itertools.product(range(alphabet_size),repeat=length)]
    z=dict(zip(coefficient_words,S.symbols('c0:'+str(len(coefficient_words)))))
    F=S.Matrix([[coefficient(u+v) for v in words] for u in words])
    G=S.Matrix([[z[u+v] for v in words] for u in words])
    d=len(words); nullity=max(d-r,0)
    V=S.Matrix(d,nullity,S.symbols('v0:'+str(d*nullity)))
    Y=S.Matrix(d,d,S.symbols('y0:'+str(d*d)))
    equations=list(G*V)+list(V.T*V-S.eye(nullity))
    slack=threshold**2*S.eye(d)-(F-G).T*(F-G)-Y*Y.T
    equations.extend(slack[i,j] for i in range(d) for j in range(i,d))
    return {'words':words,'G':G,'F':F,'V':V,'Y':Y,'equations':equations}


def strict_lower():
    t=S.symbols('t'); a=1+t-t**3; b=2*t*t+t-2
    poly=S.Poly(2*a*a-b*b,t)
    lo=S.Rational(2627,2500); hi=S.Rational(105083,100000)
    require(a.subs(t,lo)>0 and b.subs(t,lo)>0 and poly.eval(lo)>0,'Lower root guard')
    require(a.subs(t,hi)>0 and b.subs(t,hi)>0 and poly.eval(hi)<0,'Upper root guard')
    require(poly.count_roots(lo,hi)==1,'Exact root isolation')
    sqrt2_upper=S.Rational(707107,500000)
    require(sqrt2_upper**2>2,'Rational upper bound on sqrt(2)')
    perturb=sqrt2_upper/100
    lower=S.Rational(20733,20000) # 1.03665
    require(lo-perturb>lower,'Transferred finite-section lower certificate')
    require(lower>1+perturb,'Strict gap above full-target sigma3')
    return {'root_polynomial':str(poly.as_expr()),'root_interval':[str(lo),str(hi)],
            'certified_lower':str(lower),'sigma3_upper':str(1+perturb),
            'scope':'COMPLEX_AND_REAL_COMPETITORS',
            'proof':'Written finite-block angle proof + perturbation; not expanded into SOS.'}


def gramian_upper(parameters,U):
    q,p,c0,c1=map(S.Rational,parameters)
    A,alpha,beta=target()
    B=[S.Matrix([[0,1],[-q,p]]),S.zeros(2)]
    gamma=S.Matrix([1,0]); nu=S.Matrix([c0,c1])
    X=S.Matrix(2,2,(S.eye(4)-S.kronecker_product(B[0].T,B[0].T)).inv()
               *S.eye(2).reshape(4,1))
    require(X-B[0].T*X*B[0]==S.eye(2),'Exact candidate stability identity')
    stability=positive_minors(X)
    D=[S.diag(a,b) for a,b in zip(A,B)]
    c=alpha.col_join(-gamma); v=beta.col_join(nu)
    K=sum((S.kronecker_product(M,M) for M in D),S.zeros(64))
    # Solve both all-word Stein equations exactly.
    P=S.Matrix(8,8,(S.eye(64)-K.T).inv()*(c*c.T).reshape(64,1))
    Q=S.Matrix(8,8,(S.eye(64)-K).inv()*(v*v.T).reshape(64,1))
    require(P==c*c.T+sum((M.T*P*M for M in D),S.zeros(8)),'Observability identity')
    require(Q==v*v.T+sum((M*Q*M.T for M in D),S.zeros(8)),'Reachability identity')
    pmin=positive_minors(P); qmin=positive_minors(Q)
    normmin=positive_minors(U*U*P-P*Q*P)
    return {'parameters':list(map(str,(q,p,c0,c1))),'upper':str(U),
            'B_a':matrix_json(B[0]),'B_b':matrix_json(B[1]),
            'alpha':matrix_json(gamma),'beta':matrix_json(nu),
            'stability_X':matrix_json(X),'stability_minors':stability,
            'P':matrix_json(P),'Q':matrix_json(Q),'P_minors':pmin,'Q_minors':qmin,
            'norm_slack_minors':normmin,'status':'PASS'}


def target_checks():
    A,alpha,beta=target(); W=[w for length in range(4) for w in itertools.product(range(2),repeat=length)]
    require(A[0]*A[1]==A[1]*A[0],'Target commutativity')
    for length in range(5):
        for w in itertools.product(range(2),repeat=length):
            expected=(S.Integer(1) if length==2 else S.Rational(1,100) if length==3 else 0)
            if not w or any(i!=w[0] for i in w): expected=0
            require(evaluate(A,alpha,beta,w)==expected,'Exact mixed-length coefficients')
    H=S.Matrix([[evaluate(A,alpha,beta,u+v) for v in W] for u in W])
    require(H.rank()==6,'Full finite-support target rank')
    system=coherent_system(2,1,2,lambda w:evaluate(A,alpha,beta,w),S.Rational(1,2))
    require(system['G'][0,1]==system['G'][1,0],'Repeated concatenation coefficient')
    require(all(S.Poly(p).total_degree()<=2 for p in system['equations']),
            'All lower-test equations have degree at most two')
    # A rank-one finite-section ghost need not extend to a rank-one Hankel.
    ghost=S.zeros(3); ghost[2,2]=1
    require(ghost.rank()==1 and ghost[0,0]==0,'Unextendable unary ghost fixture')
    return {'target_rank':6,'lengths':[2,3],'finite_system_N1':{
             'dimension':3,'equations':len(system['equations']),'max_degree':2},
            'negative_control':'rank-one ghost has g(empty)=0 but g(a^4)=1; no nonzero one-state extension'}


def trisection_checks():
    L,U=S.symbols('L U'); a=(2*L+U)/3; b=(L+2*U)/3
    require(S.expand(U-a-S.Rational(2,3)*(U-L))==0,'Lower-update width')
    require(S.expand(b-L-S.Rational(2,3)*(U-L))==0,'Upper-update width')
    # Illustrative certified updates on this fixed target use stored actual
    # certificates. This is not an execution of the complete exhaustive search.
    lower=S.Rational(20733,20000); upper=S.Rational(6533,5000)
    require(upper-lower==S.Rational(5399,20000),'Measured certificate interval width')
    return {'width_rule':'(2/3)^j * initial_width',
            'fixed_target_certified_interval':[str(lower),str(upper)],
            'interval_width':str(upper-lower),
            'exhaustive_search_to_arbitrary_epsilon':'PROVED_NOT_IMPLEMENTED_OR_BENCHMARKED'}


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',type=Path,default=Path('certificates.json'))
    args=ap.parse_args()
    report={'version':'0.7.0','python':platform.python_version(),'sympy':S.__version__,
            'arithmetic':'EXACT','novelty':'NOT_ESTABLISHED','independent_review':'NOT_DONE',
            'universal_proof':'WRITTEN_NOT_FORMALIZED'}
    report['polynomial_refutation']=polynomial_refutation()
    print('PASS: explicit weighted-SOS refutation and negative control',flush=True)
    report['target']=target_checks(); report['strict_lower']=strict_lower()
    print('PASS: target, coherent equations, finite structural lower certificate',flush=True)
    report['upper_seed']=gramian_upper(['58127/100000','6157/5000','5443/50000','8523/25000'],S.Rational(130875,100000))
    print('PASS: seed candidate full-word error < 1.30875',flush=True)
    report['upper_improved']=gramian_upper(['290137/500000','307573/250000','106299/1000000','171563/500000'],S.Rational(6533,5000))
    print('PASS: improved candidate full-word error < 1.3066',flush=True)
    report['trisection']=trisection_checks()
    args.output.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print('Written:',args.output,flush=True)


if __name__=='__main__':main()
