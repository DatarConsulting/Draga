"""Exact compact all-word PMI model; real two-state competitors, no commutation.

Scalar polynomial optimization. The sign group acts only by state similarities.
No hypothesis that a minimizer respects alphabet exchange is imposed.
"""
import itertools
import sympy as s

def target():
    A=[]
    for j in (1,3):
        a=s.zeros(6);a[0,j]=1;a[j,j+1]=1;a[j+1,5]=1;A.append(a)
    alpha=s.eye(6)[:,0];beta=s.Matrix([0,0,1,0,1,s.Rational(1,100)])
    PA=alpha*alpha.T;QA=beta*beta.T
    for _ in range(3):
        PA=alpha*alpha.T+sum((a.T*PA*a for a in A),s.zeros(6))
        QA=beta*beta.T+sum((a*QA*a.T for a in A),s.zeros(6))
    assert PA==s.diag(1,1,1,1,1,2)
    return A,alpha,beta,PA,QA

def build(stein_equalities=False,trace_cuts=False):
    xs=[];chars=[]
    def var(name,char=0):
        a=s.Symbol(name);xs.append(a);chars.append(char);return a
    B=[]
    for a in 'ab':B.append(s.Matrix(2,2,lambda i,j:var(a+str(i)+str(j),0 if i==j else 3)))
    g=var('g',1);gamma=s.Matrix([g,0]);v=s.Matrix([var('v0',1),var('v1',2)])
    z00=var('z00');z01=var('z01',3);z11=var('z11');Z=s.Matrix([[z00,z01],[z01,z11]])
    R=s.zeros(6,2);R[0,0]=-g
    for i in range(1,6):
        for j in range(2):R[i,j]=var('r'+str(i)+str(j),1<<j)
    W=s.zeros(6,2)
    for i in range(5):
        for j in range(2):W[i,j]=var('w'+str(i)+str(j),1<<j)
    W[5,:]=(v/100).T
    u=var('u');A,alpha,beta,PA,QA=target()
    P=PA.row_join(R).col_join(R.T.row_join(s.eye(2)))
    Q=QA.row_join(3*W).col_join((3*W).T.row_join(9*Z))
    O=s.eye(2)-gamma*gamma.T-sum((b.T*b for b in B),s.zeros(2))
    C=Z-v*v.T-sum((b*Z*b.T for b in B),s.zeros(2))
    ER=R+alpha*gamma.T-sum((a.T*R*b for a,b in zip(A,B)),s.zeros(6,2))
    EW=W-beta*v.T-sum((a*W*b.T for a,b in zip(A,B)),s.zeros(6,2))
    eq=[s.expand(x) for E in (ER,EW) for x in E if x!=0]
    old=s.Rational(2839,2500)**2;upper=s.Rational(6533,5000)**2
    # An explicit redundant ball: all unscaled squared coordinates sum < 24.
    matrices=[('observability',O,[1,2]),('reachability',C,[1,2]),
              ('Z',Z,[1,2]),('I_minus_Z',s.eye(2)-Z,[1,2]),
              ('P',P,[0]*6+[1,2]),('Q',Q,[0]*6+[1,2]),
              ('error',u*P-P*Q*P,[0]*6+[1,2]),
              ('ball',s.Matrix([24-sum(x*x for x in xs)]),[0]),
              ('lower_v08',s.Matrix([u-old]),[0]),
              ('upper',s.Matrix([upper-u]),[0])]
    matrices=[(n,m.applyfunc(s.expand),c) for n,m,c in matrices]
    assert len(xs)==35 and len(eq)==20
    if stein_equalities:
        eq += [s.expand(E[i,j]) for E in (O,C) for i in range(2) for j in range(i,2)]
        matrices=[item for item in matrices if item[0] not in ('observability','reachability')]
        # Valid two-state Schur determinant cuts. For either normalized
        # radial transition (Ba +/- Bb)/sqrt(2), its determinant lies [-1,1].
        for sign in (1,-1):
            det=s.expand((B[0]+sign*B[1]).det())
            for side in (1,-1):
                matrices.append(('radial_det_'+str(sign)+'_'+str(side),s.Matrix([2+side*det]),[0]))
    if trace_cuts:
        # For PSD P,Q satisfying uP-PQP >= 0, the spectrum of PQ lies
        # in [0,u]. These redundant consequences expose the norm condition
        # at degrees two and four to a low-order moment relaxation.
        PQ=P*Q
        trace=s.expand(s.trace(PQ))
        trace_square=s.expand(sum(PQ[i,j]*PQ[j,i] for i in range(8) for j in range(8)))
        matrices += [('trace_norm',s.Matrix([8*u-trace]),[0]),
                     ('trace_square_norm',s.Matrix([s.expand(u*trace-trace_square)]),[0])]
    return {'variables':xs,'characters':chars,'B':B,'gamma':gamma,'v':v,'Z':Z,
            'R':R,'W':W,'u':u,'P':P,'Q':Q,'equalities':eq,'matrices':matrices,
            'old_squared':old,'upper_squared':upper,'ball_radius_squared':24,
            'stein_equalities':stein_equalities,'trace_cuts':trace_cuts}

def monomials(n,degree):
    return [a for d in range(degree+1) for a in itertools.combinations_with_replacement(range(n),d)]

def character(m,chars):
    c=0
    for i in m:c^=chars[i]
    return c

def polynomial(p,xs):
    out={}
    for ex,c in s.Poly(p,*xs).terms():
        if c:out[tuple(i for i,e in enumerate(ex) for _ in range(e))]=c
    return out

def multiply(a,b):return tuple(sorted(a+b))
