"""Independent consistency checks for the explicit families in the manuscript.

Exact assertions use SymPy or fractions. Floating-point Stein checks are
diagnostics, not interval certificates or proofs of global optimality.
Run: python3 verify_explicit_results.py
"""
from pathlib import Path
from fractions import Fraction as Q
from itertools import product
import json, math, platform, time
import numpy as np
import scipy
from scipy.linalg import block_diag
from scipy.optimize import brentq
import sympy as s

records=[]
start=time.time()
def record(name, kind, **data):
    records.append(dict(name=name, kind=kind, passed=True, **data))
    print('PASS',name,flush=True)

def hankel_norm(A,b,c):
    n=len(b); I=np.eye(n*n)
    L=I-sum(np.kron(a,a) for a in A)
    P=np.linalg.solve(L,np.outer(b,b).reshape(-1,order='F')).reshape(n,n,order='F')
    Qm=np.linalg.solve(L.T,np.outer(c,c).reshape(-1,order='F')).reshape(n,n,order='F')
    ev=np.linalg.eigvals(P@Qm)
    assert max(abs(ev.imag))<1e-7
    assert np.linalg.norm(L@P.reshape(-1,order='F')-np.outer(b,b).reshape(-1,order='F'))<1e-8
    return float(np.sqrt(max(0,max(ev.real))))

# Exact algebra underlying the isospectral separation.
H=s.Matrix([[1,0,1,1],[0,1,0,0],[1,0,0,0],[1,0,0,0]])
assert H.eigenvals()=={s.Integer(2):1,s.Integer(1):1,s.Integer(-1):1,s.Integer(0):1}
t=s.symbols('t'); y=t-s.Rational(1,2)
assert s.expand(4-3*(1-t*t)*(4*t-t*t)-s.Rational(51,4)*(y-s.Rational(1,17))**2-s.Rational(5,272)-3*y*y*(y+s.Rational(1,2))*(s.Rational(5,2)-y))==0
b=s.Rational(13,24); c=s.Rational(197,200); upper=s.Rational(29,25)
D=s.diag(1,1,1,1,b**4/(1-b*b)); k=s.Matrix([1,0,b,0,1]); H5=s.zeros(5); H5[:4,:4]=H
L=D*H5-c*D*k*k.T*D
minors=[]
for sign in [1,-1]:
    mm=[s.factor((upper*D+sign*L)[:j,:j].det()) for j in range(1,6)]
    assert all(v>0 for v in mm); minors.append([str(v) for v in mm])
record('Isospectral example: spectrum, polynomial identity, full-operator upper certificate','exact',minors=minors)

# Homogeneous f(w)=sum_i 1_{w=a_i^m}: all nonzero singular values are known.
fixtures=0
for d in [2,3]:
 for m in [2,3,4]:
    words=[()]+[(i,)*j for i in range(d) for j in range(1,m+1)]
    support={(i,)*m for i in range(d)}
    A=s.Matrix([[int(u+v in support) for v in words] for u in words])
    eigen=(A*A.T).eigenvals()
    assert eigen.get(s.Integer(d),0)==2
    assert eigen.get(s.Integer(1),0)==d*(m-1)
    assert eigen.get(s.Integer(0),0)==d-1
    fixtures+=1
record('Homogeneous sharp-factor fixtures: complete spectra','exact',fixtures=fixtures)

# Independent infinite-Hankel Stein norm check of the finite-parameter optimizer.
data=[]
for q in [.005,.02,.1,.25,.5,.8,.95,.99]:
    r=math.sqrt(1-q)
    def eq(u):return q*q*(1-u*u)*(u-q)**2-(1-q)**3*u*u*(q-u*u)
    u=brentq(eq,q,math.sqrt(q),xtol=1e-15)
    value=(q*q+(1-2*q)*u*u)/(q*(1-u*u))
    alpha=(1-u)/r; beta=math.sqrt(u*(1-u)/(1-q)); omega=(u-q)/(1-q)
    amp=omega*u/(beta*(1+u))
    A0=np.array([[r,0.],[0.,0.]]); A1=np.array([[0.,1.],[0.,0.]])
    As=[block_diag(A0,alpha),block_diag(A1,beta)]
    norm=hankel_norm(As,np.array([0.,1.,1.]),np.array([q,0.,-amp]))
    assert abs(norm**2-value)<2e-8,(q,norm**2,value)
    Aq=(1-q)**2*value**2+(q*q+2*q-1)*(value-q)
    quartic=Aq*Aq-4*q**3*(value-q)*(1-2*q+q*value)
    assert abs(quartic)<1e-8 and Aq>0 and 2*q/(1+q)<value<1
    data.append(dict(q=q,u=u,squared_error=norm**2,quartic_residual=quartic))
record('Exact rank-one formula versus independent infinite-operator norm','floating-point',cases=data)

# k=2 explicit boundary lift and its q^2 curvature bounds.
data=[]
for q in [.125,.05,.01,.002]:
    r=math.sqrt(1-q); yy=1/math.sqrt(2); cc=1/(2*math.sqrt(2))
    At0=np.diag([r,0.,0.]); At1=np.array([[0.,1.,0.],[0.,0.,1.],[0.,0.,0.]])
    Ag0=np.diag([r,0.]); Ag1=np.array([[0.,1.],[0.,yy]])
    A=[block_diag(At0,Ag0),block_diag(At1,Ag1)]
    norm=hankel_norm(A,np.array([0.,0.,1.,0.,1.]),np.array([q,0.,0.,-q*cc,0.]))
    lo=.75+q*q/6; hi=lo+4*q**3/3
    assert lo-2e-10<=norm**2<=hi+2e-10
    data.append(dict(q=q,squared_error=norm**2,lower=lo,upper=hi))
record('Boundary lift: finite-q curvature enclosure','floating-point',cases=data)

# Exact seven-state mixed-length family, evaluated by complete Stein systems.
data=[]
for rho in [0.,.25,.7]:
    ss=.2; eta=.01
    T=np.zeros((5,5)); T[np.arange(4),np.arange(1,5)]=1; T[4,4]=rho
    A=[block_diag(np.diag([.5,.25]),T),block_diag(np.diag([.25,.5]),np.zeros((5,5)))]
    bv=np.array([1.,1.,0.,0.,0.,-rho,1-rho*rho]); cv=np.array([ss,ss,eta,0.,0.,0.,0.])
    residual=hankel_norm([T,np.zeros((5,5))],bv[2:],cv[2:])
    assert abs(residual-abs(eta))<1e-10
    ws=[(),(0,),(0,0),(0,0,0),(0,0,0,0),(1,),(1,1)]
    def coeff(w):
        x=cv.copy()
        for a in w:x=x@A[a]
        return float(x@bv)
    finite=np.array([[coeff(u+v) for v in ws] for u in ws])
    det=np.linalg.det(finite); expected=ss**2*eta**5*(1-rho*rho)/1024
    assert abs(det-expected)<abs(expected)*1e-6
    data.append(dict(rho=rho,residual_norm=residual,rank_witness_determinant=det))
record('Mixed-length seven-to-two family: residual norm and rank witness','floating-point',cases=data)

# Exact three-state response coefficient, using the metric rather than whitening.
T=s.Matrix([[0,1],[-s.Rational(1,2),-1]]); Gm=s.Matrix([[s.Rational(8,5),s.Rational(4,5)],[s.Rational(4,5),s.Rational(12,5)]])
ell=s.Matrix([[1,0]])
assert Gm-T.T*Gm*T==ell.T*ell
gamma=s.sqrt(7)/4; x=s.Matrix([-1,-3])/(2*s.sqrt(7))
assert s.simplify((x.T*Gm*x)[0])==1
a=[(ell*T**i*x)[0] for i in range(3)]
SS=a[0]*(x.T*Gm*T**2*x)[0]+a[1]*(x.T*Gm*T*x)[0]
v=a[0]*x
M=sum((T**(2-i)*x*ell*T**i for i in range(3)),s.zeros(2))
K=(M+Gm.inv()*M.T*Gm)/2
response=s.simplify(a[1]**2-gamma*SS+gamma*(v.T*Gm*K.inv()*v)[0])
assert response==s.Rational(1,7)
record('Three-state first response coefficient b3=1/7','exact')

# Four-state response: exact arithmetic in Q[a]/(64a^4-30a^2+3).
a=s.symbols('a'); modulus=64*a**4-30*a**2+3
def red(e):
    num,den=s.cancel(e).as_numer_denom()
    num=s.rem(num,modulus,a); den=s.rem(den,modulus,a)
    return s.rem(num*s.invert(den,modulus,a),modulus,a).expand()
def matred(M):return M.applyfunc(red)
beta=6*a-s.Rational(64,3)*a**3; gg=2-8*a*a
T=s.Matrix([[0,1,0],[0,0,1],[-a,-gg,-beta]])
Qm=s.Matrix([[342+768*a*a,1116*a-4224*a**3,528-2304*a*a],[1116*a-4224*a**3,798-576*a*a,1980*a-6272*a**3],[528-2304*a*a,1980*a-6272*a**3,858-192*a*a]])/333
ell=s.Matrix([[1,0,0]])
assert matred(Qm-T.T*Qm*T-ell.T*ell)==s.zeros(3)
x0=s.Matrix([-a**3,s.Rational(21,8)*a*a-s.Rational(9,16),s.Rational(27,4)*a**3-s.Rational(15,8)*a])
assert matred(Qm*x0-(ell*T**3).T)==s.zeros(3,1)
G=red((ell*T**3*x0)[0]); assert G==s.Rational(45,64)-s.Rational(33,32)*a*a
cc=[red((ell*T**i*x0)[0]) for i in range(4)]
S0=red(sum(cc[j]*(x0.T*Qm*T**(3-j)*x0)[0] for j in range(3)))
v0=matred(cc[0]*T*x0+cc[1]*x0)
M0=matred(sum((T**(3-j)*x0*ell*T**j for j in range(4)),s.zeros(3)))
H0=matred((Qm*M0+M0.T*Qm)/2)
expectedH=s.Matrix([[s.Rational(21,32)+s.Rational(7,16)*a*a,a/2-3*a**3,s.Rational(5,8)-s.Rational(9,4)*a*a],[a/2-3*a**3,s.Rational(7,4)-s.Rational(5,2)*a*a,3*a-s.Rational(20,3)*a**3],[s.Rational(5,8)-s.Rational(9,4)*a*a,3*a-s.Rational(20,3)*a**3,3-s.Rational(19,3)*a*a]])
assert matred(H0-expectedH)==s.zeros(3)
vv=red((v0.T*Qm*matred(H0.inv())*Qm*v0)[0])
assert vv==s.Rational(171,65536)+s.Rational(5625,32768)*a*a
b4=red((cc[2]**2-S0+vv)/G)
assert b4==s.Rational(75,496)-s.Rational(103,248)*a*a
assert s.simplify(b4.subs(a*a,(15-s.sqrt(33))/64)-(855+103*s.sqrt(33))/15872)==0
record('Four-state response: Stein, cyclic vector, response matrix, exact coefficient','exact')

# Fibonacci leaves and their invariant update, using rational endpoint arithmetic.
def leaves(n,prefix=()):
    if n<=0:return [prefix]
    return leaves(n-2,prefix+(0,))+leaves(n-1,prefix+(1,))
def interval(word):
    lo,hi=Q(0),Q(1)
    for a in reversed(word):
        lo,hi=((lo/4,hi/4) if a==0 else (lo/2+Q(1,2),hi/2+Q(1,2)))
    return lo,hi
def trim(word,n):
    cost=0
    for j,a in enumerate(word):
        cost+=2 if a==0 else 1
        if cost>=n:return word[:j+1]
    raise AssertionError('not a stopped address')
fib=[0,1]
for _ in range(20):fib.append(sum(fib[-2:]))
states=transitions=0
for n in range(1,14):
    theta=Q(1,2**n); ls=leaves(n); assert len(ls)==fib[n+2]
    ints={w:interval(w) for w in ls}; endpoints=sorted(hi for lo,hi in ints.values())
    assert min(y-x for x,y in zip(endpoints,endpoints[1:]))>=3*theta/2
    for w,(lo,hi) in ints.items():
        assert hi-lo in (theta,theta/2)
        for a in [0,1]:
            v=trim((a,)+w,n); L,U=ints[v]
            image=(lo/4,hi/4) if a==0 else (lo/2+Q(1,2),hi/2+Q(1,2))
            assert L<=image[0]<=image[1]<=U
            transitions+=1
    states+=len(ls)
record('Fibonacci optimum: exact counts, endpoint separation and all leaf updates','exact',precision_levels=13,states=states,transitions=transitions)

# Difference filtration: impulse and invertible neutral padding.
for dim in range(2,8):
    J=s.zeros(dim)
    for i in range(dim-1):J[i,i+1]=1
    U=s.eye(dim)
    for n in range(dim+1):
        assert U.rank()==max(dim-n,0)
        U=U*J
    assert s.Matrix.vstack(s.eye(dim)*J,s.eye(dim)).rank()==dim
record('Suffix-window filtration: nilpotent decrease and neutral padding','exact',dimensions=list(range(2,8)))

report={'scope':'Selected explicit consistency checks; not a complete proof audit or priority certification.',
        'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'sympy':s.__version__,
        'checks':records,'passed_groups':len(records),'elapsed_seconds':time.time()-start}
Path(__file__).with_name('verification_results.json').write_text(json.dumps(report,indent=2))
print(json.dumps({'passed_groups':len(records),'elapsed_seconds':report['elapsed_seconds']}))
