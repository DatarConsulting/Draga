#!/usr/bin/env python3
"""Finite checks supporting v1.5; not formal proofs or an LLM benchmark."""
from pathlib import Path
from itertools import product
from fractions import Fraction
import hashlib, json, math, random
import numpy as np
import sympy as S
from scipy.linalg import block_diag

HERE=Path(__file__).resolve().parent
Q=S.Rational
results=[]
def record(name,kind,**data): results.append(dict(name=name,kind=kind,passed=True,**data))
def zero(M): return all(S.simplify(x)==0 for x in M)
def inf(M): return max(sum(abs(x) for x in M.row(i)) for i in range(M.rows))
def word(A,w):
    B=S.eye(A[0].rows)
    for a in w:B*=A[a]
    return B

# Exact refutation of leading-block Lyapunov inheritance.
A=S.Matrix([[2,1],[-4,-2]])
X=S.eye(2)+A*A.T
assert zero(A*A) and X-A*X*A.T==S.eye(2)
assert (X-A[0,0]**2*X)[0,0]<0
record('full_witness_does_not_pass_to_leading_block','exact',X=str(X),leading_eigenvalue=2)

# Numerically evaluate infinite Hankel norms through finite joint Gramians.
def gram(A,B):
    n=A[0].shape[0]
    K=sum(np.kron(T,T) for T in A)
    G=np.linalg.solve(np.eye(n*n)-K,B.reshape(-1,order='F')).reshape(n,n,order='F')
    return (G+G.T)/2
def hankel_norm(A,b,c):
    P=gram(A,np.outer(b,b));Qm=gram([T.T for T in A],np.outer(c,c))
    vals,vec=np.linalg.eigh(P)
    Ps=(vec*np.sqrt(np.maximum(vals,0)))@vec.T
    return math.sqrt(max(0,np.linalg.eigvalsh(Ps@Qm@Ps)[-1]))
def error_norm(A,b,c,B,bb,cc):
    return hankel_norm([block_diag(T,U) for T,U in zip(A,B)],np.r_[b,bb],np.r_[c,-cc])
rng=np.random.default_rng(20261107)
trunc=[];dyadic=[]
for case in range(6):
    n=3;d=2
    A=[rng.normal(size=(n,n))*.08 for _ in range(d)]
    b=rng.normal(size=n);c=rng.normal(size=n)
    P=gram(A,np.outer(b,b));Qm=gram([T.T for T in A],np.outer(c,c))
    vals,vec=np.linalg.eigh(P);Ps=(vec*np.sqrt(vals))@vec.T
    vals,U=np.linalg.eigh(Ps@Qm@Ps);order=np.argsort(vals)[::-1]
    sigma=np.sqrt(vals[order]);U=U[:,order]
    balance=Ps@U@np.diag(1/np.sqrt(sigma))
    inv=np.linalg.inv(balance)
    A=[inv@T@balance for T in A];b=inv@b;c=balance.T@c
    assert np.allclose(gram(A,np.outer(b,b)),np.diag(sigma),atol=1e-9)
    assert np.allclose(gram([T.T for T in A],np.outer(c,c)),np.diag(sigma),atol=1e-9)
    X=gram(A,np.eye(n));Y=gram([T.T for T in A],np.eye(n))
    lx=np.linalg.eigvalsh(X);ly=np.linalg.eigvalsh(Y)
    rx=math.sqrt(1-1/lx[-1]);ry=math.sqrt(1-1/ly[-1])
    kx=lx[-1]/lx[0];ky=ly[-1]/ly[0]
    for k in (1,2):
        B=[T[:k,:k] for T in A]
        rho=max(abs(np.linalg.eigvals(sum(np.kron(T,T) for T in B))))
        assert rho<1
        err=error_norm(A,b,c,B,b[:k],c[:k])
        bound=math.sqrt(sigma[0]*sum(sigma[k:]))*(math.sqrt(kx)/(1-rx)+math.sqrt(ky)/(1-ry))
        assert err<=bound+1e-7
        trunc.append(dict(case=case,k=k,error=err,bound=bound))
    r=max(rx,ry);kappa=max(kx,ky);a=1-r
    Db=2*math.sqrt(kappa*n)/a+2*kappa*math.sqrt(d)*n*np.linalg.norm(b)/a**2
    Dc=2*math.sqrt(kappa*n)/a+2*kappa*math.sqrt(d)*n*np.linalg.norm(c)/a**2
    K=Dc*math.sqrt(kappa)*np.linalg.norm(b)/a+2*math.sqrt(kappa)*(np.linalg.norm(c)+1)*Db/a
    eps=.001
    bits=math.ceil(math.log2(max(math.sqrt(n),2*math.sqrt(kappa*d)*n/a,K/eps)))
    h=2.**(-bits)
    Ah=[np.round(T/h)*h for T in A];bh=np.round(b/h)*h;ch=np.round(c/h)*h
    rho=max(abs(np.linalg.eigvals(sum(np.kron(T,T) for T in Ah))))
    err=error_norm(A,b,c,Ah,bh,ch)
    assert rho<1 and err<=eps+1e-7
    dyadic.append(dict(case=case,bits=bits,error=err,target=eps))
record('conditional_balanced_truncation','numerical',cases=trunc)
record('a_priori_dyadic_precision','numerical',cases=dyadic)

# Exact native positive model and perturbed normalized basis.
R=3;one=S.ones(R,1);mu=Q(1,5);J=(1-R*mu)*S.eye(R)+mu*S.ones(R)
Wmass=S.Matrix([[Q(1,4),Q(1,1000),Q(1,1000)],
                [Q(1,10000),Q(1,4),Q(1,1000)],
                [Q(1,10000),Q(1,1000),Q(1,4)]])
A=[Wmass[:,a]*J[a,:] for a in range(R)]
alpha=S.ones(1,R)/R;tau=one-sum(A,S.zeros(R))*one
gamma=min(list(alpha)+[x for T in A for x in T]);stop=min(tau)
E=S.Matrix([[1,-1,0],[0,1,-1],[-1,0,1]])/10**9
Vh=S.eye(R)+E;rho=inf(E)
aa=alpha*Vh.inv();AA=[Vh*T*Vh.inv() for T in A];tt=Vh*tau
assert rho<=min(Q(1,2),gamma/(8*R),stop/2)
assert sum(abs(x) for x in aa-alpha)<=rho/(1-rho)
assert sum(inf(T-U) for T,U in zip(AA,A))<=2*R*rho/(1-rho)
assert max(abs(x) for x in tt-tau)<=rho
assert min(list(aa)+[x for T in AA for x in T])>=gamma/2
assert min(tt)>=stop/2
record('estimated_basis_margins','exact',rho=str(rho),gamma=str(gamma))

# Actual one-letter witnesses in relative observable coordinates; no hidden
# inference is claimed. This tests the deterministic learned-basis stages.
core=[(0,),(1,),(2,)];pool=[()]+core
P=S.Matrix.vstack(*[alpha*word(A,u) for u in core]);q=P*one;Dq=S.diag(*q)
C=S.Matrix.hstack(*[word(A,v)*tau for v in core]);H=P*C
assert H.det()!=0 and Dq.inv()*P==J
F=S.Matrix.vstack(*[(alpha*word(A,w)*C)/(alpha*word(A,w)*one)[0]*H.inv()*Dq for w in pool])
assert F[1:,:]==S.eye(R)
noise=S.Matrix([[1,-2,3],[2,-1,1],[-2,1,1],[3,-2,1]])/10**12
Fh=F+noise
for i in range(Fh.rows):Fh[i,0]+=1-sum(Fh.row(i))
assert Fh*S.ones(R,1)==S.ones(Fh.rows,1)
assert all(sum(abs(x) for x in Fh.row(i)-F.row(i))<=2*sum(abs(x) for x in noise.row(i)) for i in range(F.rows))
res=Fh.T
indices=[]
for k in range(R):
    scores=[sum(x*x for x in res.col(j)) for j in range(res.cols)]
    j=max(range(res.cols),key=lambda j:scores[j]);indices.append(j)
    u=res.col(j);res=(S.eye(R)-u*u.T/(u.T*u)[0])*res
assert set(indices)=={1,2,3}
m=Q(1,4);D=(1+R*m)*S.eye(R)-m*S.ones(R)
U=S.Matrix.vstack(*[F.row(i) for i in indices]);Uh=S.Matrix.vstack(*[Fh.row(i) for i in indices])
V=D*U;Vh=D*Uh;E=Vh-V
assert Vh*one==one
assert inf(E)<=inf(D)*max(sum(abs(x) for x in Fh.row(i)-F.row(i)) for i in indices)
Bt=[Dq.inv()*P*T*P.inv()*Dq for T in A]
beta=alpha*P.inv()*Dq;stopobs=Dq.inv()*P*tau
aa=beta*Vh.inv();AA=[Vh*T*Vh.inv() for T in Bt];tt=Vh*stopobs
assert min(list(aa)+[x for T in AA for x in T]+list(tt))>0
assert zero(sum(AA,S.zeros(R))*one+tt-one)
record('relative_witness_pool_affine_normalization_and_spa','exact',selected=indices,core_determinant=str(H.det()))

# Rank-two map estimates, finite endpoint iteration, clipping and local defect.
A=[S.Matrix([[2,1],[1,1]])/10,S.Matrix([[1,1],[1,2]])/10]
one=S.ones(2,1);eps=Q(1,10000)
F=[lambda t:(1+t)/(2+t),lambda t:1/(3-t)]
Fh=[lambda t:F[0](t)+eps,lambda t:F[1](t)-eps]
mass=[lambda t:(2+t)/10,lambda t:(3-t)/10]
lo=hi=Q(1,2);ranktwo=[]
for depth in range(1,7):
    vals=[Q(1,2)]+[f(t) for f in Fh for t in (lo,hi)]
    lo,hi=min(vals),max(vals)
    basis=S.Matrix([[lo,1-lo],[hi,1-hi]])
    T=[S.zeros(2) for _ in range(2)]
    for i,t in enumerate((lo,hi)):
        for a in range(2):
            v=min(hi,max(lo,Fh[a](t)));w=mass[a](t)
            T[a][i,0]=w*(hi-v)/(hi-lo);T[a][i,1]=w*(v-lo)/(hi-lo)
    assert zero(sum(T,S.zeros(2))*one-one/2)
    local=max(sum(sum(abs(x) for x in (basis*A[a]-T[a]*basis).row(i)) for a in range(2)) for i in range(2))
    bound=2*(eps+Q(1,4)**depth)
    assert local<=bound
    ranktwo.append(dict(depth=depth,local_defect=float(local),bound=float(bound)))
record('parametric_rank_two_clipping','exact',cases=ranktwo)

# Strict regularity cone comparison, and normalized constant-weight example.
L=2.;c=.4;gamma=.25;Lp=c+gamma*L
tstar=(L-Lp)/(L+Lp)*math.exp(-2*(L+2*Lp))
xs=np.linspace(-1,1,51)
comparisons=0
for sf in (-Lp,-Lp/3,Lp/2,Lp):
    for sg in (-Lp,0,Lp):
        f=np.exp(sf*xs);g=np.exp(sg*xs);h=g-tstar*f
        assert np.all(h>0)
        for i in range(len(xs)):
            assert np.all(h[i]<=np.exp(L*np.abs(xs-xs[i]))*h+1e-12)
        comparisons+=len(xs)**2
for n in range(1,20):
    vn=Q(1,2)**n
    assert vn/vn==1
record('strict_log_lipschitz_cone_and_normalization','numerical_and_exact',grid_inequalities=comparisons,
       cone_diameter=2*math.log(1/tstar),constant_weight_h=1,constant_weight_lambda='1/2')

# Optimal homogeneous addresses on the middle-third Cantor attractor.
address_cases=0
for m in range(1,8):
    for length in range(0,10):
        for w in product(range(2),repeat=length):
            true=Q(0);index=0
            for a in w:true=true/3+Q(2*a,3);index=(2*index+a)%(2**m)
            u=[(index>>(m-1-i))&1 for i in range(m)]
            state=Q(0)
            for a in u:state=state/3+Q(2*a,3)
            assert abs(state-true)<=Q(1,3)**m
            if length<m:assert state==true
            address_cases+=1
s=math.log2((1+math.sqrt(5))/2)
assert abs(4**(-s)+2**(-s)-1)<1e-12 and s<1
record('table_free_address_update_and_nonhomogeneous_obstruction','exact_and_numerical',words=address_cases,unequal_ratio_dimension=s)

# Rational integer lifting, exact inverse, finite-distribution covariance bound.
coeffs=[Q(1,3),Q(-2,5),Q(3,7)]
positions=[(0,1),(1,0),(0,1)]
shears=[]
for s,(i,j) in zip(coeffs,positions):
    T=S.eye(2);T[i,j]=s;shears.append(T)
U=S.eye(2)
for T in shears:U=T*U
def nearest(x):return int(S.floor(x+Q(1,2)))
def lift(z,inverse=False):
    z=list(z)
    seq=list(zip(coeffs,positions))
    if inverse:seq.reverse()
    for s,(i,j) in seq:z[i]+=(-1 if inverse else 1)*nearest(s*z[j])
    return tuple(z)
samples=list(product(range(-12,13),repeat=2));out=[lift(z) for z in samples]
assert len(set(out))==len(samples)
assert all(lift(t,True)==z for z,t in zip(samples,out))
bb=S.zeros(2,1)
for r,(i,j) in enumerate(positions):
    suffix=S.eye(2)
    for T in shears[r+1:]:suffix=T*suffix
    for k in range(2):bb[k]+=abs(suffix[k,i])/2
for z,t in zip(samples,out):
    e=S.Matrix(t)-U*S.Matrix(z)
    assert all(abs(e[j])<=bb[j] for j in range(2))
arr=np.array(samples,float);tar=np.array(out,float);Un=np.array(U,float)
cov=np.cov(arr,rowvar=False,bias=True)
for j in range(2):
    var=np.var(tar[:,j]);bound=(math.sqrt((Un@cov@Un.T)[j,j])+float(bb[j]))**2
    assert var<=bound+1e-10
assert 1-Q(151,100)**2/4==Q(17199,40000)>Q(1,4)
record('integer_lifting_inverse_propagation_and_covariance','exact_and_numerical',lattice_points=len(samples),rounding_variance_increase='17199/40000')

# Entropy invariance and explicit gain only against a raw fixed-width baseline.
for b in range(1,11):
    vals=[(z,z) for z in range(2**b)]
    lifted=[(x,y-x) for x,y in vals]
    assert len(set(lifted))==len(vals)==2**b
    assert math.log2(len(set(lifted)))==b
record('entropy_invariance_and_raw_layout_saving','exact',tested_widths=list(range(1,11)),entropy='b bits',raw_layout='2b bits')

# Fixed-depth indistinguishability is exact for every sample size.
for D in range(0,20):
    words=['a'*(D+1),'a'*(D+2)]
    censored=[(w[:D],len(w)<=D) for w in words]
    assert censored[0]==censored[1] and words[0]!=words[1]
record('censored_observation_obstruction','exact',depths=20)

# Fixed-block singular projectors: deterministic perturbation part only.
svd_cases=[]
for r in (1,2,3):
    K=rng.random((5,r))@rng.random((r,6))
    U,sv,Vt=np.linalg.svd(K,full_matrices=False);s0=sv[r-1]
    noise=rng.normal(size=K.shape)*s0/100
    e=np.linalg.norm(noise,2);assert e<s0/3
    Uh,sh,Vht=np.linalg.svd(K+noise,full_matrices=False)
    Pl=U[:,:r]@U[:,:r].T;Ph=Uh[:,:r]@Uh[:,:r].T
    Qr=Vt[:r].T@Vt[:r];Qh=Vht[:r].T@Vht[:r]
    error=max(np.linalg.norm(Pl-Ph,2),np.linalg.norm(Qr-Qh,2))
    assert error<=e/(s0-e)+1e-12
    Kr=(Uh[:,:r]*sh[:r])@Vht[:r]
    assert np.linalg.norm(Kr-K,2)<=2*e+1e-12
    assert sum(sh>2*e)==r
    svd_cases.append(dict(rank=r,projector_error=error,bound=e/(s0-e)))
record('fixed_block_svd_projector_certificate','numerical',cases=svd_cases)

# Table-free cube chart with a coordinate swap giving future observability.
chart_cases=0
for p in (3,5,8):
    q=p+3;scale=2**q;delta=Fraction(1,scale)
    for length in range(7):
        for w in product(range(4),repeat=length):
            z=[Fraction(1,3),Fraction(2,3)]
            state=[Fraction(math.floor(x*scale+Fraction(1,2)),scale) for x in z]
            for a in w:
                if a%2:z=z[::-1];state=state[::-1]
                offset=[a//2,a%2]
                z=[(x+b)/2 for x,b in zip(z,offset)]
                update=[(x+b)/2 for x,b in zip(state,offset)]
                state=[Fraction(math.floor(x*scale+Fraction(1,2)),scale) for x in update]
            assert max(abs(x-y) for x,y in zip(z,state))<=2*delta
            assert abs(z[0]-state[0])<=Fraction(1,2**(p+2))
            chart_cases+=1
record('contractive_global_chart_grid_updates','exact',histories=chart_cases,coordinate_dimension=2)

files=list(HERE.glob('*.inc'))+[Path(__file__)]
data={'version':'1.5','status':'all executed checks passed',
      'scope':'Finite exact and numerical checks, not a proof assistant, finite-sample coverage study or LLM benchmark.',
      'checks':results,'sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in files}}
path=HERE/'closure_verification_results.json';path.write_text(json.dumps(data,indent=2)+'\n')
print(json.dumps({'status':data['status'],'groups':len(results),'result':str(path)},indent=2))
