#!/usr/bin/env python3
"""Reproducible algebraic checks for v1.4. Tests are not general proofs.

Requires Python 3, SymPy and NumPy. No downloaded data or Monte Carlo
coverage claim. All displayed counterexample matrices use exact rationals.
"""
from pathlib import Path
from itertools import product
import hashlib
import json
import math
import random
import sympy as S
import numpy as np

HERE = Path(__file__).resolve().parent
Q = S.Rational
checks = []

def record(name, **details):
    checks.append(dict(name=name, passed=True, **details))

def zero(M):
    return all(S.cancel(x) == 0 for x in M)

def positive(M):
    return all(x > 0 for x in M)

def infnorm(M):
    return max(sum(abs(x) for x in M.row(i)) for i in range(M.rows))

def l1(v):
    return sum(abs(x) for x in v)

def proper(alpha, A, tau):
    one = S.ones(tau.rows, 1)
    assert zero(alpha*one-S.ones(1,1))
    assert zero(sum(A, S.zeros(tau.rows))*one+tau-one)
    assert positive(alpha) and positive(tau)
    assert all(positive(B) for B in A)

def word_matrix(A, w):
    out = S.eye(A[0].rows)
    for a in w:
        out *= A[a]
    return out

def prob(alpha, A, tau, w):
    return (alpha*word_matrix(A,w)*tau)[0]

def similarity(alpha, A, tau, W):
    return alpha*W.inv(), [W*B*W.inv() for B in A], W*tau

def dilation(R,m):
    return (1+R*m)*S.eye(R)-m*S.ones(R)

# Noncommuting rank three, but a unique common word-Perron ray.
A = [S.Matrix([[4,1,1],[1,3,1],[1,2,4]])/30,
     S.Matrix([[2,3,1],[3,1,2],[1,2,3]])/30]
alpha = S.Matrix([[1,0,0]])
tau = S.Matrix([Q(3,5),Q(19,30),Q(17,30)])
one = S.ones(3,1)
assert zero((A[0]+A[1])*one+tau-one)
core = [(),(0,),(1,)]
H = S.Matrix([[prob(alpha,A,tau,u+v) for v in core] for u in core])
assert H.det() == -Q(11,50625000)
assert (A[0]*A[1]-A[1]*A[0])[0,1] == Q(1,450)
for k in range(1,6):
    for w in product(range(2), repeat=k):
        assert zero(one.T*word_matrix(A,w)-one.T/5**k)
record("common_perron_rank_three", determinant=str(H.det()), word_products=62)

# Invalid raw dilation; repaired inner-simplex dilation.
rare = S.Matrix([[Q(1,10000),Q(1,10000)],[Q(1,8),Q(1,8)]])
Pi = S.Matrix([[1,0],[Q(1,100),Q(99,100)]])
m = Q(1,98)
Wbad = dilation(2,m)*Pi
assert Wbad == S.Matrix([[9899,-99],[-1,9801]])/9800
badrow = (Wbad*rare*Wbad.inv()*S.ones(2,1))[0]
assert badrow == -Q(113851,49000000)
J = S.eye(2)/2+S.ones(2)/4
W = Wbad*J
assert W == S.Matrix([[14799,4801],[4899,14701]])/19600
transformed = W*rare*W.inv()
assert transformed == S.Matrix([[29484656149,30074228951],
                               [90086011249,91887363851]])/970200000000
assert positive(transformed)
record("raw_dilation_counterexample_and_repair", bad_row_mass=str(badrow))

# Missing column-sum factor, realized by actual one-letter prefixes.
R, mu, eta = 3,Q(1,5),Q(1,10)
Pi = S.Matrix([[1,0,0],[Q(1,10),Q(9,10),0],[0,0,1]])
J = (1-R*mu)*S.eye(R)+mu*S.ones(R)
Z = Pi*J
Qmat = S.Matrix([[Q(1,4),Q(1,1000),Q(1,1000)],
                 [Q(1,10000),Q(1,4),Q(1,1000)],
                 [Q(1,10000),Q(1,1000),Q(1,4)]])
A = [Qmat[:,a]*Z[a,:] for a in range(R)]
alpha = S.ones(1,R)/R
tau = S.ones(R,1)-sum(A,S.zeros(R))*S.ones(R,1)
assert tau == S.Matrix([Q(187,250),Q(7489,10000),Q(7489,10000)])
proper(alpha,A,tau)
for a in range(R):
    assert alpha*A[a]/(alpha*A[a]*S.ones(R,1))[0] == Z[a,:]
    rowsum=A[a]*S.ones(R,1)
    assert all(A[a][i,j] >= mu*rowsum[i] for i in range(R) for j in range(R))
H = S.Matrix([[prob(alpha,A,tau,(i,j)) for j in range(R)] for i in range(R)])
assert H.det() != 0
W = dilation(R,Q(49,100))*Z
assert W == S.Matrix([[2431,59,10],[208,2282,10],[-39,59,2480]])/2500
badrow=(W*A[0]*W.inv()*S.ones(R,1))[2]
assert badrow == -Q(94961,25000000)
proper(*similarity(alpha,A,tau,dilation(R,Q(1,4))*Z))
record("column_sum_counterexample_with_actual_prefixes", determinant=str(H.det()),
       bad_row_mass=str(badrow), corrected_upper_endpoint="5/11")

# Exact interior two-state example, independent of the preceding source.
A = [S.Matrix([[3,2],[3,3]])/20,S.Matrix([[2,2],[2,3]])/20]
alpha = S.Matrix([[Q(1,10),Q(9,10)]])
tau = S.Matrix([Q(11,20),Q(9,20)])
U = S.Matrix([[Q(30,59),Q(29,59)],[Q(20,49),Q(29,49)]])
V = S.Matrix([[2,-1],[-4,5]])*U
assert V == S.Matrix([[1760,1131],[20,2871]])/2891
av,B,tv = similarity(alpha,A,tau,V)
assert av == S.Matrix([[Q(897,5800),Q(4903,5800)]])
assert tv == S.Matrix([29539,26059])/57820
assert B[0] == S.Matrix([[24761923,2357717],[24727123,5420117]])/100606800
assert B[1] == S.Matrix([[16461862,5627438],[16427062,8689838]])/100606800
proper(av,B,tv)
assert min(list(av)+list(B[0])+list(B[1])) >= Q(1,50)
assert min(tv) >= Q(9,20)
for k in range(7):
    for w in product(range(2),repeat=k):
        assert prob(alpha,A,tau,w) == prob(av,B,tv,w)
record("interior_rank_two_exact_coefficients", words=127, gamma="1/50", stop="9/20")

# Exact sandwich and strict-margin checks for many stochastic matrices.
rng = random.Random(20261107)
grid_cases=0
for R in (2,3,4):
    for eta in (Q(1,1000),Q(1,100),Q(3,100)):
        for repetition in range(4):
            raw=S.Matrix([[rng.randint(0,10) for _ in range(R)] for _ in range(R)])
            K=S.Matrix([list(raw.row(i)/sum(raw.row(i))) for i in range(R)])
            Pi=(1-eta)*S.eye(R)+eta*K
            assert all(x >= -eta/(1-2*eta) for x in Pi.inv())
            mu=Q(1,3*R); amu=1-R*mu
            C=1+(R-1)*eta
            ml=eta/(1-2*eta); mh=mu/(amu*C); m=(ml+mh)/2
            assert ml < mh
            delta=(mh-ml)/2
            J=amu*S.eye(R)+mu*S.ones(R)
            W=dilation(R,m)*Pi*J
            Wtinv=(dilation(R,m)*Pi).inv()
            assert positive(W) and positive(Wtinv)
            assert zero(W*S.ones(R,1)-S.ones(R,1))
            assert zero(Wtinv*W-J)
            assert min(Wtinv) >= delta/(1+R*m)
            # Native transitions are positive mixtures of the inner simplex.
            A=[S.diag(*[Q(i+1,20*R) for i in range(R)])*J,
               S.diag(*[Q(R-i,20*R) for i in range(R)])*J]
            alpha=S.ones(1,R)/R
            tau=S.ones(R,1)-sum(A,S.zeros(R))*S.ones(R,1)
            av,B,tv=similarity(alpha,A,tau,W)
            proper(av,B,tv)
            wmin=min(x for M in A for x in M*S.ones(R,1))
            margin=delta/(1+R*m)*min(1,wmin)
            assert min(list(av)+[x for M in B for x in M]) >= margin
            grid_cases+=1
record("inverse_and_corrected_sandwich_grid", exact_cases=grid_cases)

# Normalization lemma and full wordwise telescoping on finite words.
norm_cases=0
for R in (2,3,4):
    alpha=S.ones(1,R)/R; gamma=Q(1,2*R); delta=gamma/(2*(R+1))
    for signs in product((-1,1), repeat=R):
        a0=alpha+S.Matrix([[delta*z for z in signs]])
        abar=a0/sum(a0)
        assert max(abs(x) for x in abar-alpha) <= (R-1)*delta/(1-R*delta) <= gamma/2
        norm_cases+=1
record("normalized_initial_error", exact_cases=norm_cases)

# Deterministic coherent perturbations of observable queries. These are
# probability tables, not a claim to have run the statistical experiment.
A=[S.Matrix([[3,2],[3,3]])/20,S.Matrix([[2,2],[2,3]])/20]
alpha=S.Matrix([[Q(1,10),Q(9,10)]])
tau=S.Matrix([Q(11,20),Q(9,20)])
R=d=2; one=S.ones(R,1); core=[(),(0,)]
P=S.Matrix.vstack(*[alpha*word_matrix(A,u) for u in core])
q=P*one; Dq=S.diag(*q); V=(Dq.inv()*P).inv()
H=S.Matrix([[prob(alpha,A,tau,u+v) for v in core] for u in core])
M=[S.Matrix([[prob(alpha,A,tau,u+(a,)+v) for v in core] for u in core]) for a in range(d)]
h=S.Matrix([[prob(alpha,A,tau,v) for v in core]])
h0=infnorm(H.inv()); qs=min(q)
vminus=infnorm(V.inv()); vplus=infnorm(V)
Ka=vminus*((2*R+1)*h0+2*R*h0*h0)
K0=2*R*h0/qs**2+(4*R*R*h0*h0+4*R*h0)/qs
numeric_rows=[]
for eps in (Q(1,10**6),Q(1,10**8)):
    # p_hat=(1-eps)*p+eps*delta_empty; all queried probabilities exact.
    Hh=(1-eps)*H; Hh[0,0]+=eps
    Mh=[(1-eps)*B for B in M]
    hh=(1-eps)*h; hh[0,0]+=eps
    qh=(1-eps)*q; qh[0]+=eps
    Dh=S.diag(*qh)
    t=max([abs(x) for x in Hh-H]+[abs(x) for x in hh-h]+[abs(x) for x in qh-q])
    assert R*h0*t <= Q(1,2) and t <= qs/2
    a0=hh*Hh.inv()*Dh*V.inv(); abar=a0/sum(a0)
    B=[V*Dh.inv()*Mb*Hh.inv()*Dh*V.inv() for Mb in Mh]
    tb=one-sum(B,S.zeros(R))*one
    proper(abar,B,tb)
    ea=Ka*t
    em=sum(infnorm(X-Y) for X,Y in zip(Mh,M))
    eA=vplus*vminus*(2*h0*em/qs+K0*t)
    assert l1(a0-alpha) <= ea
    assert sum(infnorm(X-Y) for X,Y in zip(B,A)) <= eA
    for b in (8,12,16):
        hb=Q(1,2**b)
        Ar=[X.applyfunc(lambda x:S.floor(x/hb)*hb) for X in B]
        ar=abar.applyfunc(lambda x:S.floor(x/hb)*hb)
        ar[-1]=1-sum(ar[:-1])
        tr=one-sum(Ar,S.zeros(R))*one
        proper(ar,Ar,tr)
        assert all((x/hb).q == 1 for X in [ar,tr]+Ar for x in X)
        assert l1(ar-abar)/2 <= (R-1)*hb
        assert sum(infnorm(X-Y) for X,Y in zip(Ar,B)) <= d*R*hb
        cert=ea/(1-ea)+(R-1)*hb+(eA+d*R*hb)/min(tau)
        numeric_rows.append({"mixture_epsilon":str(eps),"bits":b,"tv_bound":float(cert)})
        for k in range(7):
            for w in product(range(d),repeat=k):
                lhs=prob(ar,Ar,tr,w)-prob(alpha,A,tau,w)
                rhs=((ar-alpha)*word_matrix(Ar,w)*tr+alpha*word_matrix(A,w)*(tr-tau))[0]
                for i,a in enumerate(w):
                    rhs+=(alpha*word_matrix(A,w[:i])*(Ar[a]-A[a])*word_matrix(Ar,w[i+1:])*tr)[0]
                assert S.cancel(lhs-rhs)==0
record("coherent_query_perturbation_dyadic_output_and_telescoping",
       query_tables=2, dyadic_models=6, words_per_model=127, results=numeric_rows)

# Classical SPA application, with near-pure displacement and all thresholds.
G=np.array([[1.2,-.2],[-.3,1.3]])
H=np.array([[1-1e-7,1e-7],[1e-7,1-1e-7],[.3,.7],[.7,.3],[.5,.5]])
F=H@G
noise=np.array([[1,-1],[-1,1],[1,1],[-1,-1],[1,-1]],dtype=float)*1e-9
Fhat=F+noise
sigma=np.linalg.svd(G,compute_uv=False)[-1]; kappa=np.linalg.norm(G,2)/sigma
eF=max(np.linalg.norm(x) for x in noise); Cspa=1+80*kappa*kappa
e=eF+math.sqrt(2)*1e-7*np.linalg.norm(G,2)
assert e < sigma*min(1/(2*math.sqrt(1)),.25)/Cspa
residual=Fhat.T.copy(); selected=[]
for j in range(2):
    i=int(np.argmax(np.sum(residual**2,axis=0))); selected.append(i)
    u=residual[:,i]; residual=(np.eye(2)-np.outer(u,u)/(u@u))@residual
assert set(selected)=={0,1}
purity=(eF+Cspa*e)/sigma
assert min(H[0,0],H[1,1]) >= 1-purity
record("relative_spa_corrected_budget", selected=selected, purity_bound=purity)

claimed=101*math.log(1.5)
actual=3*math.log(103)+100*math.log(103/100)+math.log(.75)
assert claimed>40.95 and actual<16.58
record("overlapping_whole_future_kl_counterexample", false_lower_bound=claimed, actual_kl=actual)

files=[HERE/"positive_discovery.inc",HERE/"fixed_basis_sampling.inc",HERE/"discovery_audit.inc",Path(__file__)]
result={"version":"1.4","status":"all executed checks passed",
        "scope":"Exact algebra and deterministic finite checks; not a proof assistant or empirical coverage study.",
        "checks":checks,"sha256":{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in files}}
out=HERE/"discovery_verification_results.json"
out.write_text(json.dumps(result,indent=2)+"\n")
print(json.dumps({"status":result["status"],"check_groups":len(checks),"result":str(out)},indent=2))
