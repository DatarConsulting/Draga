#!/usr/bin/env python3
"""Exact algebra and finite numerical checks, not formal proof certification."""
from pathlib import Path
from itertools import product
from fractions import Fraction
import json, math
import numpy as np
import sympy as sp
from scipy.linalg import block_diag

HERE=Path(__file__).resolve().parent
results=[]
def record(name,kind,**details):
    results.append(dict(name=name,kind=kind,passed=True,**details))
Q=sp.Rational
def spectral(M): return float(np.linalg.norm(M,2))
def gram(A,F):
    n=len(A[0]);K=sum(np.kron(T.conj(),T) for T in A)
    G=np.linalg.solve(np.eye(n*n)-K,F.reshape(-1,order='F')).reshape(n,n,order='F')
    return (G+G.conj().T)/2
def psqrt(G):
    v,U=np.linalg.eigh(G)
    assert min(v)>-1e-8
    return (U*np.sqrt(np.maximum(v,0)))@U.conj().T

# Exact minimal balanced chain, and its joint reference defects.
A0=sp.zeros(3);A0[0,1]=1
A1=sp.zeros(3);A1[1,2]=1
A=[A0,A1];t=Q(1,4);b=sp.Matrix([0,0,Q(1,2)]);c=sp.Matrix([Q(1,2),0,0])
Sig=t*sp.eye(3)
assert b*b.T+sum((T*Sig*T.T for T in A),sp.zeros(3))==Sig
assert c*c.T+sum((T.T*Sig*T for T in A),sp.zeros(3))==Sig
assert all(T*U*V==sp.zeros(3) for T in A for U in A for V in A)
sum_of_norms=b[1:,0].dot(b[1:,0])+t
coupling=(A0[:1,1:]*Sig[1:,1:]*A0[:1,1:].T)[0]
assert sum_of_norms==2*t>t and coupling==t>t*t
T=[sp.diag(U,sp.zeros(1)) for U in A]
beta=b.col_join(sp.zeros(1,1));chi=c.col_join(sp.Matrix([-Q(1,2)]))
P0=sp.diag(Sig,sp.Matrix([[t]]));P0[0,3]=P0[3,0]=t
Q0=P0.copy();Q0[0,3]=Q0[3,0]=-t
assert P0*Q0==sp.diag(0,t*t,t*t,0)
DP=beta*beta.T+sum((U*P0*U.T for U in T),sp.zeros(4))-P0
S=sum((U[:1,1:]*Sig[1:,1:]*U[:1,1:].T for U in A),sp.zeros(1))
V=sum((U[1:,1:]*Sig[1:,1:]*U[:1,1:].T for U in A),sp.zeros(2,1))
Dexpected=sp.zeros(4);Dexpected[:1,3:]=-S;Dexpected[1:3,3:]=-V
Dexpected[3:,:1]=-S;Dexpected[3:,1:3]=-V.T;Dexpected[3:,3:]=-S
assert DP==Dexpected
record('balanced_chain_and_exact_reference_defect','exact',sigma=str(t),
       norm_sum=str(sum_of_norms),coupling=str(coupling),reference_cancellation=True)

# Positive-resolvent order inequality on a nilpotent complex-valued test family.
count=0
for u in [1,2,5]:
    Z=sp.Matrix([[0,u],[0,0]])
    X=sp.eye(2)+Z*Z.T
    for M in [sp.eye(2),-sp.eye(2),sp.Matrix([[0,1],[1,0]]),sp.Matrix([[0,sp.I],[-sp.I,0]])]:
        R=M+Z*M*Z.T
        assert (X+R).is_positive_semidefinite
        assert (X-R).is_positive_semidefinite
        count+=1
record('positive_stein_resolvent_order','exact',cases=count)

# Numerical tests of full infinite Hankel error via finite joint Gramians.
rng=np.random.default_rng(2026110716)
cases=[]
for complex_case in (False,True):
  for n in (2,3,4,5):
    for trial in range(3):
      d=2+trial%2
      def draw(shape):
          x=rng.normal(size=shape)
          return x+1j*rng.normal(size=shape) if complex_case else x
      A=[.08*draw((n,n)) for _ in range(d)]
      b=draw(n);c=draw(n)
      P=gram(A,np.outer(b,b.conj()));Qm=gram([Z.conj().T for Z in A],np.outer(c,c.conj()))
      Ps=psqrt(P);v,U=np.linalg.eigh(Ps@Qm@Ps);order=np.argsort(v)[::-1]
      assert min(v)>0
      sigma=np.sqrt(v[order]);U=U[:,order]
      change=Ps@U@np.diag(1/np.sqrt(sigma));inv=np.linalg.inv(change)
      A=[inv@Z@change for Z in A];b=inv@b;c=change.conj().T@c
      Sig=np.diag(sigma)
      assert np.allclose(gram(A,np.outer(b,b.conj())),Sig,atol=2e-8)
      assert np.allclose(gram([Z.conj().T for Z in A],np.outer(c,c.conj())),Sig,atol=2e-8)
      XA=gram(A,np.eye(n));YA=gram([Z.conj().T for Z in A],np.eye(n))
      for k in range(1,n):
        B=[Z[:k,:k] for Z in A]
        assert max(abs(np.linalg.eigvals(sum(np.kron(Z.conj(),Z) for Z in B))))<1
        T=[block_diag(Z,W) for Z,W in zip(A,B)]
        beta=np.r_[b,b[:k]];chi=np.r_[c,-c[:k]]
        PE=gram(T,np.outer(beta,beta.conj()))
        QE=gram([Z.conj().T for Z in T],np.outer(chi,chi.conj()))
        P0=np.zeros((n+k,n+k),dtype=complex);Q0=P0.copy()
        for F,sgn in [(P0,1),(Q0,-1)]:
            F[:n,:n]=Sig;F[n:,n:]=Sig[:k,:k]
            F[:k,n:]=sgn*Sig[:k,:k];F[n:,:k]=sgn*Sig[:k,:k]
        DP=np.outer(beta,beta.conj())+sum(Z@P0@Z.conj().T for Z in T)-P0
        DQ=np.outer(chi,chi.conj())+sum(Z.conj().T@Q0@Z for Z in T)-Q0
        KP=max(spectral(XA),spectral(gram(B,np.eye(k))))
        KQ=max(spectral(YA),spectral(gram([Z.conj().T for Z in B],np.eye(k))))
        ep=KP*spectral(DP);eq=KQ*spectral(DQ);delta=sigma[k]
        bound2=max((delta+ep)*(delta+eq),2*sigma[0]*max(ep,eq)+ep*eq)
        error2=max(0,float(np.linalg.eigvalsh(psqrt(PE)@QE@psqrt(PE))[-1]))
        assert min(np.linalg.eigvalsh(P0+ep*np.eye(n+k)-PE))>-2e-8
        assert min(np.linalg.eigvalsh(Q0+eq*np.eye(n+k)-QE))>-2e-8
        MP=spectral(sum(Z[:k,k:]@Z[:k,k:].conj().T for Z in A))
        MQ=spectral(sum(Z[k:,:k].conj().T@Z[k:,:k] for Z in A))
        cp=KP*(2*MP+math.sqrt(MP));cq=KQ*(2*MQ+math.sqrt(MQ))
        assert spectral(DP)<=(2*MP+math.sqrt(MP))*delta+2e-8
        assert spectral(DQ)<=(2*MQ+math.sqrt(MQ))*delta+2e-8
        coarse2=max(delta**2*(1+cp)*(1+cq),2*sigma[0]*delta*max(cp,cq)+cp*cq*delta**2)
        assert error2<=bound2+2e-8 and bound2<=coarse2+2e-8
        cases.append(dict(complex=complex_case,n=n,k=k,trial=trial,
                          error=math.sqrt(error2),bound=math.sqrt(bound2),
                          explicit_bound=math.sqrt(coarse2)))
record('commuting_reference_infinite_hankel_certificate','numerical',cases=cases)

# Exact homographic counterexample to the proposed strict factor.
x=sp.symbols('x',real=True)
G=[(-Q(1,2)+x/5)/(1+x/10),(Q(1,2)+x/5)/(1+x/10)]
weights=[Q(1,4)/(1+x/10)**2,Q(3,20)/(1+x/10)**2]
assert all(sp.simplify(sp.diff(f,x)-w)==0 for f,w in zip(G,weights))
gamma=Q(25,81);c=Q(2,9);L=Q(9,28)
factor=gamma+(1-gamma)*Q(1,11)
assert L==c/(1-gamma) and factor*L==Q(331,2772)<c
assert abs(sp.diff(sp.log(sum(weights)),x).subs(x,-1))==c
assert G[0].subs(x,1)<G[1].subs(x,-1)
assert Q(9,8)<Q(3,2)<2
record('proposed_cone_factor_counterexample','exact',actual_log_lipschitz=str(c),claimed=str(factor*L))

# Finite endpoint iteration needs clipping even without sampling error.
T0=sp.Matrix([[Q(1,4),0],[Q(1,8),Q(1,8)]])
T1=sp.Matrix([[0,Q(1,4)],[0,Q(1,8)]])
stop=sp.ones(2,1)-(T0+T1)*sp.ones(2,1)
assert stop==sp.Matrix([Q(1,2),Q(5,8)])
assert sp.Matrix.vstack(sp.Matrix([[0,1]]),sp.Matrix([[0,1]])*T0).det()!=0
assert sp.Matrix.hstack(stop,T0*stop).det()!=0
for j in range(1,16):
    r=1-Q(1,2**j);next_r=(1+r)/2
    assert next_r>r
    left_weight=Q(1,4)*(r-next_r)/r
    assert left_weight<0
record('finite_hull_endpoint_negative_weight','exact',depths=15)

# Covariance counterexample uses a log-concave continuous law fully in one bin.
claimed=(1/sp.sqrt(300)-Q(1,2))**2
assert sp.simplify(claimed-Q(1,12)).is_positive
record('quantized_dither_covariance_counterexample','exact',true_variance='1/12',
       claimed_lower_bound=str(claimed))

# Packed variable-length state transitions, including every startup length.
def encode(w,d):
    value=0
    for a in w:value=d*value+a+1
    return value
def decode(I,d):
    out=[]
    while I:
        out.append((I-1)%d);I=(I-1)//d
    return tuple(reversed(out))
state_count=transition_count=0
for d in (2,3,4):
  for m in range(1,6):
    start=(d**m-1)//(d-1);size=(d**(m+1)-1)//(d-1)
    ids=set()
    for length in range(m+1):
      for w in product(range(d),repeat=length):
        I=encode(w,d);ids.add(I);state_count+=1
        assert decode(I,d)==w
        for a in range(d):
            J=d*I+a+1 if I<start else start+(d*(I-start)+a)%d**m
            assert decode(J,d)==(w+(a,))[-m:]
            transition_count+=1
    assert ids==set(range(size))
record('packed_startup_encoding_and_updates','exact',states=state_count,transitions=transition_count)

# Cantor model with x0=2/3, not a fixed point of either branch.
x0=Fraction(2,3);history_count=0
for m in range(1,7):
    for length in range(m+3):
      for w in product(range(2),repeat=length):
        actual=x0
        for a in w:actual=actual/3+Fraction(2*a,3)
        rep=x0
        for a in decode(encode(w[-m:],2),2):rep=rep/3+Fraction(2*a,3)
        assert abs(actual-rep)<=Fraction(1,3**m)
        if length<=m:assert actual==rep
        history_count+=1
record('non_fixed_initial_state_address_error','exact',initial_state='2/3',histories=history_count)

payload=dict(version='1.6',all_passed=True,groups=len(results),results=results,
             limitations='Finite exact/numerical checks complement written proofs; no formal proof-assistant or statistical coverage claim.')
(HERE/'extension_verification_results.json').write_text(json.dumps(payload,indent=2)+'\n')
print(json.dumps(dict(groups=len(results),all_passed=True,
                     joint_gramian_cases=len(cases),packed_transitions=transition_count)))
