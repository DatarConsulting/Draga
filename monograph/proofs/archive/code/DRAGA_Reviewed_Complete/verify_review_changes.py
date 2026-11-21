"""Diagnostic checks for the revised DRAGA manuscripts.
These finite checks support transcription/algebra QA; they are not proofs of
optimality, all-parameter stability, or bibliographic novelty.
Requires Python 3, NumPy and SymPy. No external data is used.
"""
import json
from pathlib import Path
import numpy as np
import sympy as sp

def star(a): return np.asarray(a).conj().T

def trial(s,D,l,d,lam,w):
    m=D.shape[0];n=m+1;r=np.sqrt(1-s*s)
    row=l@np.linalg.matrix_power(D,d-1)
    gamma=np.linalg.norm(row);x=star(row)/gamma
    ev,U=np.linalg.eigh(np.eye(m)-s*s*(star(w)@w))
    F=(U*np.sqrt(ev))@star(U)
    N=np.vstack((s*w,D@F,l@F))
    z=np.vstack(([[s*lam]],np.sqrt(1-s*s*abs(lam)**2)*x,[[0.]]))
    W=np.eye(n+1)-(1-r)*(z@star(z))
    V=np.column_stack((s*z,W@N))
    B,l_s=V[:n,:],V[n:,:]
    v=np.vstack(([[r]],-s*star(N)@z))
    e0=np.eye(n)[:,[0]];A=e0@star(v)
    omega=s*l_s@np.linalg.matrix_power(B,d)@np.linalg.inv(np.eye(n)-r*A)
    phi=lambda X:star(A)@X@A+star(B)@X@B
    term=star(l_s)@l_s;Q=term.copy()
    for _ in range(d+1):term=phi(term);Q=Q+term
    C=np.kron(A.conj(),A)+np.kron(B.conj(),B)
    return dict(A=A,B=B,l=l_s,omega=omega,Q=Q,rho=float(max(abs(np.linalg.eigvals(C)))),energy=float(np.linalg.norm(star(A)@A+star(B)@B+star(l_s)@l_s-np.eye(n))))

def main():
    D2=sp.Matrix([[-sp.Rational(4,7),sp.sqrt(5)/7],[-5*sp.sqrt(5)/14,-sp.Rational(3,7)]])
    l2=sp.Matrix([[-sp.sqrt(7)/14,sp.sqrt(35)/7]])
    assert sp.simplify(D2.T*D2+l2.T*l2)==sp.eye(2)
    assert sp.simplify((l2*D2**2)*(l2*D2**2).T)[0]==sp.Rational(7,16)
    a=sp.Rational(1153)/(2*(1103+33*sp.sqrt(33)))
    b=(1103-33*sp.sqrt(33))/2048
    assert sp.simplify(a-b)==0
    assert sp.Rational(1,2)*(sp.Rational(9,7)+sp.Rational(25,23))==sp.Rational(191,161)<sp.Rational(9,7)
    cases=[('one_fast_state',np.array([[.5]]),np.array([[np.sqrt(.75)]]),2),('two_fast_states',np.array(D2,float),np.array(l2,float),3)]
    records=[]
    for name,D,l,d in cases:
        for complex_trial in (False,True):
            m=D.shape[0];lam=.3+(.1j if complex_trial else 0)
            w=np.array([[.2/(j+1)+(.07j if complex_trial else 0) for j in range(m)]])
            E=np.diag([-1]+[1]*m)
            for s in (.005,.03,.1):
                pos=trial(s,D,l,d,lam,w);neg=trial(-s,D,l,d,lam,w)
                parity=max(np.linalg.norm(neg['A']-E@pos['A']@E),np.linalg.norm(neg['B']-E@pos['B']@E),np.linalg.norm(neg['l']-pos['l']@E),np.linalg.norm(neg['omega']+pos['omega']@E))
                eta=float(min(np.linalg.eigvalsh(pos['Q'])))
                assert pos['energy']<1e-10 and parity<1e-10
                assert eta>0 and pos['rho']<1
                records.append(dict(case=name,complex_trial=complex_trial,s=s,energy_residual=pos['energy'],parity_residual=float(parity),observation_min_eigenvalue=eta,cp_spectral_radius=pos['rho']))
    result=dict(status='PASS',exact_checks=4,finite_trials=len(records),scope='Diagnostics only; not a replacement for the manuscript proofs.',trials=records)
    Path(__file__).with_name('review_validation_results.json').write_text(json.dumps(result,indent=2))
    print(json.dumps({k:v for k,v in result.items() if k!='trials'},indent=2))
if __name__=='__main__':main()
