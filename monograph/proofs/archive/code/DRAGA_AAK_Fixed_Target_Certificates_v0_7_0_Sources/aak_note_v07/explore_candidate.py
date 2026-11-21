"""Numerical discovery only; use verify_fixed_target.py for certification."""
import json
from pathlib import Path
import numpy as np
from scipy.linalg import block_diag, solve, eigh
from scipy.optimize import minimize


def target(t=0.01):
    # Shared initial/final states, two two-intermediate-state chains.
    A=[]
    for i in range(2):
        M=np.zeros((6,6)); j=1+2*i
        M[0,j]=1; M[j,5]=1; M[j,j+1]=1; M[j+1,5]=t
        A.append(M)
    return A,np.eye(6)[0],np.eye(6)[5]


AA,alpha,beta=target()


def objective(x):
    q,p,c0,c1=x
    B=np.array([[0.,1.],[-q,p]])
    radius=max(abs(np.linalg.eigvals(B)))
    if radius >= .998:
        return 10+radius+np.linalg.norm(x)
    D=[block_diag(AA[0],B),block_diag(AA[1],np.zeros((2,2)))]
    c=np.r_[alpha,-1.,0.]; v=np.r_[beta,c0,c1]
    K=sum(np.kron(M,M) for M in D)
    P=solve(np.eye(64)-K.T,np.outer(c,c).reshape(-1),check_finite=False).reshape(8,8)
    Q=solve(np.eye(64)-K,np.outer(v,v).reshape(-1),check_finite=False).reshape(8,8)
    val,V=eigh((P+P.T)/2,check_finite=False)
    R=(V*np.sqrt(np.maximum(val,0)))@V.T
    err=np.sqrt(max(eigh((R@Q@R+R@Q.T@R)/2,eigvals_only=True,check_finite=False)[-1],0))
    return float(err)


def main():
    x0=np.array([.58127,1.2314,.10886,.34092])
    results=[{'kind':'seed','parameters':x0.tolist(),'uncertified_error':objective(x0)}]
    best=x0
    for seed in [x0,x0+np.array([.01,.01,0.,.01])]:
        res=minimize(objective,seed,method='Nelder-Mead',options={'maxiter':1000,'xatol':1e-10,'fatol':1e-10})
        results.append({'kind':'Nelder-Mead','parameters':res.x.tolist(),
                        'uncertified_error':float(res.fun),'iterations':int(res.nit)})
        if objective(res.x) < objective(best): best=res.x
    Path('aak_note_v07/discovery.json').write_text(json.dumps(results,indent=2)+'\n')
    print(json.dumps(results,indent=2),flush=True)


if __name__=='__main__': main()
