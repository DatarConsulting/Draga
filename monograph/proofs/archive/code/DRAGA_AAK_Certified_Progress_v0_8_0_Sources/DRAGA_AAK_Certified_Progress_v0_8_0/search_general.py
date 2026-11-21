"""Numerical exploration, never a certificate of global optimality."""
import json, sys, time
from pathlib import Path
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import minimize

ROOT=Path(__file__).resolve().parent
def target(t=.01):
    A=[]
    for i in range(2):
        mat=np.zeros((6,6));j=1+2*i
        mat[0,j]=1;mat[j,5]=1;mat[j,j+1]=1;mat[j+1,5]=t
        A.append(mat)
    return A,np.eye(6)[0],np.eye(6)[5]
A,alpha,beta=target()
PA=np.outer(alpha,alpha); QA=np.outer(beta,beta)
for _ in range(3):
    PA=np.outer(alpha,alpha)+sum(a.T@PA@a for a in A)
    QA=np.outer(beta,beta)+sum(a@QA@a.T for a in A)

def matrices(x):
    return [np.array([[0.,1.],[-x[0],x[1]]]),np.array(x[2:6]).reshape(2,2)], np.array(x[6:8])

def objective(x):
    B,nu=matrices(x); gam=np.array([1.,0.])
    T=sum(np.kron(b,b) for b in B)
    rad=max(abs(np.linalg.eigvals(T)))
    if rad>=.99999 or not np.isfinite(rad): return 10+float(rad)+.01*np.linalg.norm(x)
    try:
        PB=np.linalg.solve(np.eye(4)-T.T,np.outer(gam,gam).ravel()).reshape(2,2)
        QB=np.linalg.solve(np.eye(4)-T,np.outer(nu,nu).ravel()).reshape(2,2)
        PC=-np.outer(alpha,gam); QC=np.outer(beta,nu)
        for _ in range(3):
            PC=-np.outer(alpha,gam)+sum(a.T@PC@b for a,b in zip(A,B))
            QC=np.outer(beta,nu)+sum(a@QC@b.T for a,b in zip(A,B))
        P=np.block([[PA,PC],[PC.T,PB]]); Q=np.block([[QA,QC],[QC.T,QB]])
        w,v=eigh((P+P.T)/2,check_finite=False)
        if w[0]<-1e-7: return 10+abs(w[0])
        R=(v*np.sqrt(np.maximum(w,0)))@v.T
        ev=eigh((R@Q@R+R@Q.T@R)/2,eigvals_only=True,check_finite=False)
        return float(np.sqrt(max(ev[-1],0)))
    except np.linalg.LinAlgError: return 100.

def main():
    rng=np.random.default_rng(260914)
    base=np.array([.580274,1.230292,0,0,0,0,.106299,.343126])
    best=base.copy(); val=objective(best); log=[]
    print('BASE',val,flush=True)
    for i in range(24):
        if i<12:
            x=base.copy(); x[2:6]+=rng.normal(size=4)*([.02,.1,.25][i%3])
            if i>=6: x[[0,1,6,7]]+=rng.normal(size=4)*.07
        else:
            x=best+rng.normal(size=8)*(.04 if i%2 else .15)
        res=minimize(objective,x,method='Nelder-Mead',options={'maxiter':2400,'xatol':2e-8,'fatol':2e-9})
        if res.fun<val: best=res.x.copy(); val=float(res.fun)
        item={'run':i,'error':float(res.fun),'x':res.x.tolist(),'iterations':int(res.nit)}
        log.append(item)
        (ROOT/'general_discovery.json').write_text(json.dumps({'best_error':val,'best':best.tolist(),'runs':log},indent=2))
        print(i,'error',res.fun,'best',val,'x',res.x,flush=True)

def random_search():
    rng=np.random.default_rng(18082008);log=[]
    for i in range(20):
        x=np.r_[rng.uniform(-.2,.3),rng.uniform(-.5,.5),rng.normal(size=4)*.18,rng.normal(size=2)*.3]
        res=minimize(objective,x,method='Nelder-Mead',options={'maxiter':2200,'xatol':1e-7,'fatol':1e-8})
        log.append({'run':i,'error':float(res.fun),'x':res.x.tolist()})
        (ROOT/'random_discovery.json').write_text(json.dumps(log,indent=2))
        print(i,res.fun,flush=True)

def symmetric_search():
    fun=lambda y:objective(np.array([y[0],y[1],0,1,-y[0],y[1],y[2],y[3]]))
    res=minimize(fun,[.29,.87,.103,.246],method='Nelder-Mead',options={'maxiter':4000,'fatol':1e-12,'xatol':1e-11})
    (ROOT/'symmetric_discovery.json').write_text(json.dumps({'error':float(res.fun),'parameters':res.x.tolist()},indent=2))
    print(res.fun,res.x,flush=True)

if __name__=='__main__':
    if '--random' in sys.argv:random_search()
    elif '--symmetric' in sys.argv:symmetric_search()
    else:main()
