"""Uncertified section exploration."""
import sys,json,itertools
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT));from search_general import matrices

def W(N):return [w for j in range(N+1) for w in itertools.product(range(2),repeat=j)]
def coeff(w):return (1. if len(w)==2 else .01 if len(w)==3 else 0.) if w and len(set(w))==1 else 0.
def hankel(x,N):
 B,nu=matrices(x); words=W(N); Pw=[]
 for w in words:
  T=np.eye(2)
  for a in w:T=T@B[a]
  Pw.append(T)
 O=np.array([p[0] for p in Pw]);C=np.array([p@nu for p in Pw]).T
 return O@C

def main():
 x=np.array(json.loads((ROOT/'general_discovery.json').read_text())['best'])
 log=[]
 for N in [2,3,4]:
  words=W(N);F=np.array([[coeff(u+v) for v in words] for u in words])
  fun=lambda y:np.linalg.norm(F-hankel(y,N),2)
  res=minimize(fun,x,method='Nelder-Mead',options={'maxiter':5000,'fatol':1e-10,'xatol':1e-8})
  item={'N':N,'error':float(res.fun),'x':res.x.tolist()};log.append(item)
  print(item,flush=True)
 # Compressed section, sqrt(2) coordinates.
 rt=np.sqrt(2);F=np.array([[0,0,0,rt],[0,1,0,.01/rt],[0,0,1,.01/rt],[rt,.01/rt,.01/rt,0]])
 J=np.zeros((7,4));J[0,0]=J[1,1]=J[2,2]=1;J[3,3]=J[6,3]=1/rt
 G=J.T@hankel(x,2)@J
 u,s,v=np.linalg.svd(G);init=np.r_[(u[:,:2]*np.sqrt(s[:2])).ravel(),(v[:2].T*np.sqrt(s[:2])).ravel()]
 def mat(y):return y[:8].reshape(4,2)@y[8:].reshape(4,2).T
 def con(y):
  g=mat(y)
  return np.array([g[0,1]-g[1,0],g[0,2]-g[2,0],g[0,3]-g[3,0],rt*g[0,3]-g[1,1]-g[2,2]])
 fun=lambda y:np.linalg.norm(F-mat(y),2)
 rng=np.random.default_rng(808)
 for i in range(8):
  seed=init if i==0 else init+rng.normal(size=16)*.2
  res=minimize(fun,seed,method='SLSQP',constraints=[{'type':'eq','fun':con}],options={'maxiter':1400,'ftol':1e-10})
  item={'compressed':True,'error':float(res.fun),'coherence_residual':float(np.linalg.norm(con(res.x))),'G':mat(res.x).tolist()}
  log.append(item);print(item,flush=True)
 (ROOT/'section_discovery.json').write_text(json.dumps(log,indent=2))

if __name__=='__main__':main()
