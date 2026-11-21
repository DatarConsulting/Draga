#!/usr/bin/env python3
"""F2 diagnostics, not proofs or global certificates. Requires numpy/scipy.
Run: python k2_first_correction_checks.py --output k2_first_correction_checks.json
Exact-in-principle infinite norms use finite Stein Gramians; no word truncation.
"""
import argparse,json,platform
from pathlib import Path
import numpy as np
import scipy
from scipy.linalg import block_diag,solve

def adj(x):return x.conj().T

def stein(letters,source,column=False):
 n=source.shape[0]
 K=sum((np.kron(a.T,adj(a)) if column else np.kron(a.conj(),a)) for a in letters)
 value=solve(np.eye(n*n)-K,source.reshape(-1,order='F')).reshape((n,n),order='F')
 value=(value+adj(value))/2
 update=sum((adj(a)@value@a if column else a@value@adj(a)) for a in letters)
 res=np.linalg.norm(value-source-update)/max(1.,np.linalg.norm(value))
 return value,float(res)

def psqrt(P):
 d,U=np.linalg.eigh((P+adj(P))/2)
 if d.min()<-1e-9*max(1.,np.abs(d).max()):raise ValueError('Nonpositive Gramian')
 return (U*np.sqrt(np.maximum(d,0)))@adj(U)

def case(q):
 r=np.sqrt(1-q);a=np.sqrt(1-q-2*q*q);c=q*(1+2*q)/2
 Fa=np.diag([r,0.,0.]);Fb=np.array([[0.,1.,0.],[0.,0.,1.],[0.,0.,0.]])
 Aa=np.diag([a,0.]);Ab=np.array([[q,1.],[-q,1.]])/np.sqrt(2)
 bf=np.array([[0.],[0.],[1.]]);cf=np.array([[q,0.,0.]])
 bg=np.array([[-q],[1.]]);cg=np.array([[c,0.]])
 letters=[block_diag(Fa,Aa),block_diag(Fb,Ab)]
 b=np.r_[bf,bg];out=np.c_[cf,-cg]
 R,rr=stein(letters,b@adj(b));Q,rq=stein(letters,adj(out)@out,column=True)
 root=psqrt(R);D=float(np.linalg.eigvalsh(root@Q@root).max())
 # Independently substitute the exact analytic block-Gramian formulas.
 den=1-r*a
 Rg=np.array([[(1+q*q)/(q*(1+q)),(1-q)/(1+q)],
              [(1-q)/(1+q),(2+q+2*q*q+q**3)/(1+q)]])
 Qg=q*(1+2*q)**2/(4*(1+q))*np.eye(2)
 X=np.array([[(1+q+q*q-q**3)/(2*den),(1-q+q*q+q**3)/2],
             [(1-q*q)/np.sqrt(2),(1+q*q)/np.sqrt(2)],[-q,1.]])
 Y=np.empty((3,2));Y[0]=[q*c/den,0.];Y[1]=Y[0]@Ab;Y[2]=Y[1]@Ab
 Re=np.block([[np.diag([1/q,1.,1.]),X],[X.T,Rg]])
 Qe=np.block([[q*np.eye(3),-Y],[-Y.T,Qg]])
 gramdiff=max(np.linalg.norm(R-Re)/max(1.,np.linalg.norm(R)),
              np.linalg.norm(Q-Qe)/max(1.,np.linalg.norm(Q)))
 exactroot=psqrt(Re);De=float(np.linalg.eigvalsh(exactroot@Qe@exactroot).max())
 return {'q':q,'actual_infinite_error_squared':D,
         'correction_divided_by_q':(D-.75)/q,
         'second_order_remainder_ratio':(D-.75+q/4)/q**2,
         'independent_exact_Gramian_error_squared':De,
         'Gramian_relative_difference':float(gramdiff),'Stein_residual':max(rr,rq),
         'passed':bool(gramdiff<2e-9 and abs(D-De)<2e-9)}

def derivative():
 y=1/np.sqrt(2)
 R0=np.zeros((5,5));R0[0,0]=R0[3,3]=1;R0[0,3]=R0[3,0]=.5
 R1=np.array([[0,0,0,0,.5],[0,1,0,y,y],[0,0,1,0,1],
              [0,y,0,-1,1],[.5,y,1,1,2]])
 Q0=np.array([[1,0,0,-.5,0],[0,1,0,0,-y/2],[0,0,1,0,-.25],
              [-.5,0,0,.25,0],[0,-y/2,-.25,0,.25]])
 Q1=np.array([[0,0,0,-.5,0],[0,0,0,-y/2,-y/2],[0,0,0,.25,-.5],
              [-.5,-y/2,.25,.75,0],[0,-y/2,-.5,0,.75]])
 v=np.eye(5)[:,0];w=v-.5*np.eye(5)[:,3]
 val=float(w@(R1@Q0+R0@Q1)@v)
 return {'lambda0':.75,'lambda_derivative':val,'passed':val==-.25}

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--output',default='k2_first_correction_checks.json');args=ap.parse_args()
 cases=[case(q) for q in [.1,.03,.01,.003,.001,.0003,.0001]]
 der=derivative();data={'status':'floating-point diagnostics, not proof or global optimizer certificates',
   'environment':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__},
   'summary':{'all_passed':der['passed'] and all(c['passed'] for c in cases),
              'infinite_norm_cases':len(cases),'analytic_first_correction':-.25,
              'max_gramian_relative_difference':max(c['Gramian_relative_difference'] for c in cases)},
   'derivative':der,'cases':cases}
 Path(args.output).write_text(json.dumps(data,indent=2)+'\n');print(json.dumps(data['summary'],indent=2))
 if not data['summary']['all_passed']:raise SystemExit('Diagnostic failure')
if __name__=='__main__':main()
