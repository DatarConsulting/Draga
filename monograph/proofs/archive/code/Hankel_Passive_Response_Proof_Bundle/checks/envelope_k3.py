import numpy as np
from scipy.linalg import expm, solve_discrete_lyapunov, sqrtm
from scipy.optimize import minimize
import json

T=np.array([[0.,1.],[-.5,-1.]])
l=np.array([[1.,0.]])
G=solve_discrete_lyapunov(T.T,l.T@l)
Gh=sqrtm(G);Gi=np.linalg.inv(Gh)
D=Gh@T@Gi; ell=l@Gi
v=(ell@np.linalg.matrix_power(D,2)).ravel()
N0=np.zeros((4,2));N0[1:3,:]=D;N0[3,:]=ell
z0=np.zeros(4);z0[1:3]=v/np.linalg.norm(v)

# z is normalized z0+delta, N=e^K N0, N^*N=I.
def parameters(x,complex_case=False):
    if complex_case:
        dz=x[:4]+1j*x[4:8]; z=z0+dz; z=z/np.linalg.norm(z)
        y=x[8:];H=np.zeros((4,4),complex);p=0
        for i in range(4):H[i,i]=y[p];p+=1
        for i in range(4):
            for j in range(i+1,4):
                H[i,j]=y[p]+1j*y[p+1];H[j,i]=np.conjugate(H[i,j]);p+=2
        N=expm(1j*H)@N0
    else:
        z=z0+x[:4];z=z/np.linalg.norm(z)
        K=np.zeros((4,4));p=4
        for i in range(4):
            for j in range(i+1,4):
                K[i,j]=x[p];K[j,i]=-x[p];p+=1
        N=expm(K)@N0
    return z,N

def val(s,x,complex_case=False):
    z,N=parameters(x,complex_case)
    W=np.eye(4)-(1-np.sqrt(1-s*s))*np.outer(z,np.conjugate(z))
    V=np.column_stack((s*z,W@N))
    B=V[:3,:];ll=V[3,:]
    J=ll@B@B@z[:3]
    return float(abs(J)**2)

if __name__=='__main__':
    rng=np.random.default_rng(2718);out={'boundary':val(0,np.zeros(10)), 'runs':[]}
    for cc in [False,True]:
        size=24 if cc else 10
        for s in [0.,.01,.02,.04,.07,.1,.2]:
            best=None
            for trial in range(4):
                x0=np.zeros(size) if trial==0 else rng.normal(size=size)*(.02 if trial<3 else .3)
                rr=minimize(lambda x:-val(s,x,cc),x0,method='BFGS',options={'gtol':1e-10,'maxiter':350})
                if best is None or rr.fun<best.fun:best=rr
            rec={'complex':cc,'s':s,'value':-best.fun,'slope':(-best.fun-7/16)/(s*s) if s else None,'success':bool(best.success),'gradnorm':float(np.linalg.norm(best.jac)),'nit':int(best.nit),'x':best.x.tolist()}
            print({k:v for k,v in rec.items() if k!='x'},flush=True);out['runs'].append(rec)
    open('research_v4_numeric/envelope_k3_results.json','w').write(json.dumps(out,indent=2))
