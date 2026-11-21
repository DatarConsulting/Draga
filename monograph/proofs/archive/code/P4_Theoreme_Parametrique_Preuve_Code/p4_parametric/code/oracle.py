"""Streaming transfer oracle with model-derived rational certificates."""
from fractions import Fraction as F
from dyadic import Dyadic
from certificate import ceil,ceil_log2

class NeedPrecision(Exception):pass

def trajectory_step(model,t,bit,c):
    A,B,C,D=model.mats[bit];S=c.S;de=model.dets[bit]
    def endpoint(x):
        den=(A+B)*x+(C+D)*(S-x)
        assert den>0
        new=c.rat(A*x+C*(S-x),den)
        weight=c.rat(de*S*S,den*den)
        assert (de*S**3).bit_length()<=2*model.k+3*c.q+2
        return new,weight
    low,g0=endpoint(t[0]);high,g1=endpoint(t[1])
    return (low[0],high[1]),(min(g0[0],g1[0]),max(g0[1],g1[1]))

def power(c,D,sigma):
    if D[0]<=0 or D[1]>c.S:raise NeedPrecision()
    low=c.log_ratio(c.S,D[1]);high=c.log_ratio(c.S,D[0])
    return c.div(c.one,c.exp(c.mul(c.rat(sigma),(low[0],high[1]))))

def transfer(model,t,n,sigma,c,first=None):
    total=(0,0)
    for word in range(1<<n):
        x=c.rat(t);D=c.one
        if first is not None:
            x,g=trajectory_step(model,x,first,c);D=c.mul(D,g)
        for j in range(n-1,-1,-1):
            x,g=trajectory_step(model,x,(word>>j)&1,c);D=c.mul(D,g)
        total=c.add(total,power(c,D,sigma))
    return total

def ratio(model,t,n,sigma,tolerance):
    q=48
    while True:
        c=Dyadic(q)
        try:
            a=transfer(model,t,n,sigma,c);b=transfer(model,t,n+1,sigma,c)
            if a[0]<=0:raise NeedPrecision()
            r=c.div(b,a)
            if F(r[1]-r[0],c.S)<=tolerance:return r,c
        except NeedPrecision:pass
        q*=2

def center_root(model,n):
    cert=model.cert;lo,hi=cert.sigma_lo,cert.sigma_hi
    target=F(1,64*model.depth_cap);U=ceil_log2(1/cert.gm)
    while True:
        mid=(lo+hi)/2;r,c=ratio(model,F(1,2),n,mid,target/512)
        if r[0]<=c.S<=r[1]:return mid
        if hi-lo<=target/(256*U*(n+1)):return mid
        if r[0]>c.S:lo=mid
        else:hi=mid

def certify_oracle(model,n,sigma):
    cert=model.cert;target=F(1,64*model.depth_cap)
    ell=cert.Jdual[1]-cert.Jdual[0]
    lip=cert.K*cert.tau**n
    M=max(1,ceil(lip*ell/target));worst=F(0);maxq=0
    for j in range(M+1):
        t=cert.Jdual[0]+ell*F(j,M)
        r,c=ratio(model,t,n,sigma,target/512)
        bound=max(abs(F(r[0],c.S)-1),abs(F(r[1],c.S)-1))
        if bound>target/8:return None
        worst=max(worst,bound);maxq=max(maxq,c.q)
    residual=2*worst+lip*ell/(2*M)
    assert residual<target
    return {'horizon':n,'sigma':str(sigma),'grid_intervals':M,
            'log_ratio_lipschitz_bound':str(lip),'max_grid_ratio_deviation':str(worst),
            'log_residual_bound':str(residual),'target':str(target),'precision_bits':maxq}

def compile_oracle(model,progress=None):
    cert=model.cert;target=F(1,64*model.depth_cap)
    nmax=0;error=2*cert.K*(cert.Jdual[1]-cert.Jdual[0])
    while error>target/1024:
        error*=cert.tau;nmax+=1
    for n in range(nmax+1):
        sigma=center_root(model,n);out=certify_oracle(model,n,sigma)
        if progress:progress({'n':n,'sigma':str(sigma),'accepted':out is not None})
        if out is not None:return sigma,n,out
    raise AssertionError('Quantitative contraction should force acceptance')

def branch(model,t,n,sigma):
    q=48;L=model.depth_cap;beta=model.cert.beta
    Q=1<<ceil_log2(1024*L/beta)
    while True:
        c=Dyadic(q)
        try:
            a=transfer(model,t,n,sigma,c,0);b=transfer(model,t,n,sigma,c,1)
            if a[0]+b[0]<=0:raise NeedPrecision()
            r=c.div(a,c.add(a,b))
            if F(r[1]-r[0],c.S)<=beta/(512*L):break
        except NeedPrecision:pass
        q*=2
    k=((r[0]+r[1])*Q+c.S)//(2*c.S)
    assert beta/2<=F(k,Q)<=1-beta/2
    return k,Q
