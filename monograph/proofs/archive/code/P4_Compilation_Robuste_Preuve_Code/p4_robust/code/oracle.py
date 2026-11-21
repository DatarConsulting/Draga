"""Streaming positive transfer sums with bounded-precision dual trajectories."""
from fractions import Fraction as F
from dyadic import Dyadic

class NeedPrecision(Exception):pass

def trajectory_step(model,t,b,c):
    A,B,C,D=model.mats[b];Q=c.S
    def endpoint(x):
        den=(A+B)*x+(C+D)*(Q-x)
        new=c.rat(A*x+C*(Q-x),den)
        weight=c.rat(model.det*Q*Q,den*den)
        assert (model.det*Q**3).bit_length()<=4*model.k+3*c.q+1
        return new,weight
    low,g0=endpoint(t[0]);high,g1=endpoint(t[1])
    return (low[0],high[1]),(min(g0[0],g1[0]),max(g0[1],g1[1]))

def power(c,D,sigma):
    if D[0]<=0:raise NeedPrecision()
    assert D[1]<=c.S
    low=c.log_ratio(c.S,D[1]);high=c.log_ratio(c.S,D[0])
    return c.div(c.one,c.exp(c.mul(c.rat(sigma),(low[0],high[1]))))

def transfer(model,tnum,tden,n,sigma,c,first=None):
    total=(0,0)
    for word in range(1<<n):
        t=c.rat(tnum,tden);D=c.one
        if first is not None:
            t,g=trajectory_step(model,t,first,c);D=c.mul(D,g)
        for j in range(n-1,-1,-1):
            t,g=trajectory_step(model,t,(word>>j)&1,c);D=c.mul(D,g)
        total=c.add(total,power(c,D,sigma))
    return total

def ratio(model,tnum,tden,n,sigma,tolerance):
    q=48
    while True:
        c=Dyadic(q)
        try:
            a=transfer(model,tnum,tden,n,sigma,c)
            b=transfer(model,tnum,tden,n+1,sigma,c)
            if a[0]<=0:raise NeedPrecision()
            r=c.div(b,a)
            if F(r[1]-r[0],c.S)<=tolerance:return r,c
        except NeedPrecision:pass
        q*=2

def center_root(model,n):
    lo,hi=F(1,3),F(3,4);p=model.p
    tol=F(1,8192*p)
    while True:
        mid=(lo+hi)/2;r,c=ratio(model,1,2,n,mid,tol)
        if r[0]<=c.S<=r[1]:return mid
        if hi-lo<=F(1,16384*p*(n+1)):return mid
        if r[0]>c.S:lo=mid
        else:hi=mid

def certify(model,n,sigma):
    p=model.p;M=64*p;target=F(1,128*p);worst=F(0);maxq=0
    for j in range(M+1):
        t=F(1,8)+F(3*j,4*M)
        r,c=ratio(model,t.numerator,t.denominator,n,sigma,F(1,8192*p))
        bound=max(abs(F(r[0],c.S)-1),abs(F(r[1],c.S)-1))
        if bound>target:return None
        worst=max(worst,bound);maxq=max(maxq,c.q)
    # Between nodes Lip(log R)<=9/4. The logarithm at a grid node is
    # bounded by 2*|R-1|. The grid is streamed, not stored.
    residual=2*worst+F(27,32*M)
    assert residual<F(1,32*p)
    return {'horizon':n,'sigma':str(sigma),'grid_intervals':M,
            'max_grid_ratio_deviation':str(worst),'log_residual_bound':str(residual),
            'claimed_log_residual_bound':str(F(1,32*p)),'precision_bits':maxq}

def compile_oracle(model,progress=None):
    nmax=0
    # 2*(27/32)*theta^n < 1/(2048p), ample space for center error and grid.
    while 4096*model.p*4**nmax>5**nmax:nmax+=1
    for n in range(nmax+1):
        sigma=center_root(model,n);cert=certify(model,n,sigma)
        if progress:progress({'n':n,'sigma':str(sigma),'accepted':cert is not None})
        if cert is not None:return sigma,n,cert
    raise AssertionError('Uniform convergence bound should force acceptance')

def branch(model,tnum,tden,n,sigma):
    q=48;p=model.p;Q=1<<(2048*p-1).bit_length()
    while True:
        c=Dyadic(q)
        try:
            a=transfer(model,tnum,tden,n,sigma,c,0)
            b=transfer(model,tnum,tden,n,sigma,c,1)
            if a[0]+b[0]<=0:raise NeedPrecision()
            r=c.div(a,c.add(a,b))
            if F(r[1]-r[0],c.S)<=F(1,1024*p):break
        except NeedPrecision:pass
        q*=2
    k=((r[0]+r[1])*Q+c.S)//(2*c.S)
    assert F(1,10)<=F(k,Q)<=F(9,10)
    return k,Q
