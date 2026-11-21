"""Streaming transfer oracle with model-derived rational certificates."""
from fractions import Fraction as F
from dyadic import Dyadic
from certificate import ceil,ceil_log2
from math import comb,isqrt

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
    gl=c.rat(model.cert.gm)[0];gh=c.rat(model.cert.grho)[1]
    # The true trajectory stays in Jdual, so this intersection is certified.
    weight=max(gl,min(g0[0],g1[0])),min(gh,max(g0[1],g1[1]))
    assert 0<=weight[0]<=weight[1]<=S
    return (low[0],high[1]),weight

def power(c,D,sigma):
    if D[0]<=0 or D[1]>c.S:raise NeedPrecision()
    sigma=F(sigma);a=sigma.numerator;b=sigma.denominator.bit_length()-1
    assert sigma>0 and sigma.denominator==1<<b,'compiled exponents must be dyadic'
    x=D
    for _ in range(b):
        lo=isqrt(x[0]*c.S);hi=isqrt(x[1]*c.S)
        if hi*hi<x[1]*c.S:hi+=1
        x=(lo,hi)
    out=c.one
    while a:
        if a&1:out=c.mul(out,x)
        a>>=1
        if a:x=c.mul(x,x)
    assert 0<=out[0]<=out[1]<=c.S
    return out

def precision_budget(model,n,degree,sigma,tolerance):
    """Conservative bit cap from the error-propagation lemma in the report.

    It is a safety assertion, not the precision normally used in practice.
    The actual outward widths decide termination. Integer sqrt avoids a
    transcendental-series precision assumption in the operational core.
    """
    H=model.cert;A=ceil(2*H.sigma_hi);b=F(sigma).denominator.bit_length()-1
    E=16+H.R**3+H.gm**(-(ceil(H.sigma_hi)+2))+1/(1-H.grho)+1/(1-H.taylor_radius)+1/H.beta+H.sigma_hi
    e=ceil_log2(E)
    needed=64+b+20*(n+degree+A+2)*e+ceil_log2(F((b+1)*(degree+1))/tolerance)
    q=48
    while q<needed:q*=2
    return q

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
    q=48;cap=precision_budget(model,n+1,0,sigma,tolerance)
    while True:
        assert q<=cap,'explicit ratio precision cap exceeded'
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

def certify_oracle_grid(model,n,sigma):
    cert=model.cert;target=F(1,64*model.depth_cap)
    ell=cert.Jdual[1]-cert.Jdual[0]
    # Always valid, without the old tau<1 condition.
    lip=2*cert.K
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

def taylor_parameters(model):
    cert=model.cert;target=F(1,64*model.depth_cap)
    r=cert.taylor_radius;A=ceil(2*cert.sigma_hi);degree=0
    factor=(1+r)**A
    assert 0<r<1 and A>=1
    while True:
        tail_ratio=r*F(degree+A+1,degree+2)
        if tail_ratio<1:
            tail=F(comb(A+degree,degree+1))*r**(degree+1)/(1-tail_ratio)
            if tail*factor<=target/65536:return degree,tail,factor
        degree+=1

def word_center(model,word,n,c):
    """At t=1/2 return derivative D_w and denominator log-slope u_w.

    Chain rule: u_(next o prefix)=u_prefix + u_next(x_prefix)*D_prefix.
    This uses bounded dyadic intervals, never exact long matrix products.
    """
    x=c.rat(F(1,2));D=c.one;u=(0,0)
    for j in range(n-1,-1,-1):
        bit=(word>>j)&1;A,B,C,E=model.mats[bit];slope=A+B-C-E
        def local(y):
            den=(A+B)*y+(C+E)*(c.S-y)
            if den<=0:raise NeedPrecision()
            return c.rat(slope*c.S,den)
        # slope/den is nonincreasing, irrespective of the sign of slope.
        low=local(x[1]);high=local(x[0]);v=(low[0],high[1])
        u=c.add(u,c.mul(v,D))
        x,g=trajectory_step(model,x,bit,c);D=c.mul(D,g)
    return D,u

def transfer_coefficients(model,n,sigma,degree,c):
    radius=max(F(1,2)-model.cert.Jdual[0],model.cert.Jdual[1]-F(1,2))
    result=[(0,0) for _ in range(degree+1)]
    for word in range(1<<n):
        D,u=word_center(model,word,n,c)
        term=power(c,D,sigma);v=c.mul(u,c.rat(-radius))
        for j in range(degree+1):
            result[j]=c.add(result[j],term)
            if j<degree:term=c.mul(c.mul(term,v),c.rat((2*sigma+j)/(j+1)))
    return result

def certify_oracle_taylor(model,n,sigma):
    """Certificate over the entire real interval, not a sampling test."""
    target=F(1,64*model.depth_cap)
    degree,tail,factor=taylor_parameters(model);q=48
    cap=precision_budget(model,n+1,degree,sigma,target/4096)
    while True:
        assert q<=cap,'explicit Taylor precision cap exceeded'
        c=Dyadic(q)
        try:
            a=transfer_coefficients(model,n,sigma,degree,c)
            b=transfer_coefficients(model,n+1,sigma,degree,c)
            if a[0][0]<=0:raise NeedPrecision()
            diffs=[c.sub(y,x) for x,y in zip(a,b)]
            uncertainty=F(sum(x[1]-x[0] for x in diffs),a[0][0])*factor
            relative_den_width=F(a[0][1]-a[0][0],a[0][0])
            if uncertainty>target/4096 or relative_den_width>target/4096:raise NeedPrecision()
            numerator=sum(max(abs(x[0]),abs(x[1])) for x in diffs)
            numerator+=tail*(a[0][1]+b[0][1])
            deviation=F(numerator,a[0][0])*factor
            if deviation>target/8:return None
            assert 2*deviation<target
            return {'method':'global_Taylor_interval','horizon':n,'sigma':str(sigma),
                    'taylor_degree':degree,'radius':str(model.cert.taylor_radius),
                    'tail_bound_per_center_mass':str(tail),'precision_bits':q,
                    'max_global_ratio_deviation':str(deviation),
                    'log_residual_bound':str(2*deviation),'target':str(target),
                    'arithmetic_width_bound':str(uncertainty)}
        except NeedPrecision:q*=2

def compile_oracle(model,progress=None):
    cert=model.cert;target=F(1,64*model.depth_cap)
    # Birkhoff cone diameter implies osc(log R_n)<=Delta*gamma^(n-1).
    nmax=1+cert.cone_half_life*ceil_log2(2048*cert.cone_diameter_bound/target)
    for n in range(nmax+1):
        sigma=center_root(model,n);out=certify_oracle_taylor(model,n,sigma)
        if out is None and n==nmax:
            out=certify_oracle_grid(model,n,sigma)
            if out is not None:out['method']='guaranteed_grid_fallback'
        if progress:progress({'n':n,'sigma':str(sigma),'accepted':out is not None,'cap':nmax})
        if out is not None:
            out['certified_horizon_cap']=nmax
            return sigma,n,out
    raise AssertionError('Projective contraction must force grid acceptance')

def branch(model,t,n,sigma):
    q=48;L=model.depth_cap;beta=model.cert.beta
    cap=precision_budget(model,n+1,0,sigma,beta/(512*L))
    Q=1<<ceil_log2(1024*L/beta)
    while True:
        assert q<=cap,'explicit branch precision cap exceeded'
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
