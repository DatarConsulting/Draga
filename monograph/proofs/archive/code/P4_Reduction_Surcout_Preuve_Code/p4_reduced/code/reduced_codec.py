"""Integer prefix codec certified on internal nodes; no tree precomputation.

All geometry, stopping, output, and transfer evaluations are inherited.
Only the scalar capacity certificate and the integer split change.
"""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'reference_code'))
from fractions import Fraction as F
from certificate import ceil,ceil_log2
from codec import Codec
from dyadic import Dyadic
from oracle import (transfer,transfer_coefficients,power,precision_budget,
                    NeedPrecision)


def branch_enclosure(model,t,n,sigma,tolerance):
    q=48;cap=precision_budget(model,n+1,0,sigma,tolerance)
    while True:
        assert q<=cap
        c=Dyadic(q)
        try:
            a=transfer(model,t,n,sigma,c,0)
            b=transfer(model,t,n,sigma,c,1)
            if a[0]+b[0]<=0:raise NeedPrecision()
            r=c.div(a,c.add(a,b))
            if F(r[1]-r[0],c.S)<=tolerance:
                return F(r[0],c.S),F(r[1],c.S),q
        except NeedPrecision:pass
        q*=2


def shape_lower_bound(model,n,sigma):
    # Every summand D_w(t)^sigma is convex. Its sum lies above its
    # tangent at 1/2. a[1] already includes the interval radius.
    q=48;tol=F(1,4096)
    cap=precision_budget(model,n,1,sigma,tol)
    while True:
        assert q<=cap
        c=Dyadic(q)
        try:
            a=transfer_coefficients(model,n,sigma,1,c)
            if a[0][0]<=0:raise NeedPrecision()
            width=F(sum(hi-lo for lo,hi in a),a[0][0])
            if width<=tol:break
        except NeedPrecision:pass
        q*=2
    slope_ratio=F(max(abs(a[1][0]),abs(a[1][1])),a[0][0])
    tangent=1-slope_ratio
    fallback=(1+model.cert.taylor_radius)**(-ceil(2*model.cert.sigma_hi))
    return max(tangent,fallback),dict(precision_bits=q,
        normalized_slope_absolute_upper=str(slope_ratio),
        tangent_lower=str(tangent),fallback_lower=str(fallback))


def positive_power_lower(model,x,sigma):
    assert 0<x<=1
    lower0=x**ceil(model.cert.sigma_hi)
    tolerance=lower0/4096
    q=48;cap=precision_budget(model,0,0,sigma,tolerance)
    # The new base x can be much smaller than a single map derivative.
    # Include its own root-chain sensitivity; do not reuse only the old
    # transfer budget for arbitrarily small invariant-interval endpoints.
    a=F(sigma).numerator;b=F(sigma).denominator.bit_length()-1
    need=max(ceil_log2(2/x)+2,
             ceil_log2(a*(4*(b+1)/x+2)/tolerance)+2)
    while cap<need:cap*=2
    while True:
        assert q<=cap
        c=Dyadic(q)
        try:
            a=power(c,c.rat(x),sigma)
            if a[0]>0 and F(a[1]-a[0],c.S)<=tolerance:
                return F(a[0],c.S),q
        except NeedPrecision:pass
        q*=2


def compile_capacity(model,n,sigma,grid_intervals=16):
    H=model.cert;L=model.depth_cap;eps=F(1,1<<model.p)
    l,r=H.J
    def z(t,x):return t*x+(1-t)*(1-x)
    # l<1/2<r makes this quadratic concave: endpoints give its minimum.
    geometric=min(4*z(t,l)*z(t,r) for t in H.Jdual)
    interior=min(F(1),geometric*(F(1,2)-2*model.point_error_bound/eps)/(r-l))
    assert interior>0
    pmin,pq=positive_power_lower(model,interior,sigma)
    vmin,shape=shape_lower_bound(model,n,sigma)
    # Old transfer-ratio and branch-rounding proofs each lose <=1/64
    # in log mass along a path. exp(-1/32) >= 31/32.
    mass=pmin*vmin*F(31,32)
    assert 0<mass<=1

    # A fixed-size grid plus a certified Lipschitz radius, not sampling.
    J0,J1=H.Jdual;raw=F(1,2);maxq=0
    for j in range(grid_intervals+1):
        t=J0+(J1-J0)*F(j,grid_intervals)
        lo,hi,q=branch_enclosure(model,t,n,sigma,H.beta/4096)
        raw=min(raw,lo,1-hi);maxq=max(maxq,q)
    interpolation=H.K*(J1-J0)/(4*grid_intervals)
    exact_beta=max(H.beta,raw-interpolation)
    # Conservatively includes the inherited midpoint and branch-grid errors.
    beta_hat=exact_beta-H.beta/(128*L)
    assert 0<beta_hat<=F(1,2)
    # N(u) >= 2^B mu_hat(u) - (3/2-beta_hat)/beta_hat.
    # If 2^B mu_hat(u)>2/beta_hat for every INTERNAL u, both children
    # receive >=1 under nearest rounding. No leaf-mass margin is needed.
    threshold=2/(beta_hat*mass)
    offset=ceil_log2(threshold)
    if (1<<offset)<=threshold:offset+=1
    B=ceil(sigma*model.p)+offset
    return B,dict(method='internal_mass_nearest_split',offset_bits=offset,
        bits=B,geometric_denominator_lower=str(geometric),
        internal_delta_over_epsilon_lower=str(interior),
        power_lower=str(pmin),power_precision_bits=pq,shape_lower=str(vmin),
        shape_certificate=shape,internal_mass_coefficient_lower=str(mass),
        branch_grid_intervals=grid_intervals,branch_grid_min_lower=str(raw),
        branch_interpolation_loss=str(interpolation),branch_exact_min_lower=str(exact_beta),
        rounded_branch_min_lower=str(beta_hat),branch_precision_bits=maxq,
        required_root_factor=str(threshold),certified_root_factor=1<<offset,
        worst_capacity_loss=str((F(3,2)-beta_hat)/beta_hat),
        fixed_model_certificate=True,tree_enumerated=False)


class ReducedCodec(Codec):
    def __init__(self,model,progress=None):
        super().__init__(model,progress)
        self.old_B=self.B
        self.B,self.capacity_certificate=compile_capacity(model,self.n,self.sigma)

    @staticmethod
    def split(N,k,Q):
        left=(2*(N-1)*k+Q)//(2*Q)
        assert 0<left<N-1,'certified capacity violated'
        return left

    def encode(self,w,truncate=False):
        model=self.model;z=model.root();lo,N=0,1<<self.B
        for ch in w:
            if ch not in ('0','1'):raise ValueError('non-binary word')
            if model.leaf(z):
                if truncate:break
                raise ValueError('path extends beyond leaf')
            b=int(ch);k,Q=self.probability(z);left=self.split(N,k,Q)
            lo+=1+(left if b else 0);N=N-1-left if b else left
            z=model.step(z,b)
        return lo,N,z

    def decode(self,i):
        if type(i) is not int or not 0<=i<1<<self.B:raise ValueError('invalid ID')
        model=self.model;z=model.root();lo,N=0,1<<self.B;w=''
        while i!=lo:
            if model.leaf(z):raise ValueError('unused ID')
            k,Q=self.probability(z);left=self.split(N,k,Q)
            b=int(i>=lo+1+left);lo+=1+(left if b else 0)
            N=N-1-left if b else left;z=model.step(z,b);w+=str(b)
        return w,z


class ReducedHybrid:
    """Static choice between the certified direct grid and the reduced codec."""
    def __init__(self,model):
        from direct_grid import DirectGrid
        grid=DirectGrid(model);fractal=ReducedCodec(model)
        self.model=model
        self.grid_bits,self.fractal_bits=grid.B,fractal.B
        self.mode='grid' if grid.B<=fractal.B else 'fractal'
        self.engine=grid if self.mode=='grid' else fractal
        self.B=self.engine.B
        self.initial=grid.initial if self.mode=='grid' else 0

    def advance(self,state,symbol):
        return self.engine.advance(state,symbol)

    def output(self,state):
        if self.mode=='grid':return self.engine.output(state)
        return F(self.engine.output(state),self.model.Q)
