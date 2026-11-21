"""Parametric bounded interval kernel. No history/cache is retained here."""
from dataclasses import dataclass
from fractions import Fraction as F
from certificate import ceil,ceil_log2,det,stochastic,check

def ceildiv(a,b):return -((-a)//b)
@dataclass(frozen=True)
class Box:
    lo:tuple
    hi:tuple
    depth:int

class Model:
    def __init__(self,p,certificate,matrices,k,parameter_error_bound):
        check(isinstance(p,int) and p>=4,'precision_p',p)
        self.p=p;self.cert=certificate;self.k=k;self.H=1<<k
        self.mats=tuple(tuple(m) for m in matrices)
        check(all(type(x) is int for m in self.mats for x in m),'integer_coefficients',k)
        rational=tuple(tuple(F(x,self.H) for x in m) for m in self.mats)
        check(stochastic(rational),'compiled_stochasticity',k)
        check(all(lo<=x<=hi for a,b,c in zip(certificate.lower,rational,certificate.upper) for lo,x,hi in zip(a,b,c)),'compiled_model_in_box',k)
        self.dets=tuple(det(m) for m in self.mats)
        self.parameter_error_bound=F(parameter_error_bound)
        eps=F(1,1<<p);C=max(certificate.parameter_TV_constant,certificate.pressure_stability_constant,F(1))
        check(self.parameter_error_bound<=eps/(8*C),'parameter_precision',str(self.parameter_error_bound))
        self.depth_cap=certificate.h*(p+2)
        L=self.depth_cap
        self.point_error_target=min(eps/64,(1-certificate.rho)*eps/32,
                                    certificate.beta/(512*L*(certificate.K+1)))
        need=(32*certificate.R**2*(L+1)+3)/self.point_error_target
        self.q=ceil_log2(need)+1;self.Q=1<<self.q
        self.point_error_bound=(32*certificate.R**2*(L+1)+3)/self.Q
        self.stop_units=self.Q>>(p+1)
        self.intermediate_cap=2*self.q+self.k+4
    @classmethod
    def from_rationals(cls,p,certificate,matrices=None):
        ms=certificate.center if matrices is None else tuple(tuple(map(F,m)) for m in matrices)
        check(stochastic(ms),'source_stochasticity','column sums')
        check(certificate.radius>0,'quantization_interior','a positive coefficient-box radius is required')
        check(all(abs(x-c)<=certificate.radius/2 for source,center in zip(ms,certificate.center) for x,c in zip(source,center)),
              'source_in_inner_box','source must lie in the half-radius box')
        C=max(certificate.parameter_TV_constant,certificate.pressure_stability_constant,F(1))
        k=max(p+4,ceil_log2(24*C*(1<<p)))
        k_limit=max(k,ceil_log2(6/certificate.radius))
        while k<=k_limit:
            H=1<<k;ns=[[0]*4,[0]*4]
            for col in (0,1):
                slots=[(0,col),(0,col+2),(1,col),(1,col+2)]
                for b,j in slots[:3]:ns[b][j]=(ms[b][j]*H).numerator//(ms[b][j]*H).denominator
                b,j=slots[3];ns[b][j]=H-sum(ns[c][i] for c,i in slots[:3])
            if all(lo<=F(x,H)<=hi for a,b,c in zip(certificate.lower,ns,certificate.upper) for lo,x,hi in zip(a,b,c)):
                assert max(abs(F(x,H)-y) for n,m in zip(ns,ms) for x,y in zip(n,m))<=F(3,H)
                return cls(p,certificate,ns,k,F(3,H))
            k+=1
        raise AssertionError('the certified interior margin forces quantization to terminate')
    def root(self):return Box((self.Q,0,0,self.Q),(self.Q,0,0,self.Q),0)
    def step(self,z,bit):
        assert bit in (0,1) and z.depth<=self.depth_cap
        a,b,c,d=self.mats[bit]
        def mul(x):return x[0]*a+x[1]*c,x[0]*b+x[1]*d,x[2]*a+x[3]*c,x[2]*b+x[3]*d
        low=mul(z.lo);high=mul(z.hi);W=sum(high);Q=self.Q
        assert max((*low,*high,W)).bit_length()<=self.q+self.k+4
        assert (max(high)*Q).bit_length()<=self.intermediate_cap
        ans=Box(tuple(v*Q//W for v in low),tuple(ceildiv(v*Q,W) for v in high),z.depth+1)
        assert min(ans.lo)>0 and max(ans.hi)<=Q and sum(ans.hi)<=Q+3
        return ans
    def point(self,z,x):
        x=F(x);i,j=x.numerator,x.denominator;a,b,c,d=z.lo;A,B,C,D=z.hi;Q=self.Q
        nl=i*a+(j-i)*b;nu=i*A+(j-i)*B
        dl=i*(a+c)+(j-i)*(b+d);du=i*(A+C)+(j-i)*(B+D)
        l,r=self.cert.J
        # The intersection uses outward dyadic endpoints for arbitrary rational J.
        return max((l*Q).numerator//(l*Q).denominator,nl*Q//du),min(ceil(r*Q),ceildiv(nu*Q,dl))
    def enclosure(self,z):return self.point(z,self.cert.J[0])[0],self.point(z,self.cert.J[1])[1]
    def leaf(self,z):
        lo,hi=self.enclosure(z);return hi-lo<=self.stop_units
    def dual(self,z):
        lo=(z.lo[0]+z.lo[2])*self.Q//sum(z.hi)
        hi=ceildiv((z.hi[0]+z.hi[2])*self.Q,sum(z.lo))
        mid=F((lo+hi)//2,self.Q)
        return max(self.cert.Jdual[0],min(self.cert.Jdual[1],mid))
    def output(self,z):
        lo,hi=self.enclosure(z) if self.leaf(z) else self.point(z,F(1,2))
        center=(lo+hi)//2;a,b,c,d=self.mats[0]
        v=((a+c)*center+(b+d)*(self.Q-center))//self.H
        assert 0<=v<=self.Q
        return v
