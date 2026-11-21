"""Bounded-precision positive matrix products. All decisions are integer."""
from dataclasses import dataclass
from fractions import Fraction as F

def ceildiv(a,b):return -((-a)//b)
def quantize(x,lo,hi,k):
    x=F(x);assert lo<=x<=hi
    Q=1<<k;v=(2*x.numerator*Q+x.denominator)//(2*x.denominator)
    v=max(ceildiv(lo.numerator*Q,lo.denominator),min(hi.numerator*Q//hi.denominator,v))
    assert abs(F(v,Q)-x)<=F(1,Q)
    return v

@dataclass(frozen=True)
class Box:
    lo:tuple
    hi:tuple
    depth:int

class Model:
    """Only immutable model constants; no history or oracle cache."""
    def __init__(self,p,snum,anum):
        assert p>=4
        self.p=p;self.k=p+12;self.q=p+p.bit_length()+20
        self.Q=1<<self.q;self.H=1<<self.k;self.snum=snum;self.anum=anum
        assert F(8,25)<=F(snum,self.H)<=F(17,50)
        assert F(37,50)<=F(anum,self.H)<=F(19,25)
        h=self.H;s=snum;a=anum
        m0=((h-a)*(h-s),a*s,(h-a)*s,a*(h-s))
        m1=(a*(h-s),(h-a)*s,a*s,(h-a)*(h-s))
        self.mats=(m0,m1);self.det=m0[0]*m0[3]-m0[1]*m0[2]
        assert self.det>0 and self.det==m1[0]*m1[3]-m1[1]*m1[2]
    @classmethod
    def from_rationals(cls,p,s,a):
        k=p+12
        return cls(p,quantize(s,F(8,25),F(17,50),k),quantize(a,F(37,50),F(19,25),k))
    def root(self):return Box((self.Q,0,0,self.Q),(self.Q,0,0,self.Q),0)
    def step(self,z,b):
        assert b in (0,1) and z.depth<=self.p
        a,b0,c,d=self.mats[b]
        def mul(x):return (x[0]*a+x[1]*c,x[0]*b0+x[1]*d,x[2]*a+x[3]*c,x[2]*b0+x[3]*d)
        lo=mul(z.lo);hi=mul(z.hi);W=sum(hi);Q=self.Q
        assert max((*lo,*hi,W)).bit_length()<=self.q+2*self.k+4
        assert (max(hi)*Q).bit_length()<=2*self.q+2*self.k+4
        ans=Box(tuple(v*Q//W for v in lo),tuple(ceildiv(v*Q,W) for v in hi),z.depth+1)
        assert min(ans.lo)>0 and max(ans.hi)<=Q and sum(ans.hi)<=Q+3
        return ans
    def point(self,z,num,den):
        a,b,c,d=z.lo;A,B,C,D=z.hi;Q=self.Q
        nl=num*a+(den-num)*b;nu=num*A+(den-num)*B
        dl=num*(a+c)+(den-num)*(b+d);du=num*(A+C)+(den-num)*(B+D)
        return max(3*Q//8,nl*Q//du),min(5*Q//8,ceildiv(nu*Q,dl))
    def enclosure(self,z):return self.point(z,3,8)[0],self.point(z,5,8)[1]
    def leaf(self,z):
        lo,hi=self.enclosure(z);return hi-lo<=2*(self.Q>>self.p)
    def dual(self,z):
        lo=(z.lo[0]+z.lo[2])*self.Q//sum(z.hi)
        hi=ceildiv((z.hi[0]+z.hi[2])*self.Q,sum(z.lo))
        return max(self.Q//8,min(7*self.Q//8,(lo+hi)//2))
    def output(self,z):
        lo,hi=self.enclosure(z) if self.leaf(z) else self.point(z,1,2)
        center=(lo+hi)//2
        v=(self.anum*self.Q+(self.H-2*self.anum)*center)//self.H
        assert 0<=v<=self.Q
        return v
