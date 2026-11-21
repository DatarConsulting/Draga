"""Outward integer intervals; inherited series with explicit remainder bounds."""
from fractions import Fraction as F

class Dyadic:
    def __init__(self,q):
        self.q=q;self.S=1<<q;self.one=(self.S,self.S)
        self.ln2=self.log_unit(self.rat(2))
    def rat(self,x,y=1):
        if isinstance(x,F):x,y=x.numerator,x.denominator
        a=x*self.S;return a//y,-((-a)//y)
    def add(self,a,b):return a[0]+b[0],a[1]+b[1]
    def sub(self,a,b):return a[0]-b[1],a[1]-b[0]
    def mul(self,a,b):
        v=[a[0]*b[0],a[0]*b[1],a[1]*b[0],a[1]*b[1]]
        return min(v)//self.S,-((-max(v))//self.S)
    def div(self,a,b):
        assert b[0]>0
        lo=min(x*self.S//y for x in a for y in b)
        hi=max(-((-x*self.S)//y) for x in a for y in b)
        return lo,hi
    def log_unit(self,a):
        # 1 <= a <= 2, log via atanh, with a positive geometric tail.
        z=self.div(self.sub(a,self.one),self.add(a,self.one));z2=self.mul(z,z)
        term=z;total=(0,0);K=self.q+8
        for k in range(K):
            total=self.add(total,self.div(term,self.rat(2*k+1)));term=self.mul(term,z2)
        total=self.mul(self.rat(2),total)
        tail=self.div(self.mul(self.rat(2),term),self.mul(self.rat(2*K+1),self.sub(self.one,z2)))
        return total[0],total[1]+tail[1]
    def log_ratio(self,num,den):
        assert num>=den>0
        e=max(0,num.bit_length()-den.bit_length())
        if num<den<<e:e-=1
        assert (den<<e)<=num<(den<<(e+1))
        m=self.rat(num,den<<e)
        return self.add(self.mul(self.rat(e),self.ln2),self.log_unit(m))
    def exp(self,a):
        assert a[0]>=0
        k=max(0,a[1].bit_length()-self.q)
        v=self.div(a,self.rat(1<<k));assert v[1]<=self.S
        total=self.one;term=self.one;N=self.q+8
        for j in range(1,N+1):
            term=self.div(self.mul(term,v),self.rat(j));total=self.add(total,term)
        tail=self.div(self.div(self.mul(term,v),self.rat(N+1)),self.sub(self.one,self.div(v,self.rat(N+2))))
        total=total[0],total[1]+tail[1]
        for _ in range(k):total=self.mul(total,total)
        return total
    def power_inverse(self,num,den,s):
        return self.div(self.one,self.exp(self.mul(self.rat(s),self.log_ratio(num,den))))

