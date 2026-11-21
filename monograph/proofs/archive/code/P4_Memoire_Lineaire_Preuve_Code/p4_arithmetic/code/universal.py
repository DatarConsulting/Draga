"""Certified arbitrary-p arithmetic code. Integers only; no count table.

Finite transfer sums are streamed. The implementation deliberately rebuilds
each short product, avoiding a stack of long rational projective states.
The universal backend is a correctness prototype, not the fast C backend.
"""
from fractions import Fraction as F
from fast import append,leaf

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

def horizon(p,factor=16):
    n=0
    while factor*p*3**n>5**n:n+=1
    return n

def transfer(tnum,tden,n,s,ctx,first=None):
    """U_n(t), or the b-prefixed contribution to U_(n+1)(t)."""
    acc=(0,0);depth=n+(first is not None)
    for word in range(1<<n):
        z=(1,0,0,1)
        if first is not None:z=append(z,first)
        for j in range(n-1,-1,-1):z=append(z,(word>>j)&1)
        a,b,c,d=z;den=tnum*(a+b)+(tden-tnum)*(c+d)
        term=ctx.power_inverse(den*den,9**depth*tden*tden,s)
        acc=ctx.add(acc,term)
    return acc

def choose_s(p):
    return refine_s(p,F(50501,100000),F(25251,50000))

def refine_s(p,lo,hi):
    """Precondition: [lo,hi] lies in [1/2,51/100] and contains d."""
    assert F(1,2)<=lo<hi<=F(51,100)
    while hi-lo>F(1,8*p):
        s=(lo+hi)/2;n=horizon(p,128);q=64
        while True:
            c=Dyadic(q);a=transfer(1,2,n,s,c);b=transfer(1,2,n+1,s,c)
            if a[0]<=0:q*=2;continue
            r=c.div(b,a);rlo,rhi=F(r[0],c.S),F(r[1],c.S)
            if rhi-rlo<=F(1,2048*p) and F(3,4)<=rlo<=rhi<=F(5,4):break
            q*=2
        eta=F(1,128*p);lower=rlo*(1-eta);upper=rhi*(1+eta)
        if upper<1:hi=s
        elif lower>1:lo=s
        else:
            assert 2*(upper-lower)<F(1,16*p)
            return s
    return (lo+hi)/2

class Codec:
    def __init__(self,p):
        assert p>=4
        self.p=p;self.s=choose_s(p)
        self.n=horizon(p);self.B=-((-self.s.numerator*p)//self.s.denominator)+5
        self.qprob=(2048*p-1).bit_length();self.Q=1<<self.qprob
    def branch(self,S,T):
        q=64
        while True:
            c=Dyadic(q);a=transfer(S,S+T,self.n,self.s,c,0);b=transfer(S,S+T,self.n,self.s,c,1)
            if a[0]+b[0]<=0:q*=2;continue
            r=c.div(a,c.add(a,b));lo,hi=F(r[0],c.S),F(r[1],c.S)
            if hi-lo<=F(1,1024*self.p):break
            q*=2
        mid=(lo+hi)/2;k=(mid.numerator*self.Q*2+mid.denominator)//(2*mid.denominator)
        assert max(abs(F(k,self.Q)-lo),abs(F(k,self.Q)-hi))<=F(1,256*self.p)
        assert F(1,5)<=F(k,self.Q)<=F(4,5)
        return k
    def encode(self,w):
        lo,N,z=0,1<<self.B,(1,0,0,1)
        for b in map(int,w):
            assert not leaf(z,self.p)
            k=self.branch(z[0]+z[2],z[1]+z[3]);left=(N-1)*k//self.Q
            assert 0<left<N-1
            lo+=1+(left if b else 0);N=N-1-left if b else left;z=append(z,b)
        return lo,N
    def decode(self,i):
        if not 0<=i<1<<self.B:raise ValueError('invalid ID')
        lo,N,z,w=0,1<<self.B,(1,0,0,1),''
        while i!=lo:
            if leaf(z,self.p):raise ValueError('unused ID')
            k=self.branch(z[0]+z[2],z[1]+z[3]);left=(N-1)*k//self.Q
            assert 0<left<N-1
            b=int(i>=lo+1+left);lo+=1+(left if b else 0);N=N-1-left if b else left
            w+=str(b);z=append(z,b)
        return w
    def advance(self,i,b):
        assert b in (0,1)
        raw=str(b)+self.decode(i);z=(1,0,0,1);w=''
        for c in raw:
            z=append(z,int(c));w+=c
            if leaf(z,self.p):break
        return self.encode(w)[0]
    def output(self,i):
        z=(1,0,0,1)
        for b in map(int,self.decode(i)):z=append(z,b)
        a,b,c,d=z
        if leaf(z,self.p):
            nl,nr=3*a+5*b,5*a+3*b;dl,dr=nl+3*c+5*d,nr+5*c+3*d
            return F(3*dl*dr-nl*dr-nr*dl,4*dl*dr)
        return F(3*(a+b+c+d)-2*(a+b),4*(a+b+c+d))
