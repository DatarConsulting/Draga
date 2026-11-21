"""Outward integer intervals. The operational power uses integer square roots."""
from fractions import Fraction as F

class Dyadic:
    def __init__(self,q):
        self.q=q;self.S=1<<q;self.one=(self.S,self.S)
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
