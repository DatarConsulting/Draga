#!/usr/bin/env python3
"""Exact dyadic interval certificate: L_s h - h has a uniform strict sign.
All decisions use Python integers. Floats appear only in readable diagnostics.
"""
from fractions import Fraction as F
from pathlib import Path
import json
BITS=128; SCALE=1<<BITS
class I:
    def __init__(self,lo,hi=None): self.lo=lo; self.hi=lo if hi is None else hi
    @staticmethod
    def rational(x):
        x=F(x); a=x.numerator*SCALE; b=x.denominator
        return I(a//b,-((-a)//b))
    def __add__(self,b): return I(self.lo+b.lo,self.hi+b.hi)
    def __neg__(self): return I(-self.hi,-self.lo)
    def __sub__(self,b): return self+-b
    def __mul__(self,b):
        vs=[self.lo*b.lo,self.lo*b.hi,self.hi*b.lo,self.hi*b.hi]
        return I(min(vs)//SCALE,-((-max(vs))//SCALE))
    def __truediv__(self,b):
        assert b.lo>0
        vs=[F(self.lo*SCALE,b.lo),F(self.lo*SCALE,b.hi),F(self.hi*SCALE,b.lo),F(self.hi*SCALE,b.hi)]
        l,h=min(vs),max(vs)
        return I(l.numerator//l.denominator,-((-h.numerator)//h.denominator))
    def __pow__(self,n):
        z=I.rational(1)
        for _ in range(n):z=z*self
        return z
one=I.rational(1)
def log_i(v):
    z=(v-one)/(v+one); assert 0<=z.lo<=z.hi<SCALE
    z2=z*z; term=z; total=I.rational(0); K=32
    for k in range(K):
        total=total+term/I.rational(2*k+1); term=term*z2
    partial=I.rational(2)*total
    tail=I.rational(2)*term/(I.rational(2*K+1)*(one-z2))
    return I(partial.lo,partial.hi+tail.hi)
def exp_i(v):
    assert 0<=v.lo<=v.hi<2*SCALE
    total=one; term=one; N=40
    for k in range(1,N+1):
        term=term*v/I.rational(k); total=total+term
    tail=(term*v/I.rational(N+1))/(one-v/I.rational(N+2))
    return I(total.lo,total.hi+tail.hi)
def h(x):
    t=x-I.rational(F(1,2)); t2=t*t
    return one+I.rational(F(99577,100000))*t2+I.rational(F(125691,100000))*t2*t2

def residual(x,s):
    x=F(x); f0=F(1,3)+x/(3*(3-2*x)); f1=F(1,3)+x/(1+2*x)
    val=I.rational(0)
    for den,y in [(3-2*x,f0),(1+2*x,f1)]:
        weight=one/exp_i(I.rational(2*s)*log_i(I.rational(den)))
        val=val+weight*h(I.rational(y))
    return val-h(I.rational(x))

def main():
    grid=512; step=F(1,4*grid)
    # Analytically |(L_s h-h)''| <=24 throughout J and s in [1/2, 51/100].
    interpolation=24*step**2/8
    rows=[]
    for s,sign in [(F(50501,100000),1),(F(25251,50000),-1)]:
        vals=[residual(F(3,8)+k*step,s) for k in range(grid+1)]
        lo=min(v.lo for v in vals); hi=max(v.hi for v in vals)
        low=F(lo,SCALE)-interpolation; high=F(hi,SCALE)+interpolation
        assert low>0 if sign==1 else high<0
        rows.append({'s':str(s),'sign':sign,'grid_lower_units':lo,'grid_upper_units':hi,
                     'uniform_lower':str(low),'uniform_upper':str(high),
                     'uniform_lower_decimal':float(low),'uniform_upper_decimal':float(high)})
    result={'passed':True,'dimension_lower':'50501/100000','dimension_upper':'25251/50000',
            'precision_bits':BITS,'grid_intervals':grid,'second_derivative_bound':24,
            'interpolation_error':str(interpolation),'h_coefficients':['1','99577/100000','125691/100000'],
            'log_terms':32,'exp_terms':40,'certificates':rows}
    p=Path(__file__).resolve().parents[1]/'results/PRESSURE.json';p.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
