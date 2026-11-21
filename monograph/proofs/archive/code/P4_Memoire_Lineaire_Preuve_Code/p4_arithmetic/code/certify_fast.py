"""Whole-domain certificate for the exact C fixed-point expression.
Integer interval propagation encloses every integer t24 in each cell.
Comparison uses outward dyadic intervals for the ideal probability at its
cell midpoint and its proven Lipschitz constant 3/10. No float decision.
"""
import sys,json,time
from pathlib import Path
from fractions import Fraction as F
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'reference'))
from certify_pressure import I,SCALE,log_i,exp_i,h
from fast import Q,P,LN2Q,trunc
class R:
    def __init__(self,a,b=None):self.a=a;self.b=a if b is None else b
    def __add__(self,o):
        o=o if isinstance(o,R) else R(o);return R(self.a+o.a,self.b+o.b)
    __radd__=__add__
    def __neg__(self):return R(-self.b,-self.a)
    def __sub__(self,o):return self+-o
    def __rsub__(self,o):return -self+o
    def __mul__(self,o):
        o=o if isinstance(o,R) else R(o);v=[self.a*o.a,self.a*o.b,self.b*o.a,self.b*o.b]
        assert max(map(abs,v))<1<<63
        return R(min(v),max(v))
    __rmul__=__mul__
    def __floordiv__(self,o):
        o=o if isinstance(o,R) else R(o);assert o.a>0
        v=[trunc(x,y) for x in [self.a,self.b] for y in [o.a,o.b]]
        return R(min(v),max(v))
    def square(self):
        return R(0 if self.a<=0<=self.b else min(self.a**2,self.b**2),max(self.a**2,self.b**2))
def weight(t,b):
    den=5*Q+2*t if b else 7*Q-2*t
    root=R(3*Q*Q)//den;child=(3 if b else 1)*Q*(Q+t)//den
    z=(den-6*Q)*Q//(den+6*Q);z2=z.square()//Q;z3=z2*z//Q;z5=z3*z2//Q
    logv=LN2Q+2*(z+z3//3+z5//5)
    x=2*1003*logv//200000;e=Q-x+x.square()//(2*Q);a=root*e//Q
    v=(child-Q//2)//3;v2=v.square()//Q;v4=v2.square()//Q
    H=Q+99577*v2//100000+125691*v4//100000
    return a*H//Q
def ideal(t):
    ws=[]
    for b in [0,1]:
        den=5+2*t if b else 7-2*t
        child=(3 if b else 1)*(1+t)/den
        a=I.rational(1)/exp_i(I.rational(F(101003,100000))*log_i(I.rational(den/3)))
        ws.append(a*h(I.rational((child+1)/3)))
    return ws[0]/(ws[0]+ws[1])
def main():
    start=time.perf_counter();cells=16384;span=3*Q//4;worst=F(0);min_k=P;max_k=0
    ln=log_i(I.rational(2));assert ln.lo//(SCALE//Q)==LN2Q==ln.hi//(SCALE//Q)
    for i in range(cells):
        a=Q//8+span*i//cells;b=Q//8+span*(i+1)//cells
        t=R(a,b);w0=weight(t,0);w1=weight(t,1);v=w0*P//(w0+w1)
        min_k=min(min_k,v.a);max_k=max(max_k,v.b)
        mid=F(a+b,2*Q);pi=ideal(mid)
        # True t lies within [a/Q,(b+1)/Q], because the input was floored.
        margin=F(3,10)*F(b-a+2,2*Q)
        err=max(abs(F(v.a,P)-F(pi.hi,SCALE)),abs(F(v.b,P)-F(pi.lo,SCALE)))+margin
        worst=max(worst,err)
        assert err<F(1,8000),(i,float(err))
        if i%4096==0:print('certificate cells',i,flush=True)
    assert F(min_k,P)>=F(2,5) and F(max_k,P)<=F(3,5)
    out={'passed':True,'domain':['1/8','7/8'],'cells':cells,'t_fraction_bits':24,
         'probability_fraction_bits':16,'ideal_lipschitz_bound':'3/10',
         'max_error_upper_rational':str(worst),'max_error_upper_decimal':float(worst),
         'claimed_error_bound':'1/8000','branch_min_bound':str(F(min_k,P)),
         'branch_max_bound':str(F(max_k,P)),'all_modeled_signed_products_below_2_63':True,
         'seconds_local':time.perf_counter()-start}
    (ROOT/'results/FAST_CERTIFICATE.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
