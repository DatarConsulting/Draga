from reduced_codec import ReducedHybrid
from compile_comparison import ROOT,CASES
from certificate import certify,point
from model import Model
from audit_reduction import out,save
from fractions import Fraction as F

rows=[]
for name,p,mode in [('S',4,'grid'),('S',16,'fractal'),('A',24,'fractal')]:
    ms=dict(CASES)[name];H=certify(ms,J=(F(1,4),F(3,4)))
    h=ReducedHybrid(Model.from_rationals(p,H))
    assert h.mode==mode
    x=F(1,2);state=h.initial;peak=abs(h.output(state)-out(ms,x))
    for ch in '01101001':
        b=int(ch);state=h.advance(state,b);x=point(ms[b],x)
        err=abs(h.output(state)-out(ms,x));assert err<=F(1,1<<p)
        peak=max(peak,err)
    rows.append(dict(model=name,p=p,mode=mode,bits=h.B,steps=8,
                     max_TV_over_epsilon=float(peak*(1<<p))))
save('HYBRID.json',dict(passed=True,rows=rows))
print(rows)
