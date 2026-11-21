"""Short interface smoke checks where the selector actually chooses the codec."""
from audit_direct_grid import ROOT, CASES, out, save
from direct_grid import HybridCodec
from certificate import certify, point
from model import Model
from fractions import Fraction as F

rows=[]
for name,p in [('S',128),('S',256),('A',256)]:
    ms=dict(CASES)[name]
    c=certify(ms,J=(F(1,4),F(3,4)))
    h=HybridCodec(Model.from_rationals(p,c))
    assert h.mode=='fractal'
    state=h.initial; x=F(1,2); peak=abs(h.output(state)-out(ms,x))
    for ch in '01101001':
        b=int(ch); x=point(ms[b],x); state=h.advance(state,b)
        e=abs(h.output(state)-out(ms,x)); assert e<=F(1,1<<p)
        peak=max(peak,e)
    rows.append(dict(model=name,p=p,bits=h.B,steps=8,max_TV_over_epsilon=float(peak*(1<<p))))
    print(rows[-1],flush=True)
save('HYBRID_BRANCH.json',dict(passed=True,rows=rows,scope='Eight updates per case, from the initial state; interface smoke tests only, no deep-state or all-history verification.'))
