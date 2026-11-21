"""Replay a rounding-boundary fixture found by exact parameter bisection."""
from fractions import Fraction as F
from pathlib import Path
import json
from model import Model
from codec import Codec
from check_products import mult,project
p=16;m=Model(p,89478485,203538167);w='000001';z=m.root();exact=(1,0,0,1)
for b in map(int,w):
    assert not m.leaf(z)
    z=m.step(z,b);exact=mult(exact,m.mats[b])
width=project(exact,5,8)-project(exact,3,8);threshold=F(2,1<<p)
assert width<=threshold and not m.leaf(z)
assert all(m.leaf(m.step(z,b)) for b in (0,1))
co=Codec(m);checks=0
for v in [w,w+'0',w+'1']:
    i,_,_=co.encode(v);assert co.decode(i)[0]==v
    for b in (0,1):
        raw=str(b)+v;j,_,box=co.encode(raw,truncate=True)
        assert co.advance(i,b)==j
        if v!=w:assert m.leaf(box)
        checks+=1
out={'passed':True,'p':p,'stored_s':str(F(m.snum,m.H)),'stored_a':str(F(m.anum,m.H)),
     'word':w,'exact_width':str(width),'threshold':str(threshold),
     'upper_width':str(F(m.enclosure(z)[1]-m.enclosure(z)[0],m.Q)),
     'exact_rule_would_stop':True,'certified_rule_continues':True,'both_children_stop':True,
     'codec_transition_checks':checks,
     'note':'A real rounding-boundary case: the conservative tree deliberately keeps one more symbol.'}
(Path(__file__).resolve().parents[1]/'results/BOUNDARY.json').write_text(json.dumps(out,indent=2)+'\n')
print({k:v for k,v in out.items() if k!='exact_width'})
