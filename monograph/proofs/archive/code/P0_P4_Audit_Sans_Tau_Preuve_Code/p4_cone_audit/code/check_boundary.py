"""A deterministic rounding-boundary fixture; exact stored source coefficients."""
from fractions import Fraction as F
from pathlib import Path
import json
from certificate import certify,hmm,point
from model import Model
from codec import Codec
ROOT=Path(__file__).resolve().parents[1]

def main():
    p=8;word='0110';J=(F(7,20),F(13,20));target=F(1,1<<(p+1))
    def matrices(a):return hmm(F(5,16),F(11,32),F(1,4),a)
    def width(a):
        ms=matrices(a);lo,hi=J
        for ch in reversed(word):lo=point(ms[int(ch)],lo);hi=point(ms[int(ch)],hi)
        return hi-lo
    low,high=F(45,64),F(49,64)
    assert width(low)>target>width(high)
    for _ in range(40):
        mid=(low+high)/2
        if width(mid)>target:low=mid
        else:high=mid
    ms=matrices(high);cert=certify(ms,J=J)
    k=max(x.denominator.bit_length()-1 for mat in ms for x in mat)
    assert all((x*(1<<k)).denominator==1 for mat in ms for x in mat)
    model=Model(p,cert,tuple(tuple(int(x*(1<<k)) for x in mat) for mat in ms),k,F(0))
    z=model.root()
    for ch in word:
        assert not model.leaf(z);z=model.step(z,int(ch))
    lo,hi=model.enclosure(z)
    assert width(high)<=target<F(hi-lo,model.Q)
    assert not model.leaf(z) and all(model.leaf(model.step(z,b)) for b in (0,1))
    co=Codec(model)
    for w in [word,word+'0',word+'1']:
        i,_,_=co.encode(w);assert co.decode(i)[0]==w
        for b in (0,1):
            new=co.advance(i,b);expected=co.encode(str(b)+w,truncate=True)[0];assert new==expected
    result={'passed':True,'p':p,'word':word,'emission_parameter':str(high),'coefficient_bits':k,'q':model.q,
            'exact_width':str(width(high)),'threshold':str(target),'upper_width':str(F(hi-lo,model.Q)),
            'exact_would_stop':True,'rounded_continues':True,'both_children_stop':True,'codec_transitions':6,
            'scope':'Stress fixture uses exactly dyadic coefficients at higher input precision, explicitly counted.'}
    (ROOT/'results/BOUNDARY.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='exact_width'}),flush=True)
if __name__=='__main__':main()
