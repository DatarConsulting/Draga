"""Independent rational belief-coordinate checks. Enumeration is test-only."""
from fractions import Fraction as F
from pathlib import Path
import json,time,random
from certificate import certify
from cases import CASES
from model import Model
from codec import Codec

ROOT=Path(__file__).resolve().parents[1]

def belief(m):
    a,b,c,d=m;return a-b,b,a+c-b-d,b+d
def mul(x,y):
    a,b,c,d=x;A,B,C,D=y
    return a*A+b*C,a*B+b*D,c*A+d*C,c*B+d*D
def value(m,x):a,b,c,d=m;return (a*x+b)/(c*x+d)
def output(ms,x):a,b,c,d=ms[0];return (a+c)*x+(b+d)*(1-x)
def run_case(case,p,certificate=None):
    name,ms,J=case;start=time.perf_counter();cert=certificate or certify(ms,J=J)
    model=Model.from_rationals(p,cert,ms);codec=Codec(model)
    preparation=time.perf_counter()-start
    stored=tuple(tuple(F(x,model.H) for x in m) for m in model.mats)
    true_maps=tuple(map(belief,ms));stored_maps=tuple(map(belief,stored))
    eps=F(1,1<<p);rows=[];I=(F(1),F(0),F(0),F(1))
    def visit(w,z,zt,zh):
        lo,hi=model.enclosure(z);a,b=value(zh,J[0]),value(zh,J[1])
        assert F(lo,model.Q)<=a<=b<=F(hi,model.Q)
        assert max(a-F(lo,model.Q),F(hi,model.Q)-b)<=model.point_error_bound
        leaf=model.leaf(z);out=F(model.output(z),model.Q)
        xs=(value(zt,J[0]),value(zt,J[1])) if leaf else (value(zt,F(1,2)),)
        tv=max(abs(out-output(ms,x)) for x in xs)
        assert tv<=eps
        rows.append({'w':w,'z':z,'true':zt,'leaf':leaf,'tv':tv,'out':model.output(z)})
        if not leaf:
            assert len(w)<model.depth_cap
            for bit in (0,1):visit(w+str(bit),model.step(z,bit),mul(zt,true_maps[bit]),mul(zh,stored_maps[bit]))
    visit('',model.root(),I,I)
    words={r['w']:r for r in rows};ids={};minwidth=1<<codec.B
    for row in rows:
        i,width,z=codec.encode(row['w']);ids[row['w']]=i;minwidth=min(width,minwidth)
        assert z==row['z'] and codec.decode(i)[0]==row['w']
        assert codec.output(i)==row['out']
    assert len(set(ids.values()))==len(rows) and minwidth>8
    transitions=0
    for row in rows:
        for b in (0,1):
            raw=str(b)+row['w'];dest=raw
            if words['']['leaf']:dest=''
            else:
                for j in range(1,len(raw)+1):
                    if words[raw[:j]]['leaf']:dest=raw[:j];break
            target=words[dest]
            xs=(value(row['true'],J[0]),value(row['true'],J[1])) if row['leaf'] else (value(row['true'],F(1,2)),)
            transformed=[value(true_maps[b],x) for x in xs]
            lower,upper=value(target['true'],J[0]),value(target['true'],J[1])
            assert all(lower<=x<=upper for x in transformed)
            if row['leaf']:assert target['leaf']
            if not target['leaf']:assert transformed[0]==value(target['true'],F(1,2))
            i=codec.advance(ids[row['w']],b);assert i==ids[dest]
            assert codec.output(i)==target['out'];transitions+=1
    rng=random.Random(981+p);invalid=0;valid=set(ids.values())
    for _ in range(12):
        i=rng.randrange(1<<codec.B)
        try:w,_=codec.decode(i);assert i in valid and ids[w]==i
        except ValueError:assert i not in valid;invalid+=1
    ratio=max(r['tv']/eps for r in rows)
    result={'name':name,'p':p,'states':len(rows),'transitions':transitions,'outputs':len(rows)+transitions,
            'max_depth':max(len(r['w']) for r in rows),'depth_cap':model.depth_cap,'bits':codec.B,
            'constant_offset':cert.offset_bits,'k':model.k,'q':model.q,'determinants_distinct':model.dets[0]!=model.dets[1],
            'minimum_integer_width':minwidth,'max_TV_over_epsilon':str(ratio),
            'max_TV_over_epsilon_decimal':float(ratio),'sampled_unused_rejected':invalid,
            'oracle':codec.certificate,'preparation_seconds':preparation,'seconds':time.perf_counter()-start}
    (ROOT/'results'/('certificate_'+name+'.json')).write_text(json.dumps(cert.json(),indent=2)+'\n')
    return result

def main():
    schedule=[(case,p) for p in (4,8) for case in CASES]+[(CASES[i],12) for i in (1,3)]
    rows=[]
    for case,p in schedule:
        r=run_case(case,p);rows.append(r)
        print(json.dumps({k:r[k] for k in ['name','p','states','bits','max_TV_over_epsilon_decimal','seconds']}),flush=True)
        result={'passed':len(rows)==len(schedule),'complete':len(rows)==len(schedule),'configs':rows,
                'total_states':sum(x['states'] for x in rows),'total_transitions':sum(x['transitions'] for x in rows),
                'total_outputs':sum(x['outputs'] for x in rows),
                'scope':'All valid states and both transitions; unused IDs sampled, not exhausted. Exact true and quantized model references.'}
        (ROOT/'results/GRAPHS.json').write_text(json.dumps(result,indent=2)+'\n')
    print('ALL GRAPHS PASSED',flush=True)
if __name__=='__main__':main()
