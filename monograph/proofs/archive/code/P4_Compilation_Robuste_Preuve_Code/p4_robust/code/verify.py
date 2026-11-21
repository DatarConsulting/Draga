"""Exhaustive small graphs, independent rational belief-coordinate oracle.
All enumeration and collected metrics are exclusively in this verifier.
"""
from pathlib import Path
from fractions import Fraction as F
import json,time,random
from model import Model
from codec import Codec
ROOT=Path(__file__).resolve().parents[1]

def matrices(s,a):
    return ((1-s-a,s*a,1-2*a,a),(a-s,s*(1-a),2*a-1,1-a))
def append(z,m):
    a,b,c,d=z;A,B,C,D=m
    return a*A+b*C,a*B+b*D,c*A+d*C,c*B+d*D
def point(z,x):
    a,b,c,d=z;return (a*x+b)/(c*x+d)
def bounds(z):return point(z,F(3,8)),point(z,F(5,8))
def bayes(x,b,s,a):return point(matrices(s,a)[b],x)
def predict(x,a):return a+(1-2*a)*x

def verify_one(p,s,a):
    start=time.perf_counter();m=Model.from_rationals(p,s,a);prep=time.perf_counter()
    co=Codec(m);prep=time.perf_counter()-prep
    sh,ah=F(m.snum,m.H),F(m.anum,m.H);mt=matrices(s,a);mh=matrices(sh,ah)
    eps=F(1,1<<p);delta=max(abs(s-sh),abs(a-ah));rows=[]
    def visit(w,z,zt,zh):
        true=bounds(zt);hat=bounds(zh);lo,hi=m.enclosure(z);err=F(1568*z.depth+1,m.Q)
        assert F(lo,m.Q)<=hat[0]<=hat[1]<=F(hi,m.Q)
        assert hat[0]-F(lo,m.Q)<=err and F(hi,m.Q)-hat[1]<=err
        leaf=m.leaf(z);out=F(m.output(z),m.Q)
        tv=max(abs(out-predict(x,a)) for x in true) if leaf else abs(out-predict(point(zt,F(1,2)),a))
        assert tv<=eps
        # Independent analytic TV budget including coefficient quantization.
        assert tv<=F(13,25)*eps+F(1,256)*eps+F(9,4)*delta
        row={'word':w,'box':z,'true':true,'hat':hat,'center':point(zt,F(1,2)),
             'leaf':leaf,'output':m.output(z),'tv':tv}
        rows.append(row)
        if not leaf:
            assert len(w)<p
            for b in (0,1):visit(w+str(b),m.step(z,b),append(zt,mt[b]),append(zh,mh[b]))
    visit('',m.root(),(F(1),F(0),F(0),F(1)),(F(1),F(0),F(0),F(1)))
    byword={z['word']:z for z in rows};ids={};minwidth=1<<co.B
    for z in rows:
        i,N,box=co.encode(z['word']);ids[z['word']]=i;minwidth=min(minwidth,N)
        assert box==z['box'] and co.decode(i)[0]==z['word']
        assert co.output(i)==z['output']
    assert len(set(ids.values()))==len(rows) and minwidth>4
    edges=0
    for z in rows:
        w=z['word'];i=ids[w]
        for b in (0,1):
            raw=str(b)+w;new=raw
            for j in range(1,len(raw)+1):
                if byword[raw[:j]]['leaf']:new=raw[:j];break
            dest=byword[new]
            if z['leaf']:
                assert dest['leaf']
                assert dest['true'][0]<=bayes(z['true'][0],b,s,a)<=bayes(z['true'][1],b,s,a)<=dest['true'][1]
            else:
                x=bayes(z['center'],b,s,a)
                assert dest['true'][0]<=x<=dest['true'][1]
                if not dest['leaf']:assert x==dest['center']
            ni=co.advance(i,b);assert ni==ids[new] and co.output(ni)==dest['output'];edges+=1
    # A deterministic sample of invalid IDs, without any cache in the kernel.
    valid=set(ids.values());rng=random.Random(314159+p);invalid=0
    for _ in range(12):
        i=rng.randrange(1<<co.B)
        try:w,_=co.decode(i);assert i in valid and ids[w]==i
        except ValueError:assert i not in valid;invalid+=1
    ratio=max(z['tv']/eps for z in rows)
    return {'p':p,'true_s':str(s),'true_a':str(a),'stored_s':str(sh),'stored_a':str(ah),
            'parameter_error':str(delta),'parameter_fraction_bits':m.k,'matrix_fraction_bits':m.q,
            'coefficient_max_integer_bits':max(v.bit_length() for mat in m.mats for v in mat),
            'matrix_max_integer_bits':max(v.bit_length() for z in rows for v in (*z['box'].lo,*z['box'].hi)),
            'matrix_temporary_cap_bits':2*m.q+2*m.k+4,
            'states':len(rows),'transitions':edges,'outputs':len(rows)+edges,'bits':co.B,
            'minimum_index_bits':(len(rows)-1).bit_length(),'minimum_integer_interval':minwidth,
            'max_depth':max(map(len,byword)),'max_TV_over_epsilon_exact':str(ratio),
            'max_TV_over_epsilon_decimal':float(ratio),'sampled_unused_rejected':invalid,
            'oracle_certificate':co.certificate,'preparation_seconds':prep,'total_seconds':time.perf_counter()-start}

def main():
    rows=[]
    for p in [8,12,16]:
        for s in [F(8,25),F(1,3),F(17,50)]:
            for a in [F(37,50),F(3,4),F(19,25)]:
                row=verify_one(p,s,a);rows.append(row)
                print(json.dumps({k:row[k] for k in ['p','true_s','true_a','states','bits','max_TV_over_epsilon_decimal','total_seconds']}),flush=True)
                out={'passed':True,'configs':rows,'total_states':sum(x['states'] for x in rows),
                     'total_transitions':sum(x['transitions'] for x in rows),'total_outputs':sum(x['outputs'] for x in rows),
                     'scope':'All states and both outgoing transitions of each rounded stopping tree; no cache in the tested codec.'}
                (ROOT/'results/VERIFIED.json').write_text(json.dumps(out,indent=2)+'\n')
    print('ALL 27 GRAPHS PASSED',flush=True)
if __name__=='__main__':main()
