"""Longer histories, model immutability and warm-state Python memory profiles."""
from fractions import Fraction as F
from pathlib import Path
import json,random,time,tracemalloc,statistics,sys
from certificate import certify
from cases import CASES
from model import Model
from codec import Codec
from verify import belief,value,output,run_case
ROOT=Path(__file__).resolve().parents[1]

def deep_size(obj,seen=None):
    if seen is None:seen=set()
    if id(obj) in seen:return 0
    seen.add(id(obj));n=sys.getsizeof(obj)
    if isinstance(obj,dict):n+=sum(deep_size(k,seen)+deep_size(v,seen) for k,v in obj.items())
    elif isinstance(obj,(list,tuple)):n+=sum(deep_size(x,seen) for x in obj)
    elif isinstance(obj,F):n+=deep_size(obj.numerator,seen)+deep_size(obj.denominator,seen)
    elif hasattr(obj,'__dict__'):n+=deep_size(vars(obj),seen)
    return n

def signature(codec):
    return repr((vars(codec),vars(codec.model)))

def main():
    rows=[]
    for index,p in [(1,16),(3,24),(4,32)]:
        name,ms,J=CASES[index];cert=certify(ms,J=J);m=Model.from_rationals(p,cert)
        start=time.perf_counter();co=Codec(m);prep=time.perf_counter()-start
        before=signature(co);rng=random.Random(527+p);state=0;x=F(1,2);worst=F(0);timings=[]
        for _ in range(64):
            bit=rng.randrange(2);x=value(belief(ms[bit]),x)
            t=time.perf_counter();state=co.advance(state,bit);prediction=F(co.output(state),m.Q);timings.append(time.perf_counter()-t)
            tv=abs(prediction-output(ms,x));assert tv<=F(1,1<<p);worst=max(worst,tv*(1<<p))
        assert before==signature(co)
        word,z=co.decode(state);assert m.leaf(z)
        compiled=deep_size(co)
        tracemalloc.start()
        for bit in [0,1,1,0,0,1,0,1]:state=co.advance(state,bit);co.output(state)
        _,peak=tracemalloc.get_traced_memory();tracemalloc.stop()
        assert before==signature(co)
        row={'name':name,'p':p,'history_steps':64,'bits':co.B,'k':m.k,'q':m.q,'oracle_horizon':co.n,
             'max_TV_over_epsilon':str(worst),'max_TV_over_epsilon_decimal':float(worst),
             'initial_profile_leaf_depth':len(word),'compiled_python_bytes':compiled,
             'online_python_peak_bytes':peak,'model_plus_online_peak_bytes':compiled+peak,
             'memory_profile_steps':8,'median_update_and_output_seconds':statistics.median(timings),
             'preparation_seconds':prep,'model_unchanged':True,'oracle':co.certificate}
        rows.append(row);print(json.dumps({k:row[k] for k in ['name','p','bits','oracle_horizon','max_TV_over_epsilon_decimal','model_plus_online_peak_bytes']}),flush=True)
        (ROOT/'results/STREAMS.json').write_text(json.dumps({'passed':len(rows)==3,'rows':rows,
            'scope':'Finite histories and warm-state tracemalloc allocations, not an exhaustive peak or interpreter RSS.'},indent=2)+'\n')
    # A general positive instrument without a shared hidden transition matrix.
    name,ms,J=CASES[1];changed=[list(m) for m in ms]
    changed[0][0]+=F(1,1000);changed[0][2]-=F(1,1000)
    assert changed[0][0]/(changed[0][0]+changed[0][2])!=changed[1][0]/(changed[1][0]+changed[1][2])
    row=run_case(('instrument_E',tuple(tuple(m) for m in changed),J),8)
    (ROOT/'results/GENERAL_INSTRUMENT.json').write_text(json.dumps({'passed':True,'config':row,'shared_hidden_transition':False},indent=2)+'\n')
    print('STREAMS AND GENERAL INSTRUMENT PASSED',flush=True)
if __name__=='__main__':main()
