"""Finite streams at precisions beyond old C limit, with Python memory audit."""
from fractions import Fraction as F
from pathlib import Path
import json,time,sys,tracemalloc,statistics,gc
from model import Model
from codec import Codec
from verify import bayes,predict

def size(obj,seen=None):
    if seen is None:seen=set()
    if id(obj) in seen:return 0
    seen.add(id(obj));z=sys.getsizeof(obj)
    if isinstance(obj,dict):return z+sum(size(k,seen)+size(v,seen) for k,v in obj.items())
    if isinstance(obj,(tuple,list)):return z+sum(size(x,seen) for x in obj)
    if isinstance(obj,F):return z+size(obj.numerator,seen)+size(obj.denominator,seen)
    if hasattr(obj,'__dict__'):return z+size(vars(obj),seen)
    return z

def main():
    root=Path(__file__).resolve().parents[1];rows=[]
    for p,s,a in [(8,F(8,25),F(37,50)),(16,F(1,3),F(3,4)),
                   (32,F(17,50),F(19,25)),(64,F(1,3),F(3,4))]:
        gc.collect();tracemalloc.start();start=time.perf_counter()
        m=Model.from_rationals(p,s,a);co=Codec(m)
        seconds=time.perf_counter()-start;_,prep_peak=tracemalloc.get_traced_memory();tracemalloc.stop()
        model_bytes=size(co);before=repr(vars(m))+repr(vars(co));state=0
        times=[];tvmax=F(0);x=F(1,2);state=0;seed=0x123456789abcdef
        for j in range(80):
            seed^=seed<<13;seed&=(1<<64)-1;seed^=seed>>7;seed^=seed<<17;seed&=(1<<64)-1;b=seed&1
            start=time.perf_counter();state=co.advance(state,b);out=co.output(state);times.append(time.perf_counter()-start)
            x=bayes(x,b,s,a);tv=abs(F(out,m.Q)-predict(x,a));assert tv<=F(1,1<<p);tvmax=max(tvmax,tv)
        # Measure after the 80-step warm-up, so the ID represents a full leaf.
        # Reference fractions and verifier histories are outside this phase.
        assert m.leaf(co.decode(state)[1])
        tracemalloc.start()
        for j in range(12):state=co.advance(state,j%2);out=co.output(state)
        _,work_peak=tracemalloc.get_traced_memory();tracemalloc.stop()
        assert before==repr(vars(m))+repr(vars(co)) # no history-dependent cache/counter
        row={'p':p,'s':str(s),'a':str(a),'model_k':m.k,'matrix_q':m.q,'bits':co.B,'horizon':co.n,
             'stream_steps':80,'max_TV_over_epsilon':str(tvmax*(1<<p)),
             'max_TV_over_epsilon_decimal':float(tvmax*(1<<p)),
             'compiled_model_python_bytes':model_bytes,'preparation_tracemalloc_peak_bytes':prep_peak,
             'online_tracemalloc_peak_bytes':work_peak,'model_plus_online_peak_bytes':model_bytes+work_peak,
             'preparation_seconds_with_tracing':seconds,'median_advance_and_output_seconds':statistics.median(times),
             'model_unchanged_after_stream':True,'memory_profile_warmup_steps':80,
             'memory_profile_measured_steps':12,'certificate':co.certificate}
        rows.append(row);print(json.dumps({k:row[k] for k in ['p','bits','horizon','model_plus_online_peak_bytes','median_advance_and_output_seconds']}),flush=True)
        (root/'results/STREAM_MEMORY.json').write_text(json.dumps({'passed':True,'rows':rows,
             'note':'Tracemalloc counts Python allocations made in each phase, not OS/interpreter RSS. Compiled model size is separately included online. Preparation timings include tracing overhead; online timings do not.'},indent=2)+'\n')
    print('STREAMS AND MEMORY AUDIT PASSED',flush=True)
if __name__=='__main__':main()
