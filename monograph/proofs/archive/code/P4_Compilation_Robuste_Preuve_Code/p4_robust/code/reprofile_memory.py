"""Reprofile existing certified configurations after a full-depth warm-up.
This verifier-only loader reuses the compiled certificates already recorded
by check_stream_memory.py. It does not add a certificate-bypass API to Codec.
"""
from pathlib import Path
from fractions import Fraction as F
import json,tracemalloc,gc
from model import Model
from codec import Codec
from check_stream_memory import size
root=Path(__file__).resolve().parents[1];path=root/'results/STREAM_MEMORY.json';data=json.loads(path.read_text())
for row in data['rows']:
    m=Model.from_rationals(row['p'],F(row['s']),F(row['a']))
    co=Codec.__new__(Codec);co.model=m;co.sigma=F(row['certificate']['sigma'])
    co.n=row['horizon'];co.B=row['bits'];co.certificate=row['certificate']
    assert co.B==-((-co.sigma.numerator*m.p)//co.sigma.denominator)+8
    assert F(co.certificate['log_residual_bound'])<F(1,32*m.p)
    state=0;seed=0x123456789abcdef
    for j in range(80):
        seed^=seed<<13;seed&=(1<<64)-1;seed^=seed>>7;seed^=seed<<17;seed&=(1<<64)-1
        state=co.advance(state,seed&1)
    _,z=co.decode(state);assert m.leaf(z);depth=z.depth;del z
    gc.collect();tracemalloc.start()
    for j in range(12):state=co.advance(state,j%2);out=co.output(state)
    _,peak=tracemalloc.get_traced_memory();tracemalloc.stop()
    row['compiled_model_python_bytes']=size(co);row['online_tracemalloc_peak_bytes']=peak
    row['model_plus_online_peak_bytes']=size(co)+peak
    row['memory_profile_warmup_steps']=80;row['memory_profile_measured_steps']=12
    row['profile_initial_leaf_depth']=depth
    path.write_text(json.dumps(data,indent=2)+'\n')
    print({k:row[k] for k in ['p','profile_initial_leaf_depth','model_plus_online_peak_bytes']},flush=True)
