from pathlib import Path
import json,time
from uniform import Codec,RollingCodec,compile_lib
ROOT=Path(__file__).resolve().parents[1]
lib=compile_lib();rows=[]
for p,words in [(4,['','0','1']),(8,['','010','1010']),(23,['','010101010101','001010101010'])]:
    start=time.perf_counter();r=RollingCodec(p);c=Codec(p,'c',lib)
    assert r.size==c.size
    for word in words:
        index=r.encode(word)
        assert index==c.encode(word) and r.decode(index)==word
    row={'p':p,'words_checked':words,'codes_match':True,
         'rolling_counter_cells':(max(r.cost)+1)*r.J,'elapsed_seconds':time.perf_counter()-start}
    rows.append(row);print(row,flush=True)
(ROOT/'results/ROLLING.json').write_text(json.dumps({'passed':True,'configs':rows},indent=2)+'\n')
