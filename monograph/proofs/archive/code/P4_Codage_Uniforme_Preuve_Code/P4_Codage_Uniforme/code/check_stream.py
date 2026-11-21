from fractions import Fraction as F
from pathlib import Path
import random,json,time
from uniform import Codec
from verify import bayes
ROOT=Path(__file__).resolve().parents[1]
rows=[]
for p in [4,23,31]:
    start=time.perf_counter();c=Codec(p,'python');x=F(1,2);state=0
    eps=F(49,132*(1<<p));worst=F(0);rng=random.Random(1823+p)
    for i in range(600):
        b=(i&1) if i<200 else 0 if i<400 else rng.randrange(2)
        x=bayes(x,b);state,out=c.advance(state,b)
        error=abs(out-(F(3,4)-x/2));assert error<=eps
        worst=max(worst,error/eps)
    row={'p':p,'steps':600,'coding_words':c.size,'coding_bits':c.bits,
         'max_error_over_epsilon_display':float(worst),'exact_checks_passed':True,'seconds':time.perf_counter()-start}
    rows.append(row);print(row,flush=True)
(ROOT/'results/STREAM.json').write_text(json.dumps({'passed':True,'configs':rows},indent=2)+'\n')
