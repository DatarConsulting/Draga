from pathlib import Path
from fractions import Fraction as F
import json,random,sys
from portable import advance
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'reference'))
from oracle import bayes
rows=[]
for p in [4,19,21,23,30,64,256,1024]:
    eps=F(49,132*(1<<p));state=1;x=F(1,2)
    rng=random.Random(324123+p);steps=max(400,p+100)
    maxbits=maxcoord=maxdepth=0;maxratio=F(0)
    for i in range(steps):
        # Includes long alternating and constant histories as well as random.
        bit=(i&1) if i<steps//3 else 0 if i<2*steps//3 else rng.randrange(2)
        x=bayes(x,bit)
        state,out,stats=advance(p,state,bit)
        error=abs(out-(F(3,4)-x/2))
        assert error<=eps
        maxratio=max(maxratio,error/eps)
        maxbits=max(maxbits,state.bit_length());maxcoord=max(maxcoord,stats['coordinate_bits'])
        maxdepth=max(maxdepth,stats['processed_bits'])
    row={'p':p,'steps':steps,'max_state_bits':maxbits,'max_coordinate_bits':maxcoord,
         'max_processed_bits':maxdepth,'max_error_over_epsilon_display':float(maxratio),'exact_error_checks_passed':True}
    rows.append(row);print(row,flush=True)
(ROOT/'results/PORTABLE.json').write_text(json.dumps({'passed':True,'configs':rows,
    'note':'Exact comparisons on test histories; all-history guarantee follows from invariant-interval proof.'},indent=2)+'\n')
