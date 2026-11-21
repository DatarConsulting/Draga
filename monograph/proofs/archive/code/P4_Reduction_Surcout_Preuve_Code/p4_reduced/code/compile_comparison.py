from reduced_codec import ReducedCodec
from certificate import certify,hmm
from model import Model
from direct_grid import DirectGrid
from fractions import Fraction as F
from pathlib import Path
import json,time

ROOT=Path(__file__).resolve().parents[1]
CASES=[('S',hmm('1/4','1/4','6/25','19/25')),
       ('A',hmm('6/25','13/50','6/25','3/4'))]

def run():
    rows=[]
    for name,ms in CASES:
        H=certify(ms,J=(F(1,4),F(3,4)))
        for p in (4,8,12,16,24,64,128,256):
            t=time.perf_counter();m=Model.from_rationals(p,H);c=ReducedCodec(m)
            row=dict(model=name,p=p,old_bits=c.old_B,new_bits=c.B,
                     new_offset=c.capacity_certificate['offset_bits'],grid_bits=DirectGrid(m).B,
                     sigma=str(c.sigma),oracle_horizon=c.n,seconds=time.perf_counter()-t,
                     capacity=c.capacity_certificate)
            rows.append(row)
            print(json.dumps({k:row[k] for k in row if k!='capacity'}),flush=True)
            (ROOT/'results/COMPARISON.json').write_text(json.dumps(dict(passed=len(rows)==16,rows=rows),indent=2)+'\n')
if __name__=='__main__':run()
