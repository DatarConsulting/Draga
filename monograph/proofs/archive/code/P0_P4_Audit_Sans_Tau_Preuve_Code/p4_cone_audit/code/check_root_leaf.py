"""Admissible narrow invariant interval: the root itself is a stopping leaf."""
from fractions import Fraction as F
from pathlib import Path
import json
from certificate import certify
from verify import run_case
ROOT=Path(__file__).resolve().parents[1]
def main():
    mats=[]
    for u,c in [(F(3,4),F(1,2)-F(1,5000)),(F(1,4),F(1,2)+F(1,5000))]:
        v=1-u;A=u*c+F(1,40);B=v*c-F(1,40)
        mats.append((A,B,u-A,v-B))
    J=(F(999,2000),F(1001,2000));ms=tuple(mats)
    cert=certify(ms,J=J,radius=F(1,10**8))
    row=run_case(('root_leaf',ms,J),4,cert)
    assert row['states']==1 and row['transitions']==2 and row['max_depth']==0
    (ROOT/'results/ROOT_LEAF.json').write_text(json.dumps({'passed':True,'config':row},indent=2)+'\n')
    print('ROOT LEAF PASSED',flush=True)
if __name__=='__main__':main()
