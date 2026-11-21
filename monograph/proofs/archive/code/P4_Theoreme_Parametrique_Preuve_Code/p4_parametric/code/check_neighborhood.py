"""Reuse one unchanged certificate for two distinct true models in its inner box."""
from fractions import Fraction as F
from pathlib import Path
import json
from certificate import certify
from cases import CASES
from verify import run_case
ROOT=Path(__file__).resolve().parents[1]
def main():
    _,base,J=CASES[1];cert=certify(base,J=J);before=cert.json();rows=[]
    for sign in (-1,1):
        ms=[list(m) for m in base];delta=sign*cert.radius/4
        ms[0][0]+=delta;ms[1][2]-=delta
        ms[0][1]-=delta;ms[1][3]+=delta
        name='box_minus' if sign<0 else 'box_plus'
        r=run_case((name,tuple(tuple(m) for m in ms),J),8,cert);rows.append(r)
        assert cert.json()==before
        print(name,r['states'],r['max_TV_over_epsilon_decimal'],flush=True)
    (ROOT/'results/NEIGHBORHOOD.json').write_text(json.dumps({'passed':True,'same_certificate_reused':True,'configs':rows},indent=2)+'\n')
if __name__=='__main__':main()
