"""High-precision independent numerical sanity checks, NOT formal certificates.
Core computations remain integer intervals; mpmath is used only in this test.
"""
from pathlib import Path
from fractions import Fraction as F
import json
import mpmath as mp
from certificate import certify
from cases import CASES
from model import Model
from oracle import ratio
ROOT=Path(__file__).resolve().parents[1]
mp.mp.dps=90
def real(x):x=F(x);return mp.mpf(x.numerator)/x.denominator
def V(model,t,n,s):
    if n==0:return mp.mpf(1)
    result=mp.mpf(0)
    for mat,de in zip(model.mats,model.dets):
        a,b,c,d=map(mp.mpf,mat);den=(a+b)*t+(c+d)*(1-t)
        g=mp.mpf(de)/den**2;nextt=(a*t+c*(1-t))/den
        result+=mp.power(g,s)*V(model,nextt,n-1,s)
    return result

def main():
    configs=json.loads((ROOT/'results/GRAPHS.json').read_text())['configs'];checks=0
    for row in configs:
        case=next(x for x in CASES if x[0]==row['name']);cert=certify(case[1],J=case[2]);model=Model.from_rationals(row['p'],cert)
        n=row['oracle']['horizon'];s=F(row['oracle']['sigma'])
        for j in range(5):
            t=cert.Jdual[0]+F(j,4)*(cert.Jdual[1]-cert.Jdual[0])
            interval,c=ratio(model,t,n,s,F(1,1<<36))
            ref=V(model,real(t),n+1,real(s))/V(model,real(t),n,real(s))
            assert mp.mpf(interval[0])/c.S<=ref<=mp.mpf(interval[1])/c.S
            checks+=1
    result={'passed':True,'checks':checks,'decimal_digits':90,
            'scope':'Independent floating high-precision sanity test only; rigorous enclosures follow the interval proof.'}
    (ROOT/'results/ORACLE_REFERENCE.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)
if __name__=='__main__':main()
