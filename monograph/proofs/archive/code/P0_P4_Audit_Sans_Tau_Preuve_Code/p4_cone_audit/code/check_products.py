"""Independent exact projective products, and relative-error induction replay."""
from fractions import Fraction as F
from pathlib import Path
from math import gcd
from functools import reduce
import json,random,time
from certificate import certify,point
from cases import CASES
from model import Model
ROOT=Path(__file__).resolve().parents[1]

def mul(x,y):
    a,b,c,d=x;A,B,C,D=y
    return a*A+b*C,a*B+b*D,c*A+d*C,c*B+d*D

def main():
    started=time.perf_counter();rows=[];relative_checks=0
    for case in CASES:
        name,ms,J=case;c=certify(ms,J=J)
        for p in [4,8,16,32,64,128,256]:
            m=Model.from_rationals(p,c);z=m.root();exact=(1,0,0,1);X=tuple(map(F,exact))
            rng=random.Random(333+p);maxbits=0;refbits=0;checks=0
            for depth in range(1,m.depth_cap+2):
                b=rng.randrange(2)
                if p<=8:
                    W=sum(mul(z.hi,m.mats[b]));X=tuple(v*m.Q/W for v in mul(X,m.mats[b]))
                z=m.step(z,b);exact=mul(exact,m.mats[b]);g=reduce(gcd,exact);exact=tuple(v//g for v in exact)
                maxbits=max(maxbits,max(v.bit_length() for v in (*z.lo,*z.hi)))
                refbits=max(refbits,max(v.bit_length() for v in exact))
                for x in [J[0],F(1,2),J[1]]:
                    lo,hi=m.point(z,x);v=point(exact,x)
                    assert F(lo,m.Q)<=v<=F(hi,m.Q)
                    assert max(v-F(lo,m.Q),F(hi,m.Q)-v)<=(32*c.R**2*depth+3)/m.Q
                    checks+=1
                if p<=8:
                    eta=8*c.R**2*depth/m.Q
                    for low,xx,high in zip(z.lo,X,z.hi):
                        assert (1-eta)*xx<=F(low,m.Q)<=xx<=F(high,m.Q)<=(1+eta)*xx
                        relative_checks+=1
            eps=F(1,1<<p)
            assert c.rho**m.depth_cap*(J[1]-J[0])<=eps/4
            assert c.rho*eps/2+2*m.point_error_bound<eps/2
            assert 2**c.offset_bits*c.mass_constant>8+4/c.beta
            assert c.K>=c.sigma_hi*c.Ddual+c.grho*c.K
            assert max(c.parameter_TV_constant,c.pressure_stability_constant)*m.parameter_error_bound<=eps/8
            rows.append({'name':name,'p':p,'q':m.q,'k':m.k,'depth':m.depth_cap+1,'point_checks':checks,
                         'stored_max_bits':maxbits,'exact_reduced_max_bits':refbits,'temporary_cap_bits':m.intermediate_cap})
    result={'passed':True,'rows':rows,'exact_point_checks':sum(r['point_checks'] for r in rows),
            'relative_error_component_checks':relative_checks,'seconds':time.perf_counter()-started,
            'scope':'Exact large products and exact scales are verifier-only; the filter stores neither.'}
    (ROOT/'results/PRODUCTS.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='rows'}),flush=True)
if __name__=='__main__':main()
