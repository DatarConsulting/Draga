"""Long products against an exact rational oracle, used only in tests."""
from fractions import Fraction as F
from pathlib import Path
import random,json,time,math
from model import Model

def mult(z,m):
    a,b,c,d=z;A,B,C,D=m
    return a*A+b*C,a*B+b*D,c*A+d*C,c*B+d*D
def project(z,x,y):
    a,b,c,d=z
    return F(x*a+(y-x)*b,x*(a+c)+(y-x)*(b+d))
def main():
    start=time.perf_counter();rows=[]
    for p in [4,8,16,32,64,128,256]:
        for s,a in [(F(8,25),F(37,50)),(F(1,3),F(3,4)),(F(17,50),F(19,25))]:
            m=Model.from_rationals(p,s,a);z=m.root();ex=(1,0,0,1);reduced=ex;rng=random.Random(p)
            maxq=0;pointchecks=0
            for n in range(1,p+2):
                b=rng.randrange(2);z=m.step(z,b);ex=mult(ex,m.mats[b])
                reduced=mult(reduced,m.mats[b]);common=math.gcd(*reduced)
                reduced=tuple(v//common for v in reduced)
                maxq=max(maxq,max(v.bit_length() for v in (*z.lo,*z.hi)))
                for x,y in [(3,8),(1,2),(5,8)]:
                    exact=project(ex,x,y);lo,hi=m.point(z,x,y)
                    assert F(lo,m.Q)<=exact<=F(hi,m.Q)
                    assert max(exact-F(lo,m.Q),F(hi,m.Q)-exact)<=F(1568*n+1,m.Q)
                    pointchecks+=1
                if n>=p:assert m.leaf(z)
            row={'p':p,'s':str(s),'a':str(a),'k':m.k,'q':m.q,'tested_depth':p+1,
                 'stored_matrix_max_bits':maxq,'matrix_intermediate_cap_bits':2*m.q+2*m.k+4,
                 'exact_reference_max_bits':max(x.bit_length() for x in ex),
                 'gcd_reduced_reference_max_bits':max(x.bit_length() for x in reduced),'point_checks':pointchecks}
            assert project(reduced,1,2)==project(ex,1,2)
            rows.append(row)
        print('products p',p,'passed',flush=True)
    out={'passed':True,'rows':rows,'total_exact_point_checks':sum(x['point_checks'] for x in rows),
         'note':'Exact reference products grow quadratically and are test-only; operational products remain capped.',
         'seconds':time.perf_counter()-start}
    (Path(__file__).resolve().parents[1]/'results/PRODUCTS.json').write_text(json.dumps(out,indent=2)+'\n')
    print('PASSED',out['total_exact_point_checks'],flush=True)
if __name__=='__main__':main()
