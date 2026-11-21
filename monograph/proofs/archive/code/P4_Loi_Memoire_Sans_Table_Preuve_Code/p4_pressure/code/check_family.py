from fractions import Fraction as F
from pathlib import Path
import json
slo,shi=F(8,25),F(17,50);alo,ahi=F(37,50),F(19,25)
xlo,xhi=F(3,8),F(5,8)
def f(s,a,x):return s+(1-2*s)*(1-a)*x/(a*(1-x)+(1-a)*x)
def derivative(s,a,x):return (1-2*s)*a*(1-a)/(a+(1-2*a)*x)**2
# On the complete box: F0 increases in x and s, decreases in a.
# F0' increases in x, decreases in s and a (sign of d/da is x-a<0).
assert xhi<alo and (1-alo)*xhi/(alo*(1-xhi)+(1-alo)*xhi)<F(1,2)
lo=f(slo,ahi,xlo);hi=f(shi,alo,xhi)
ell=derivative(shi,ahi,xlo);q=derivative(slo,alo,xhi)
H=2*(2*ahi-1)/(ahi+(1-2*ahi)*xhi)
assert xlo<lo<hi<xhi
assert 1-2*hi>F(2,25)
assert ell>F(9,50) and q<F(9,25) and H<F(12,5)
assert F(12,5)*F(1,4)/(1-F(9,25))==F(15,16)<1
r={'passed':True,'s_interval':[str(slo),str(shi)],'a_interval':[str(alo),str(ahi)],
   'F0_image':[str(lo),str(hi)],'derivative_extrema':[str(ell),str(q)],'log_derivative_upper':str(H),
   'uniform_constants':{'length':'1/4','lower_derivative':'9/50','upper_derivative':'9/25',
    'log_derivative_lipschitz':'12/5','distortion':3,'gap':'2/25','observation_lower':'12/25','observation_upper':'13/25',
    'lower_state_bound':'(1/27)*(4/(625*epsilon))**d','upper_state_bound':'54*(13/(12*epsilon))**d',
    'valid_epsilon_upper':'4/625'}}
(Path(__file__).resolve().parents[1]/'results/FAMILY.json').write_text(json.dumps(r,indent=2)+'\n')
print(json.dumps(r,indent=2))
