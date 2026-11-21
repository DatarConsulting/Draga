"""Exact rational verification of all constants in the accompanying proof."""
from fractions import Fraction as F
from pathlib import Path
import json
slo,shi=F(8,25),F(17,50);alo,ahi=F(37,50),F(19,25)
zlo=(1+6*slo)/8;zhi=1-zlo
Dlo=ahi+(1-2*ahi)*zhi;Dhi=ahi+(1-2*ahi)*zlo
Glo=(1-ahi)*zlo/Dhi;Ghi=(1-alo)*zhi/(alo+(1-2*alo)*zhi)
glo=(1-2*shi)*ahi*(1-ahi)/Dhi**2
ghi=(1-2*slo)*alo*(1-alo)/Dlo**2
Hlog=2*(2*ahi-1)*(1-2*slo)/Dlo
assert F(1,8)<Glo<Ghi<F(7,8)
assert glo>F(7,40) and ghi<F(3,8) and Hlog<F(7,8)
assert F(7,40)>F(1,8) and F(3,8)**3<F(1,16)
assert F(7,40)**3>F(1,256)
K=F(9,8);q=F(3,8);D=F(3,4);A=F(3,4)*F(7,8)
assert A+q*K==F(69,64)<K
assert q+D*K/2==F(51,64)<F(4,5)
assert F(27,32)*2/F(4096)+F(1,4096)<F(3,4096)
assert F(1,64)+F(27,2048)<F(1,32)
assert F(3,40)<(1-ahi)*slo
assert ahi*(1-slo)<F(21,40)
assert F(21,40)/F(3,40)==7
assert F(1571,1<<20)<F(1,512)
assert F(18,25)+F(1,256)<2
assert F(1,4)+F(1,256)<2
assert (2-F(1,256))*F(3,50)*F(247,64)>F(1,4)
assert F(15,7)*3<7
assert F(3,4096)+F(9,16)*F(1,2048)<F(1,512)
assert F(27,64)+F(3,32)+F(1,32)+F(1,32)==F(37,64)<F(2,3)
assert 32*F(3,4)-20>0
assert F(4,3)**16>2
Fa=F(9,25)*F(1,4)/F(87,200)**2
modelTV=1+F(13,25)*(1+Fa)/(1-F(9,25))
assert modelTV==F(29589,13456)<F(9,4)
assert F(13,25)+F(1,256)+F(9,4*4096)<1
assert F(25,8)+F(25,6)<8
assert F(25,4)+(2*ahi-1)/(ahi*(1-ahi))<10
assert 26*F(3,4)<20
assert F(9)+F(1,16)+F(10,4096)<10
for p in range(4,10001):
    qbits=p+p.bit_length()+20;Q=1<<qbits
    assert F(392*(p+1),Q)<F(1,2)
    assert F(1568*(p+1)+3,Q)<F(1,512*(1<<p))
    assert p*4<=1<<p
out={'passed':True,'scope':'Exact checks of proof constants, with finite extra loop p=4..10000; analytic arguments are in the report.',
     'dual_derivative_bounds':['7/40','3/8'],'dual_log_derivative_bound':'7/8',
     'cone_K':'9/8','contraction_theta':'4/5','pressure_root_interval':['1/3','3/4'],
     'parameter_TV_bound':'9/4 * delta','pressure_root_stability_bound':'40 * delta',
     'matrix_relative_error':'392*n/2^q','point_error':'(1568*n+3)/2^q',
     'precision_q':'p+ceil(log2(p+1))+20','parameter_k':'p+12',
     'dyadic_precision_theorem':'B <= d_true*p+10','arbitrary_precision_theorem':'B <= d_true*log2(1/epsilon)+11'}
(Path(__file__).resolve().parents[1]/'results/PROOF_CHECKS.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
