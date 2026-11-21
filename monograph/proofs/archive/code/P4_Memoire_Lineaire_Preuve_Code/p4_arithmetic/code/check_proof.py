from fractions import Fraction as F
from pathlib import Path
import json
q=F(16,49);K=F(3,5);theta=K;D=F(3,4)
assert F(4,1)*F(51,100)/F(21,4)<F(2,5)
assert F(2,5)+q*K<K
assert q+D*K/2<theta
assert F(9,7)**4<F(3,2)**3
assert F(9,7)**4*F(33,32)**3<F(29,20)**3
assert F(1,8)+F(9,320)<F(1,4)
assert F(9,40)+F(1,16)+F(1,4)+F(1,32)==F(91,160)<F(2,3)
assert F(16,81)**3>F(3,8)**5
assert F(8,7)<F(5,4)
assert F(12103,10692)>1
assert 2*(F(1,2048)+F(5,256))<F(1,16)
assert (1-F(30,3200))*(1-F(30,100000))>F(7,8)
assert q**20/4<=F(49,33*(1<<30))
assert 7**20<1<<64
assert max(49*(16*7**20)**2,528*9**20*2**30)<1<<128
assert 2*4*7**20<1<<60
out={'passed':True,'method':'Exact rational checks of constants used in the written proof; not a formal proof assistant.',
     'cone_K':'3/5','contraction_theta':'3/5','log_mass_loss_bound':'91/160',
     'universal_bits':'ceil(s*p)+5 <= d*L_p+6',
     'finite_bits':'ceil(101003*p/200000)+3, 4<=p<=30',
     'finite_rounding_loss_bound':5,'universal_rounding_loss_bound':10}
(Path(__file__).resolve().parents[1]/'results/PROOF_CHECKS.json').write_text(json.dumps(out,indent=2)+'\n')
print(out)
