from fractions import Fraction as F
from pathlib import Path
import json
q=F(16,49)
assert 1/(2*(1-q))<F(3,4)
assert F(32,21)<2 # log(2)>1/2 implies the weight Lipschitz bound 2
assert F(12103,10692)**6>2 # true suffix cost < p-1/6
assert F(16,49)<F(1,2) and F(16,81)>F(1,8)
assert F(247,64)*4*F(16,81)*F(49,132)==F(12103,10692)
assert F(27)**6<2**29
assert F(50,3)**6<2**25
assert F(132,49)**3>16
assert F(29,6)+1+F(3,5)*(F(1,2)+F(25,6)-F(4,3))<8
assert F(16,49)**20/4<=F(49,33*2**30)
assert 7**20<2**64
assert max(49*(16*7**20)**2,528*9**20*2**30)<2**128
for p in range(4,10001):
    m=8*(p+1);R=m*p+(m+5)//6
    assert F(p)-F(1,6)+F(5*p,2*m)<F(R,m)
    assert F(R,m)+F(5*p,2*m)<p+F(1,2)
    assert R<m*(p+1)
out={'passed':True,'rational_inequalities':'verified exactly',
     'budget_inequalities_sampled_p':[4,10000],
     'note':'The report proves the budget inequalities for all p; finite loop is an extra check.',
     'theorem':'B <= d log2(1/epsilon_p)+8',
     'scope':'fixed central model; polynomial preprocessing and immutable count data are charged'}
(Path(__file__).resolve().parents[1]/'results/PROOF_CHECKS.json').write_text(json.dumps(out,indent=2)+'\n')
print(out)
