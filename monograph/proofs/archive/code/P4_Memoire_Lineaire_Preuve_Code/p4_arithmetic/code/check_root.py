"""Replay of the executed root-refinement test (roughly 40 s locally)."""
import time,json
from fractions import Fraction as F
from pathlib import Path
from universal import refine_s
start=time.perf_counter()
s=refine_s(16,F(1,2),F(51,100))
assert max(abs(s-F(50501,100000)),abs(s-F(25251,50000)))<=F(1,256)
out={'passed':True,'p':16,'test_initial_bracket':['1/2','51/100'],
     'returned_s':str(s),'allowed_error':'1/256','seconds_local':time.perf_counter()-start,
     'note':'Exercises the root interval/refinement branch with a wider valid starting bracket. Not a large-p performance validation.'}
(Path(__file__).resolve().parents[1]/'results/ROOT_REFINEMENT.json').write_text(json.dumps(out,indent=2)+'\n')
print(out)
