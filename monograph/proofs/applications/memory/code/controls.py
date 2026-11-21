"""Boundary obstructions and diagnostic floating-point baselines."""
import json
from pathlib import Path
import sympy as sp
import numpy as np
from fractions import Fraction as F
ROOT=Path(__file__).resolve().parents[1]
def run():
    I=sp.eye(2);P=sp.Matrix([[0,1],[1,0]]);records=[]
    for c in [sp.Rational(3,5),sp.Rational(1),sp.Rational(4,5)]:
        M=c*sp.ones(2);inverse=(I-M.T).inv();X=I/(1-2*c*c)
        assert X-c*I*X*(c*I).T-c*P*X*(c*P).T==I
        records.append({'scale':str(c),'counting_Stein_X':[[str(x) for x in row] for row in X.tolist()],
          'counting_stable':bool(1-2*c*c>0),'comparison_inverse':[[str(x) for x in row] for row in inverse.tolist()],
          'comparison_certificate_accepts':all(x>=0 for x in inverse),'pointwise_contraction_infinity':str(c),
          'reason':'c=3/5 is counting-stable but entrywise-max resolvent rejects. For 1-2c^2<=0, the trace identity rules out positive Stein X.'})
    (ROOT/'results/boundary_controls.json').write_text(json.dumps(records,indent=2))
    rows=[]
    to_float=lambda a:np.array([[float(F(x)) for x in row] for row in a])
    for path in sorted((ROOT/'models').glob('*.json')):
        m=json.loads(path.read_text());AA=[to_float(a) for a in m['source_A']];ab=to_float(m['source_beta']);ax=to_float(m['source_alpha'])
        for basis,f in m['frames'].items():
          if not f['accepted']:continue
          for typ in [np.float16,np.float32]:
            A=[to_float(a).astype(typ) for a in f['unrounded_A']];b=to_float(f['unrounded_beta']).astype(typ);x=to_float(f['unrounded_alpha']).astype(typ);src=ax.copy()
            mx=float(np.max(abs(src@ab-(x@b).astype(float))));nodes=1;finite=True
            for depth in range(1,17):
                x=np.concatenate([x@a for a in A]);src=np.concatenate([src@a for a in AA]);y=x@b
                finite=finite and bool(np.isfinite(x).all() and np.isfinite(y).all());mx=max(mx,float(np.max(abs(src@ab-y.astype(float)))));nodes+=len(x)
            rows.append({'id':m['id'],'basis':basis,'format':np.dtype(typ).name,'persistent_state_bytes':2*np.dtype(typ).itemsize,
              'nodes_through_length16':nodes,'max_observed_abs_error':mx,'all_finite':finite})
    (ROOT/'results/float_controls.json').write_text(json.dumps({'records':rows,'scope':'NumPy matmul and typed state storage, diagnostic only; internal accumulation depends on NumPy/backend. No all-horizon FP16/FP32 proof.'},indent=2))
    print('PASS controls',len(rows),flush=True)
if __name__=='__main__':run()
