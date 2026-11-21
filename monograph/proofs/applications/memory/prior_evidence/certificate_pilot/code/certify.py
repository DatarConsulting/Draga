"""Exact dyadic stability and Hankel-error upper certificates for saved heads.

Run after training: python code/certify.py
Synthetic checks only: python code/certify.py --self-test

The numerical Stein solves propose candidates. Exact rational residual bounds
and exact LDL positivity make the final saved-coefficient claims rigorous.
These certificates do not cover execution roundoff or establish AAK optimality.
"""
import os
for _key in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS']:
    os.environ[_key]='1'
import argparse, hashlib, json, math, time
from fractions import Fraction as F
from pathlib import Path
import numpy as np
import torch
import benchmark

ROOT=Path(__file__).resolve().parents[1]

class Dyadic:
    """Exact ndarray of integers divided by one power of two."""
    def __init__(self, nums, exponent=0):
        self.n=np.asarray(nums,dtype=object); self.e=int(exponent)
    @classmethod
    def floats(cls,a):
        a=np.asarray(a,dtype=np.float64)
        if not np.isfinite(a).all(): raise ValueError('nonfinite coefficient')
        ratios=[float(x).as_integer_ratio() for x in a.flat]
        powers=[d.bit_length()-1 for _,d in ratios]
        e=max(powers,default=0)
        return cls(np.array([n << (e-p) for (n,_),p in zip(ratios,powers)],dtype=object).reshape(a.shape),e)
    @property
    def T(self): return Dyadic(self.n.T,self.e)
    def __add__(self,other):
        e=max(self.e,other.e)
        return Dyadic(self.n*(1<<(e-self.e))+other.n*(1<<(e-other.e)),e)
    def __sub__(self,other):
        e=max(self.e,other.e)
        return Dyadic(self.n*(1<<(e-self.e))-other.n*(1<<(e-other.e)),e)
    def __matmul__(self,other): return Dyadic(self.n@other.n,self.e+other.e)
    def scalar(self): return F(int(self.n.item()),1<<self.e)
    def abs_sum(self): return F(sum(abs(int(x)) for x in self.n.flat),1<<self.e)
    def fractions(self):
        return [[F(int(x),1<<self.e) for x in row] for row in self.n]
    def export(self):
        return {'shape':list(self.n.shape),'denominator_power_of_two':self.e,
                'numerators':np.array([str(int(x)) for x in self.n.flat],dtype=object).reshape(self.n.shape).tolist()}

def rational(q): return {'numerator':str(q.numerator),'denominator':str(q.denominator)}

def float_upper(q):
    value=float(q)
    if F.from_float(value)<q: value=float(np.nextafter(value,np.inf))
    assert F.from_float(value)>=q
    return value

def sqrt_upper(q,decimal_places=12):
    if q<0: raise ValueError('negative square bound')
    scale=10**decimal_places; n=q.numerator*scale*scale; d=q.denominator
    k=math.isqrt(n//d)
    if k*k*d<n: k+=1
    eta=F(k,scale)
    assert eta*eta>=q
    return eta

def ldl_positive(m):
    """Exact no-pivot LDL; positivity of all pivots proves SPD."""
    n=len(m); l=[[F(0) for _ in range(n)] for _ in range(n)]; ds=[]
    for i in range(n):
        l[i][i]=F(1)
        d=m[i][i]-sum(l[i][j]*l[i][j]*ds[j] for j in range(i))
        if d<=0: return False,ds+[d]
        ds.append(d)
        for k in range(i+1,n):
            l[k][i]=(m[k][i]-sum(l[k][j]*l[i][j]*ds[j] for j in range(i)))/d
    return True,ds

def stable(r,q):
    r=np.asarray(r); n=r.shape[-1]
    rd=[Dyadic.floats(a) for a in r]
    g=Dyadic(np.zeros((n,n),dtype=object))
    for a in rd: g=g+a.T@a
    m=g.fractions()
    m=[[ (q if i==j else F(0))-m[i][j] for j in range(n)] for i in range(n)]
    ok,pivots=ldl_positive(m)
    return {'passed':ok,'q':rational(q),'ldl_pivots':[rational(x) for x in pivots],
            'minimum_pivot_float':float(min(pivots)) if pivots else None,
            'numeric_sum_RtR_largest_eigenvalue':float(np.linalg.eigvalsh(sum(a.T@a for a in r)).max()) if n else 0.}

def apply_exact(e,x):
    y=Dyadic(np.zeros_like(x.n,dtype=object))
    for a in e: y=y+a@x@a.T
    return y

def error_certificate(r,b,c,k,q):
    r=np.asarray(r,dtype=np.float64); b=np.asarray(b,dtype=np.float64); c=np.asarray(c,dtype=np.float64)
    n=len(b); assert 0<=k<=n
    full_stability=stable(r,q); small_stability=stable(r[:,:k,:k],q)
    if not full_stability['passed'] or not small_stability['passed']:
        raise AssertionError('exact sufficient stability check failed')
    e=np.zeros((len(r),n+k,n+k),dtype=np.float64)
    e[:,:n,:n]=r; e[:,n:,n:]=r[:,:k,:k]
    B=np.concatenate([b,b[:k]]); C=np.concatenate([c,-c[:k]])
    cp=sum(np.kron(a,a) for a in e)
    system=np.eye((n+k)**2)-cp
    g=np.linalg.solve(system,np.outer(B,B).reshape(-1)).reshape(n+k,n+k)
    g=(g+g.T)/2
    s=np.linalg.solve(system,g.reshape(-1)).reshape(n+k,n+k)
    s=(s+s.T)/2
    assert np.array_equal(g,g.T) and np.array_equal(s,s.T)
    ed=[Dyadic.floats(a) for a in e]
    bd=Dyadic.floats(B[:,None]); cd=Dyadic.floats(C[:,None])
    gd=Dyadic.floats(g); sd=Dyadic.floats(s)
    d1=bd@bd.T+apply_exact(ed,gd)-gd
    d2=gd+apply_exact(ed,sd)-sd
    assert np.array_equal(d1.n,d1.n.T) and np.array_equal(d2.n,d2.n.T)
    a1,a2=d1.abs_sum(),d2.abs_sum()
    cnorm=(cd.T@cd).scalar(); center=(cd.T@sd@cd).scalar()
    correction=cnorm*(a1/(1-q)**2+a2/(1-q))
    upper=center+correction
    eta=sqrt_upper(upper)
    # Independently valid, coarser geometric certificate.
    bnorm=(bd.T@bd).scalar()
    generic=cnorm*bnorm/(1-q)**2
    if upper>generic:
        # Both are proved; retaining the smaller cannot invalidate the claim.
        used=generic; eta=sqrt_upper(used); bound_source='generic_geometric'
    else: used=upper; bound_source='exact_Stein_residual'
    return {
        'original_rank':n,'reduced_rank':k,'full_stability':full_stability,'projected_stability':small_stability,
        'numeric_hankel_hilbert_schmidt_squared':float(C@s@C),
        'exact_Stein_candidate_quadratic':rational(center),
        'exact_residual1_entrywise_l1':rational(a1),'exact_residual2_entrywise_l1':rational(a2),
        'exact_error_correction':rational(correction),'exact_Stein_upper_squared':rational(upper),
        'exact_generic_upper_squared':rational(generic),'selected_upper_squared':rational(used),
        'selected_bound_source':bound_source,'eta_rational':rational(eta),'eta_float':float_upper(eta),
        'eta_squared_ge_upper_exactly':eta*eta>=used,
        'stein_candidates':{'G0':gd.export(),'S0':sd.export()},
    },eta

def load_checkpoint(path):
    ck=torch.load(path,map_location='cpu',weights_only=True)
    a=ck['args']; vocab=ck['vocab']; outputs=len(vocab) if a['task']=='text' else 6
    model=benchmark.make_model(a['arm'],len(vocab),outputs,a['task'],a['seed'])
    model.load_state_dict(ck['state_dict'])
    if 'effective_matrices' not in ck: raise RuntimeError('Checkpoint lacks frozen effective matrices')
    model.frozen_matrices=ck['effective_matrices'].detach().clone()
    model.eval()
    return model,ck

def certify_checkpoint(path,k=4):
    model,ck=load_checkpoint(path)
    with torch.no_grad():
        r=model.matrices().detach().cpu().numpy().copy()
        b=model.write.detach().cpu().numpy().copy()
        c=model.head_read.detach().cpu().numpy().copy()
        v=model.output_vectors.weight.detach().cpu().numpy().copy()
    assert all(x.dtype==np.float32 for x in [r,b,c,v])
    q=F(99,100) if model.arm=='multi_cert' else F(999,1000)
    exports={'R_modes_heads':Dyadic.floats(r).export(),'b_heads':Dyadic.floats(b).export(),
             'c_heads':Dyadic.floats(c).export(),'v_output_by_head':Dyadic.floats(v).export()}
    # Retain decimal arrays too; every value exactly represents its FP32 source.
    values={'R':r.tolist(),'b':b.tolist(),'c':c.tolist(),'v':v.tolist()}
    results=[]; aggregate=F(0)
    for h in range(model.heads):
        result,eta=error_certificate(r[:,h],b[h],c[h],k,q)
        vinf=max(F.from_float(float(x)) for x in np.abs(v[:,h]))
        aggregate+=vinf*eta
        result.update(head=h,v_infinity_exact=rational(vinf))
        results.append(result)
    delta_bounds={t:sqrt_upper(aggregate*aggregate*t) for t in [32,64,96,128]}
    tv_bounds={t:min(F(1),bound/2) for t,bound in delta_bounds.items()}
    result={
        'checkpoint':path.name,'checkpoint_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
        'task':ck['args']['task'],'arm':model.arm,'seed':ck['args']['seed'],
        'q':rational(q),'projection_rank_per_head':k,'source_state_bytes':256,
        'projected_state_bytes':4*(model.main_dim+model.heads*k),
        'effective_coefficients_exact_dyadic':exports,'effective_coefficients_decimal':values,
        'heads':results,'all_stability_checks_passed':all(h['full_stability']['passed'] and h['projected_stability']['passed'] for h in results),
        'all_error_bounds_exactly_rounded_up':all(h['eta_squared_ge_upper_exactly'] for h in results),
        'uniform_logit_coefficient_A_exact':rational(aggregate),
        'uniform_logit_bound':'A * sqrt(t), because every scalar injection has absolute value at most one',
        'uniform_TV_bound':'min(1, A * sqrt(t) / 2)',
        'uniform_logit_bounds_exact_upward':{str(t):rational(bound) for t,bound in delta_bounds.items()},
        'uniform_TV_bounds_exact_upward':{str(t):rational(bound) for t,bound in tv_bounds.items()},
        'uniform_logit_bounds_float':{str(t):float_upper(bound) for t,bound in delta_bounds.items()},
        'uniform_TV_bounds_float':{str(t):float_upper(bound) for t,bound in tv_bounds.items()},
        'scope':'Exact real arithmetic on frozen effective FP32 coefficients; same routes and injections; zero initial states. No execution-roundoff guarantee. Feasible coordinate projection, not AAK optimum.'}
    return result

def self_test():
    # One-letter scalar f(n)=b*c*r^n has exact Hankel HS squared=(bc)^2/(1-r^2)^2.
    r=np.array([[[.5]]]); b=np.array([.75]); c=np.array([.25]); q=F(99,100)
    out,eta=error_certificate(r,b,c,0,q)
    expected=F(1,16)
    upper=F(int(out['selected_upper_squared']['numerator']),int(out['selected_upper_squared']['denominator']))
    assert upper>=expected and float(upper-expected)<1e-9
    # Identical projection must have zero true error; certified numerical residual may remain tiny.
    out0,eta0=error_certificate(r,b,c,1,q)
    assert float(eta0)<1e-5
    rng=np.random.default_rng(817)
    rr=rng.normal(size=(3,8,8)).astype(np.float32)*.04
    bb=rng.normal(size=8).astype(np.float32); cc=rng.normal(size=8).astype(np.float32)
    out1,eta1=error_certificate(rr,bb,cc,4,q)
    # Finite explicit Hankel is a lower bound on the infinite operator norm.
    words=[()]+[(a,) for a in range(3)]+[(a,b) for a in range(3) for b in range(3)]
    def value(word,k):
        x=bb[:k].astype(np.float64)
        for a in reversed(word): x=rr[a,:k,:k].astype(np.float64)@x
        return cc[:k]@x
    H=np.array([[value(u+v,8)-value(u+v,4) for v in words] for u in words])
    assert np.linalg.norm(H,2)<=float(eta1)+1e-12
    print(json.dumps({'self_test':'passed','scalar_exact_HS_squared':str(expected),'scalar_certificate':out['eta_float'],
                      'identity_projection_eta':float(eta0),'random_finite_Hankel_norm':float(np.linalg.norm(H,2)),
                      'random_certificate_eta':float(eta1)}))

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--self-test',action='store_true')
    parser.add_argument('--checkpoint',type=Path);parser.add_argument('--include-single',action='store_true')
    parser.add_argument('--rank',type=int,default=4); args=parser.parse_args()
    torch.set_num_threads(1)
    if args.self_test: self_test();return
    if args.checkpoint: paths=[args.checkpoint]
    else:
        paths=sorted((ROOT/'checkpoints').glob('*_multi_cert_seed*.pt'))
        if len(paths)!=10: raise RuntimeError(f'Expected 10 multi_cert checkpoints; found {len(paths)}')
        if args.include_single: paths+=sorted((ROOT/'checkpoints').glob('*_single_seed*.pt'))+sorted((ROOT/'checkpoints').glob('*_single_param_seed*.pt'))
    target=ROOT/'certificates';target.mkdir(exist_ok=True)
    summary=[]
    for path in paths:
        start=time.perf_counter();result=certify_checkpoint(path,args.rank)
        result['certificate_seconds']=time.perf_counter()-start
        output=target/(path.stem+'.json');output.write_text(json.dumps(result,indent=2)+'\n')
        row={k:result[k] for k in ['checkpoint','task','arm','seed','projected_state_bytes','all_stability_checks_passed','all_error_bounds_exactly_rounded_up','uniform_TV_bounds_float','certificate_seconds']}
        row['eta_per_head']=[h['eta_float'] for h in result['heads']]
        row['certificate_file']=output.name;summary.append(row)
        print(json.dumps(row),flush=True)
    (target/'summary.json').write_text(json.dumps({'runs':summary,'all_passed':all(r['all_stability_checks_passed'] and r['all_error_bounds_exactly_rounded_up'] for r in summary),'scope':'Frozen-coefficient exact-real sufficient stability and feasible Hankel-error upper certificate; no AAK optimum or execution-roundoff certificate.'},indent=2)+'\n')

if __name__=='__main__': main()
