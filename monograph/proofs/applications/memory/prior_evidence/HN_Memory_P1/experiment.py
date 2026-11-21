#!/usr/bin/env python3
"""Approximate invariant flag, rational comparison certificate, held-out tests.

The flag is inherited from P0; it is not a new Harder--Narasimhan filtration
of the perturbed representation. Four state coordinates are retained.
"""
from pathlib import Path
import csv
import hashlib
import itertools
import json
import platform
import time
import numpy as np
import sympy as sp
from codec import candidates,pack,unpack

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results'
EPSILONS=('0','0.000001','0.0001','0.001','0.01')
BUDGETS=(16,24,32)
HORIZONS=(1,8,32)
GAINS=list(itertools.permutations((8.,2.,.5,.125)))
DIAGS=np.array([[.98,.90,.65,.25],[.95,.75,.30,-.20]])
B=np.array([.1,.25,.5,.8])


def save_json(path,obj):
    path.write_text(json.dumps(obj,indent=2,ensure_ascii=False)+'\n')


def rational_setup(ep):
    e=sp.Rational(ep)
    a0=sp.diag(*map(sp.Rational,('98/100','9/10','65/100','1/4')))
    a1=sp.diag(*map(sp.Rational,('95/100','3/4','3/10','-1/5')))+e*(sp.ones(4)-sp.eye(4))
    K=a0+e*(sp.ones(4)-sp.eye(4))
    b=sp.Matrix([sp.Rational(1,10),sp.Rational(1,4),sp.Rational(1,2),sp.Rational(4,5)])
    inv=(sp.eye(4)-K).inv()
    ranges=inv*b
    if not all(x>0 for x in ranges):
        return {'epsilon':ep,'status':'NOT_CERTIFIED_BY_THIS_COMPARISON',
                'reason':'The rational solution R=(I-K)^(-1)b is not strictly positive; no stability conclusion.'}
    assert K*ranges+b==ranges
    assert all((K*ranges)[i]<ranges[i] for i in range(4))
    assert all(v>=0 for v in inv)
    # Defect: lower-left (4-k)-by-k constant matrix epsilon.
    defects=[]
    for k in (1,2,3):
        P=sp.diag(*([1]*k+[0]*(4-k)))
        L=(sp.eye(4)-P)*a1*P
        assert sum(v*v for v in L)==e*e*k*(4-k)
        assert L.rank()==(0 if e==0 else 1)
        defects.append(str(e*sp.sqrt(k*(4-k))))
    T1=sp.Matrix(4,4,lambda i,j:a1[i,j] if i<=j else 0)
    E=a1-T1
    for k in (1,2,3):
        P=sp.diag(*([1]*k+[0]*(4-k)))
        assert (sp.eye(4)-P)*T1*P==sp.zeros(4)
    structural_bound=inv*E*ranges
    record={'epsilon':ep,'status':'EXACT_RATIONAL_WITNESS',
            'R':[str(x) for x in ranges],
            'weighted_contraction':str(max((K*ranges)[i]/ranges[i] for i in range(4))),
            'flag_defects_spectral_norm':defects,'max_flag_defect':str(2*e),
            'triangular_defect_infinity_norm':str(3*e),
            'structural_state_bound':[str(x) for x in structural_bound],
            'KR_plus_b_equals_R':True,'inverse_nonnegative':True,
            'triangular_transition_preserves_flag':True,
            'original_nontrivial_common_invariant_memory_subspaces':bool(e==0)}
    return {'record':record,'A':np.array([a0.tolist(),a1.tolist()],dtype=float),
            'K':np.array(K,dtype=float),'inverse':np.array(inv,dtype=float),
            'R':np.array(ranges,dtype=float).reshape(4),
            'D':np.array(E,dtype=float),'rational':(a0,a1,K,b,ranges,inv,T1)}


def make_data(seed,n,length):
    rng=np.random.default_rng(seed)
    tokens=rng.integers(0,2,size=(length,n),dtype=np.int8)
    numerators=rng.integers(-1024,1025,size=(length,n),dtype=np.int16)
    return tokens,numerators.astype(float)/1024,numerators


def apply_transition(q,a,epsilon,triangular=False,absolute=False):
    diag=np.abs(DIAGS[a]) if absolute else DIAGS[a]
    y=q*diag
    if triangular:
        other=np.flip(np.cumsum(np.flip(q,axis=-1),axis=-1),axis=-1)-q
    else:
        other=np.sum(q,axis=-1,keepdims=True)-q
    return y+epsilon*other*(a[...,None]==1)


def evaluate(data,setup,bits,budget,triangular,burn=64):
    tokens,u,_=data
    length,n=tokens.shape
    k=len(bits); eps=float(setup['record']['epsilon']); ranges=setup['R']
    mid=2**(bits-1)-1; step=ranges/mid; delta=step/2
    D=setup['D'] if triangular else np.zeros((4,4))
    source=delta+(D@ranges)[None,:]
    global_v=source@setup['inverse'].T
    global_v=np.minimum(global_v,2*ranges)
    payload=pack(np.broadcast_to(mid[:,None,:],(k,n,4)),bits,budget)
    x=np.zeros((n,4)); bound=np.zeros((k,n,4))
    mse=np.zeros((k,n,4)); bound_mse=np.zeros((k,n,4)); energy=np.zeros((n,4))
    maxerr=np.zeros((k,4)); boundviol=0; globalviol=0; clip=0; count=0
    worst_signed_excess=-1e100
    for t in range(length):
        a=tokens[t]
        q=(unpack(payload,bits).astype(float)-mid[:,None,:])*step[:,None,:]
        pre=apply_transition(q,a,eps,triangular)+B*u[t,None,:,None]
        clip+=int(np.count_nonzero(np.abs(pre)>ranges+1e-12))
        codes=np.rint(pre/step[:,None,:]).astype(np.int64)+mid[:,None,:]
        codes=np.clip(codes,0,2*mid[:,None,:])
        payload=pack(codes,bits,budget)
        qnew=(unpack(payload,bits).astype(float)-mid[:,None,:])*step[:,None,:]
        eta=qnew-pre
        leak=(np.cumsum(q,axis=-1)-q)*eps*(a[None,:,None]==1) if triangular else 0
        # This a posteriori bound needs the codec state and its local defect,
        # never the native/reference state x.
        bound=apply_transition(bound,a,eps,absolute=True)+np.abs(leak-eta)
        x=apply_transition(x,a,eps)+B*u[t,:,None]
        error=np.abs(qnew-x[None,:,:])
        excess=error-bound
        worst_signed_excess=max(worst_signed_excess,float(excess.max()))
        boundviol+=int(np.count_nonzero(excess>1e-10))
        globalviol+=int(np.count_nonzero(error>global_v[:,None,:]+1e-10))
        maxerr=np.maximum(maxerr,error.max(axis=1))
        if t>=burn:
            mse+=error**2; bound_mse+=bound**2; energy+=x**2; count+=1
    assert payload.nbytes==k*n*budget//8
    assert clip==0 and boundviol==0 and globalviol==0
    return {'mse':mse/count,'bound_mse':bound_mse/count,'energy':energy/count,
            'global_v':global_v,'max_error':maxerr,
            'checks':{'scored_steps_per_sequence':count,'clipping_count':clip,
                      'aposteriori_bound_violations_at_1e-10':boundviol,
                      'global_bound_violations_at_1e-10':globalviol,
                      'max_signed_aposteriori_excess_float64':worst_signed_excess,
                      'scalar_coordinate_checks':k*n*length*4}}


def structural_probe(data,setup,burn=64):
    tokens,u,_=data; n=tokens.shape[1]
    x=np.zeros((n,4)); q=x.copy(); v=x.copy(); errsq=x.copy(); maxerr=np.zeros(4)
    eps=float(setup['record']['epsilon']); count=0
    glob=setup['inverse']@setup['D']@setup['R']
    for t,a in enumerate(tokens):
        defect=(np.cumsum(q,axis=-1)-q)*eps*(a[:,None]==1)
        v=apply_transition(v,a,eps,absolute=True)+np.abs(defect)
        q=apply_transition(q,a,eps,triangular=True)+B*u[t,:,None]
        x=apply_transition(x,a,eps)+B*u[t,:,None]
        err=np.abs(x-q)
        assert np.all(err<=v+1e-11) and np.all(err<=glob+1e-11)
        maxerr=np.maximum(maxerr,err.max(axis=0))
        if t>=burn:errsq+=err**2;count+=1
    gains=np.array(GAINS[0]); reverse=np.array(GAINS[-1])
    return {'epsilon':setup['record']['epsilon'],
            'max_state_error_observed':float(maxerr.max()),
            'uniform_state_error_bound':float(glob.max()),
            'max_output_error_observed_gains_initial':float(np.max(gains*maxerr)),
            'uniform_output_error_bound_gains_initial':float(np.max(gains*glob)),
            'test_mse_initial_gains':float(np.mean(errsq/count@gains**2)/4),
            'test_mse_reversed_gains':float(np.mean(errsq/count@reverse**2)/4)}


def main():
    OUT.mkdir(exist_ok=True); start=time.time()
    cal=make_data(20261111,64,256);test=make_data(20261112,128,512)
    np.savez_compressed(OUT/'datasets.npz',calibration_tokens=cal[0],calibration_input_numerators=cal[2],
                        test_tokens=test[0],test_input_numerators=test[2],input_denominator=1024)
    rng=np.random.default_rng(20261113); bootstrap=rng.integers(0,128,size=(500,128))
    rows=[];records=[];valid=[];struct=[];evidence={}
    for ep in EPSILONS:
        setup=rational_setup(ep);records.append(setup['record'])
        struct.append(structural_probe(test,setup))
        for budget in BUDGETS:
            bits=candidates(budget)
            uniform=int(np.flatnonzero(np.all(bits==budget//4,axis=1))[0])
            for triangular in (False,True):
                mode='triangular' if triangular else 'full'
                ca=evaluate(cal,setup,bits,budget,triangular)
                te=evaluate(test,setup,bits,budget,triangular)
                prefix=f'{ep}_{budget}_{mode}'
                evidence[prefix+'_bits']=bits
                for name in ('mse','bound_mse','energy','global_v','max_error'):
                    evidence[prefix+'_test_'+name]=te[name]
                evidence[prefix+'_calibration_mse']=ca['mse']
                valid.append({'epsilon':ep,'budget':budget,'transition':mode,
                              'calibration':ca['checks'],'test':te['checks']})
                for case,g in enumerate(GAINS):
                    g2=np.square(g)
                    cal_mse=np.mean(ca['mse']@g2/4,axis=1)
                    calidx=int(np.argmin(cal_mse))
                    boundidx=int(np.argmin(te['global_v']**2@g2/4))
                    base=te['mse'][uniform]@g2/4
                    for pol,idx in [('uniform',uniform),('calibrated',calidx),('global_certificate',boundidx)]:
                        seq=te['mse'][idx]@g2/4
                        reductions=100*(1-seq[bootstrap].mean(axis=1)/base[bootstrap].mean(axis=1))
                        rows.append({'epsilon':ep,'budget_bits':budget,'persistent_state_bytes':budget//8,
                            'transition':mode,'case':case,'gains':list(g),'policy':pol,'bits':bits[idx].tolist(),
                            'test_mse':float(seq.mean()),
                            'test_nmse_pct':float(100*seq.mean()/np.mean(te['energy']@g2/4)),
                            'reduction_vs_uniform_same_transition_pct':float(100*(1-seq.mean()/base.mean())),
                            'paired_bootstrap_95ci_reduction':[float(x) for x in np.quantile(reductions,[.025,.975])],
                            'rms_aposteriori_bound':float(np.sqrt(np.mean(te['bound_mse'][idx]@g2/4))),
                            'uniform_rms_output_bound':float(np.sqrt(te['global_v'][idx]**2@g2/4)),
                            'observed_rms_error':float(np.sqrt(seq.mean()))})
            print('epsilon',ep,'budget',budget,'done',round(time.time()-start,1),'s',flush=True)
    summary=[]
    for ep,budget,mode,pol in itertools.product(EPSILONS,BUDGETS,('full','triangular'),('calibrated','global_certificate')):
        rr=[r for r in rows if (r['epsilon'],r['budget_bits'],r['transition'],r['policy'])==(ep,budget,mode,pol)]
        vals=np.array([r['reduction_vs_uniform_same_transition_pct'] for r in rr])
        summary.append({'epsilon':ep,'budget_bits':budget,'transition':mode,'policy':pol,
                        'better_cases':int(np.sum(vals>1e-9)),'equal_cases':int(np.sum(np.abs(vals)<=1e-9)),
                        'worse_cases':int(np.sum(vals< -1e-9)),
                        'median_reduction_pct':float(np.median(vals)),
                        'worst_reduction_pct':float(vals.min())})
    records.append(rational_setup('0.05'))
    for name,obj in [('certificates',records),('results',rows),('summary',summary),('structural_only',struct),('validation',valid)]:
        save_json(OUT/(name+'.json'),obj)
    with (OUT/'results.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    np.savez_compressed(OUT/'per_sequence_evidence.npz',**evidence)
    save_json(OUT/'protocol.json',{'version':'P1-1.0','date':'2026-11-01','epsilons':EPSILONS,
        'budgets_bits':BUDGETS,'gains':GAINS,'selection':'calibration-only post-storage output MSE; exact ties lexicographic',
        'metric_change_from_P0':'P1 scores the deployed recurrent output after every quantization; no quantization-free future probe',
        'calibration':{'sequences':64,'length':256,'burn':64,'seed':20261111},
        'test':{'sequences':128,'length':512,'burn':64,'seed':20261112},
        'input':'iid discrete uniform integers -1024..1024 divided by 1024; external controls shared by both systems',
        'bootstrap':'500 paired resamples of whole independent test sequences, seed 20261113, no multiplicity correction',
        'gain_cases_share_trajectories':True,'test_not_used_for_selection':True,
        'certificate_selector':'minimize global analytic bound, independent of calibration and test data',
        'persistent_memory':'2, 3, or 4 packed bytes for four state coordinates; same for full/triangular codecs',
        'certificate_memory':'global certificate static per codec: zero additional per-stream state; online 4-component aposteriori audit costs 32 float64 bytes per stream if deployed',
        'temporary_buffers':'float64 decoding and calculation buffers not included in persistent state; no speed or peak-RAM claims',
        'weights':'triangular variant removes six off-diagonal transition entries; no reduction in state dimension',
        'arithmetic_scope':'proof is real arithmetic; matrix witnesses and exact rational replay are exact; bulk tests float64 with explicit 1e-10 verification tolerance',
        'environment':{'python':platform.python_version(),'numpy':np.__version__,'sympy':sp.__version__},
        'elapsed_seconds':time.time()-start})
    print('COMPLETE',len(rows),'rows',round(time.time()-start,2),'seconds',flush=True)


if __name__=='__main__':main()
