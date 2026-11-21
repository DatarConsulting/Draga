#!/usr/bin/env python3
"""Calibration-only common value bases; weighted Stein proxy, not a Qwen certificate."""
import os
for _key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[_key]='1'
import argparse, hashlib, json, time
from pathlib import Path
import numpy as np

DOCS=('cal_river','cal_stage')
OFFSETS=(1,9)
RHO=0.95
HEADS=16
DIM=128
RANK=96

def digest(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

def savejson(p,obj):p.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')

def principal_basis(M,r):
    vals,vecs=np.linalg.eigh((M+M.T)/2)
    ix=np.argsort(vals)[::-1]; vals=vals[ix]; B=vecs[:,ix[:r]]
    for j in range(r):
        if B[np.argmax(abs(B[:,j])),j]<0:B[:,j]*=-1
    return vals,B

def apply_phi(Q,k,beta,alpha):
    # mode,head,key input; Q=head,key,key; exact-real algebra of A^T Q A
    qk=np.einsum('hij,mhj->mhi',Q,k,optimize=True)
    kk=np.einsum('mhi,mhj->mhij',k,k,optimize=True)
    dot=np.einsum('mhi,mhi->mh',k,qk,optimize=True)
    cross=np.einsum('mhi,mhj->mhij',k,qk,optimize=True)
    vals=Q[None,:,:,:]-beta[:,:,None,None]*(cross+cross.swapaxes(-1,-2))
    vals+=((beta**2)*dot)[:,:,None,None]*kk
    return RHO**2*np.mean(alpha[:,:,None,None]**2*vals,axis=0)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--calibration-root',type=Path,default=Path('qwen_test/calibration'))
    ap.add_argument('--output',type=Path,default=Path('qwen_test/bases'))
    args=ap.parse_args(); args.output.mkdir(parents=True,exist_ok=True)
    plan={'version':1,'purpose':'Old-memory predictive basis using a fixed weighted command tree',
          'calibration_docs':list(DOCS),'calibration_cut_token1':95,'calibration_horizon':32,
          'calibration_control_offsets':list(OFFSETS),'action_modes':4,
          'mode_probabilities':[.25]*4,'discount':RHO,'rank':RANK,
          'basis_shared_by_heads':HEADS,'layout':'C contiguous float32 [value=128,rank=96]',
          'selection':'one fixed rank, no heldout logits or states used',
          'heldout_distribution_shift':'calibration cut95 versus evaluation cut428',
          'scope':'This weighted control-tree Gram optimizes squared future raw-readout loss. It is not the ordinary scalar word-Hankel norm or a certificate for Qwen logits.',
          'query_readout_scale':'captured q_conv_predelta is divided by sqrt(128) inside native GDN readout',
          'arithmetic':'float64 solve, float32 basis deployment; diagnostic residual, not interval certified'}
    savejson(args.output/'fit_plan.json',plan)
    receipts=[]; states=[]; keys=[]; queries=[]; betas=[]; alphas=[]; recurrence_checks=[]
    for doc in DOCS:
        folder=args.calibration_root/(doc+'_c95')
        meta=json.loads((folder/'metadata.json').read_text())
        if meta['cut_after_token_number1']!=95 or meta['horizon']!=32:raise ValueError('calibration scope')
        data={}
        for name,shape in [('native_state',(16,128,128)),('q',(33,16,128)),('k',(33,16,128)),
                           ('v',(33,16,128)),('g',(33,16)),('beta',(33,16)),('readout',(33,16,128))]:
            p=folder/(name+'.f32'); v=np.fromfile(p,dtype='<f4').astype(np.float64).reshape(shape)
            if not np.isfinite(v).all():raise ValueError('nonfinite '+str(p))
            data[name]=v;receipts.append({'path':str(p.resolve()),'sha256':digest(p),'bytes':p.stat().st_size})
        states.append(data['native_state'])
        for off in OFFSETS:
            keys.append(data['k'][off]);queries.append(data['q'][off]/np.sqrt(DIM));betas.append(data['beta'][off]);alphas.append(np.exp(data['g'][off]))
        # Validate key-side orientation and captured log-decay on first next read.
        S=data['native_state']; alpha=np.exp(data['g'][1]); beta=data['beta'][1]; k=data['k'][1]
        v=data['v'][1]; q=data['q'][1]/np.sqrt(DIM)
        decayed=alpha[:,None,None]*S
        prediction=np.einsum('hvk,hk->hv',decayed,k)
        updated=decayed+beta[:,None,None]*(v-prediction)[:,:,None]*k[:,None,:]
        raw=np.einsum('hvk,hk->hv',updated,q)
        observed=data['readout'][1]
        relative_error=float(np.linalg.norm(raw-observed)/np.linalg.norm(observed))
        if relative_error>1e-5:raise RuntimeError('captured recurrence mismatch '+doc)
        recurrence_checks.append({'doc':doc,'relative_readout_l2_error':relative_error,
                                  'max_readout_abs_error':float(np.max(abs(raw-observed)))})
    k=np.asarray(keys); q=np.asarray(queries); beta=np.asarray(betas); alpha=np.asarray(alphas)
    keynorm2=np.sum(k*k,axis=-1)
    normA=alpha*np.maximum(1,abs(1-beta*keynorm2))
    contractions=RHO**2*np.mean(normA**2,axis=0)
    unweighted_bounds=np.sum(normA**2,axis=0)
    if np.any(contractions>=1):raise RuntimeError('weighted contraction failure')
    C=np.mean(np.einsum('mhi,mhj->mhij',q,q,optimize=True),axis=0)
    Q=np.zeros_like(C);started=time.perf_counter()
    for iteration in range(1,2001):
        nextQ=C+apply_phi(Q,k,beta,alpha)
        change=float(np.linalg.norm(nextQ-Q));Q=nextQ
        if change<=1e-13*max(float(np.linalg.norm(Q)),1e-30):break
    else:raise RuntimeError('Stein iteration did not converge')
    residual=C+apply_phi(Q,k,beta,alpha)-Q
    residual_norm=np.linalg.norm(residual,ord=2,axis=(-2,-1))
    gram_error_estimate=residual_norm/(1-contractions)
    M_energy=sum(np.einsum('hvk,hwk->vw',S,S,optimize=True) for S in states)/(len(states)*HEADS)
    M_stein=sum(np.einsum('hvi,hij,hwj->vw',S,Q,S,optimize=True) for S in states)/(len(states)*HEADS)
    np.savez(args.output/'calibration_grams.npz',energy=M_energy,stein=M_stein,Q=Q,C=C,
             control_k=k,control_q=q,control_beta=beta,control_alpha=alpha,states=np.asarray(states))
    bases={};basis_stats={}
    for label,M in [('energy',M_energy),('stein',M_stein)]:
        vals,B=principal_basis(M,RANK);B32=B.astype('<f4');path=args.output/(label+'_r96.f32');B32.tofile(path)
        bases[label]=B
        P=B@B.T;P32=B32.astype(np.float64)@B32.astype(np.float64).T
        basis_stats[label]={'basis_file':str(path.resolve()),'sha256':digest(path),
            'basis_bytes':path.stat().st_size,'float64_orthogonality_error':float(np.linalg.norm(B.T@B-np.eye(RANK),2)),
            'float32_orthogonality_error':float(np.linalg.norm(B32.astype(float).T@B32.astype(float)-np.eye(RANK),2)),
            'eigenvalues':vals.tolist(),'own_objective_discarded_fraction':float(np.sum(vals[RANK:])/np.sum(vals)),
            'energy_residual':float(np.trace((np.eye(DIM)-P)@M_energy)),
            'stein_residual':float(np.trace((np.eye(DIM)-P)@M_stein)),
            'float32_energy_residual':float(np.trace((np.eye(DIM)-P32).T@M_energy@(np.eye(DIM)-P32))),
            'float32_stein_residual':float(np.trace((np.eye(DIM)-P32).T@M_stein@(np.eye(DIM)-P32)))}
    summary={'plan':plan,'input_receipts':receipts,'recurrence_checks':recurrence_checks,
        'keys_norm2_minmax':[float(keynorm2.min()),float(keynorm2.max())],
        'alpha_minmax':[float(alpha.min()),float(alpha.max())],'beta_minmax':[float(beta.min()),float(beta.max())],
        'weighted_contraction_bounds_by_head':contractions.tolist(),
        'unweighted_sum_norm_squared_bounds_by_head':unweighted_bounds.tolist(),
        'unweighted_stability_note':'A bound above one does not prove instability; no unweighted counting-Hankel stability claim is made.',
        'iterations':iteration,'fit_seconds':time.perf_counter()-started,
        'stein_residual_spectral_by_head':residual_norm.tolist(),
        'numerical_gram_error_estimate_by_head':gram_error_estimate.tolist(),
        'basis_statistics':basis_stats,
        'payload_byte_accounting':{'dense_fp16':HEADS*DIM*DIM*2,'reduced_fp16':HEADS*DIM*RANK*2,
             'common_fp32_basis':DIM*RANK*4,'reduced_plus_basis':HEADS*DIM*RANK*2+DIM*RANK*4,
             'ratio_vs_dense_fp16':(HEADS*DIM*RANK*2+DIM*RANK*4)/(HEADS*DIM*DIM*2),
             'not_counted':'native FP32 caches, temporaries, encoder work; no actual allocated-memory reduction claimed'}}
    savejson(args.output/'fit_summary.json',summary)
    print(json.dumps({'status':'complete','iterations':iteration,'recurrence_checks':recurrence_checks,
                     'contraction_max':float(contractions.max()),'basis_statistics':{x:{k:v for k,v in row.items() if k!='eigenvalues'} for x,row in basis_stats.items()},
                     'payload_bytes':summary['payload_byte_accounting']},indent=2))

if __name__=='__main__':main()
