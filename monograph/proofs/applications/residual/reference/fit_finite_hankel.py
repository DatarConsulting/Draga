#!/usr/bin/env python3
"""Supplementary candidate selection using a complete finite scalar word-Hankel.

The source is a constructed H=2 finite-support series from calibration controls.
Its exact-real objective differs from both Qwen TV and infinite predictive energy.
"""
import os
for _k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[_k]='1'
from pathlib import Path
import numpy as np
import json,csv,hashlib,itertools,time

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'hankel_finite'
OUT.mkdir(exist_ok=True)

def sha(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(n,v):(OUT/n).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')

plan={'source':'Finite-support scalar word series from calibration old memory',
      'calibration':['cal_river_c95','cal_stage_c95'],'control_offsets':[1,9],
      'action_alphabet_size':4,'read_label_count':2*16*128,
      'max_action_length':2,'action_weight':'19/40 = discount19/20 times sqrt(1/4)',
      'query':'average of four captured q/sqrt128 vectors per head',
      'source_rule':'f(read_(doc,head,j) a_word)=q_head^T A_word S_doc,head^T e_j for |word|<=2, zero otherwise',
      'candidate_rule':'Replace each S by B B^T S; same common rank96 B',
      'selection':'41 FP32 bases top96 of (1-theta) M_energy/trace(M_energy) + theta M_stein/trace(M_stein), theta=i/40; reuse existing endpoint bytes',
      'no_heldout_data_read':True,'global_optimality_claim':False,
      'arithmetic':'float64 numeric full finite norm evaluation; not interval-certified',
      'rank_scope':'rank96 refers to value codec; all constructed word-Hankels have ordinary rank at most22 because their only nonzero column blocks have21 and1 columns',
      'model_scope':'Full infinite word-Hankel of a declared finite-support surrogate, not the full Qwen sequence function'}
dump('PREDECLARED.json',plan)
t0=time.perf_counter()
d=np.load(ROOT/'bases/calibration_grams.npz')
S=d['states'];k=d['control_k'];beta=d['control_beta'];alpha=d['control_alpha'];q=d['control_q'].mean(axis=0)
Me=d['energy'];Ms=d['stein']
words=[()]+[(a,) for a in range(4)]+list(itertools.product(range(4),repeat=2))
I=np.eye(128)
Araw=(alpha[:,:,None,None]*(I[None,None]-beta[:,:,None,None]*np.einsum('mhi,mhj->mhij',k,k))).astype('<f4')
A=Araw.astype(np.float64)*(19/40)
Vu=[];labels=[]
for doc in range(2):
 for h in range(16):
  for u in words:
   row=q[h].copy()
   # A_u^T q = A_last^T ... A_first^T q; applying in word order.
   for a in u:row=A[a,h].T@row
   V=np.zeros((128,len(words)))
   for j,v in enumerate(words):
    if len(u)+len(v)>2:continue
    col=row.copy()
    for a in v:col=A[a,h].T@col
    V[:,j]=S[doc,h]@col
   Vu.append(V);labels.append((doc,h,list(u)))
V=np.asarray(Vu)
flat=V.transpose(1,0,2).reshape(128,-1)
source=np.einsum('biv,biw->vw',V,V,optimize=True)
source=(source+source.T)/2
source_norm=float(np.sqrt(np.linalg.eigvalsh(source)[-1]))
def objective(B):
 # Compute residual directly to avoid subtractive cancellation in source - projected.
 residual=flat-B@(B.T@flat)
 R=residual.reshape(128,len(V),21).transpose(1,0,2).reshape(-1,21)
 G=R.T@R
 return float(np.sqrt(max(0,np.linalg.eigvalsh((G+G.T)/2)[-1]))),G
def basis(theta):
 if theta in (0.,1.):
  label='energy' if theta==0. else 'stein'
  return np.fromfile(ROOT/'bases'/(label+'_r96.f32'),dtype='<f4').reshape(128,96).astype(float)
 M=(1-theta)*Me/np.trace(Me)+theta*Ms/np.trace(Ms)
 val,vec=np.linalg.eigh((M+M.T)/2);B=vec[:,-96:][:,::-1]
 for j in range(96):
  if B[np.argmax(abs(B[:,j])),j]<0:B[:,j]*=-1
 return B.astype('<f4').astype(float)
rows=[];best=None
for i in range(41):
 theta=i/40;B=basis(theta);norm,G=objective(B)
 R=I-B@B.T
 row={'grid_index':i,'theta':theta,'hankel_error':norm,'relative_hankel_error':norm/source_norm,
      'energy_residual':float(np.trace(R.T@Me@R)),
      'stein_residual':float(np.trace(R.T@Ms@R))}
 rows.append(row)
 if best is None or norm<best[0]:best=(norm,theta,B.copy(),G.copy())
norm,theta,B,G=best
B32=B.astype('<f4')
endpoint_reused=None
if theta in (0.,1.):
 endpoint_reused='energy' if theta==0. else 'stein'
 # Reuse the existing endpoint artifact, avoiding arbitrary eigensolver roundoff
 # creating a nominally distinct deployment arm for the same selected subspace.
 B32=np.fromfile(ROOT/'bases'/(endpoint_reused+'_r96.f32'),dtype='<f4').reshape(128,96)
p=OUT/'hankel_r96.f32';B32.tofile(p)
deployed_norm,G32=objective(B32.astype(float))
# Independent full nonzero block SVD on selected candidate.
residual=(flat-B@(B.T@flat)).reshape(128,len(V),21).transpose(1,0,2).reshape(-1,21)
singular=np.linalg.svd(residual,compute_uv=False)
assert abs(float(singular[0])-norm)<=1e-10*max(1,norm)
# The disjoint empty-prefix/read-column block has l2 norm equal to first column.
other_block=float(np.linalg.norm(residual[:,0]));assert other_block<=norm*(1+1e-12)
endpoints={}
for name in ('energy','stein'):
 p0=ROOT/'bases'/(name+'_r96.f32');B0=np.fromfile(p0,dtype='<f4').reshape(128,96).astype(float)
 z,_=objective(B0);endpoints[name]={'basis_sha256':sha(p0),'deployed_hankel_error':z,'relative':z/source_norm}
np.savez(OUT/'source_model.npz',states=S,unweighted_actions_f32=Araw,query=q,
         energy_gram=Me,stein_gram=Ms,source_block_gram=source,selected_error_gram=G,
         selected_deployed_error_gram=G32,basis64=B,basis32=B32)
with (OUT/'grid.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
summary={'plan':plan,'source_input_sha256':sha(ROOT/'bases/calibration_grams.npz'),
         'source_action_parameter_precision':'Unweighted A matrices rounded once to FP32; exact-real source multiplies these dyadic matrices by rational19/40',
         'source_norm_numeric':source_norm,'nonzero_main_block_shape':list(residual.shape),
         'selected_theta':theta,'selected_hankel_error':norm,'selected_relative_error':norm/source_norm,
         'deployed_basis_hankel_error':deployed_norm,'other_disjoint_block_norm':other_block,
         'svd_vs_gram_discrepancy':abs(float(singular[0])-norm),
         'basis_path':str(p.resolve()),'basis_sha256':sha(p),'basis_bytes':p.stat().st_size,
         'canonical_existing_endpoint_reused':endpoint_reused,
         'endpoints':endpoints,'improvement_vs_energy':1-deployed_norm/endpoints['energy']['deployed_hankel_error'],
         'improvement_vs_stein':1-deployed_norm/endpoints['stein']['deployed_hankel_error'],
         'seconds':time.perf_counter()-t0,'grid_candidates':len(rows),
         'claim':'Best of41 evaluated interpolation bases only. No global optimum or Qwen probability certificate.'}
dump('summary.json',summary)
print(json.dumps(summary,indent=2))
