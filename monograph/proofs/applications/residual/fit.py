#!/usr/bin/env python3
"""Calibration-only reconstruction of the archived finite Hankel selector."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[key]='1'
import json,time,hashlib,itertools
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parent
def digest(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    t=time.perf_counter(); out=ROOT/'bases'; out.mkdir(exist_ok=True)
    d=np.load(ROOT/'inputs/calibration_grams.npz')
    S=d['states']; k=d['control_k']; beta=d['control_beta']; alpha=d['control_alpha']; q=d['control_q'].mean(0)
    Me=d['energy']; Ms=d['stein']; I=np.eye(128)
    A=(alpha[:,:,None,None]*(I[None,None]-beta[:,:,None,None]*np.einsum('mhi,mhj->mhij',k,k))).astype('<f4').astype(float)*(19/40)
    words=[()]+[(a,) for a in range(4)]+list(itertools.product(range(4),repeat=2))
    blocks=[]
    for doc in range(2):
        for h in range(16):
            for u in words:
                row=q[h].copy()
                for a in u: row=A[a,h].T@row
                V=np.zeros((128,21))
                for j,v in enumerate(words):
                    if len(u)+len(v)>2: continue
                    col=row.copy()
                    for a in v: col=A[a,h].T@col
                    V[:,j]=S[doc,h]@col
                blocks.append(V)
    V=np.asarray(blocks); flat=V.transpose(1,0,2).reshape(128,-1)
    def obj(B):
        R=(flat-B@(B.T@flat)).reshape(128,len(V),21).transpose(1,0,2).reshape(-1,21)
        return float(np.sqrt(max(0,np.linalg.eigvalsh(R.T@R)[-1])))
    ranks=json.loads((ROOT/'contract.json').read_text())['rank_sweep']
    best={r:None for r in ranks}; rows=[]; energy={}
    for i in range(41):
        theta=i/40
        M=(1-theta)*Me/np.trace(Me)+theta*Ms/np.trace(Ms)
        _,B=np.linalg.eigh((M+M.T)/2); B=B[:,::-1].copy()
        for j in range(128):
            if B[np.argmax(abs(B[:,j])),j]<0: B[:,j]*=-1
        for r in ranks:
            Br=B[:,:r].astype('<f4')
            if r==96 and i in (0,40):
                fn='energy_r96.f32' if i==0 else 'hankel_r96.f32'
                Br=np.fromfile(ROOT/'inputs'/fn,dtype='<f4').reshape(128,96)
            value=obj(Br.astype(float)); rows.append({'rank':r,'theta':theta,'hankel_error':value})
            if i==0: energy[r]=Br.copy()
            if best[r] is None or value<best[r][0]: best[r]=(value,theta,Br.copy())
        print('Hankel selection',i+1,'/41',flush=True)
    summaries=[]
    for r in ranks:
        value,theta,B=best[r]
        if r==96:
            assert theta==1.0
            assert np.array_equal(B,np.fromfile(ROOT/'inputs/hankel_r96.f32',dtype='<f4').reshape(128,96))
        for label,base in [('energy',energy[r]),('hankel',B)]:
            p=out/f'{label}_{r}.f32'; base.T.copy().tofile(p)
            summaries.append({'method':label,'rank':r,'theta':0.0 if label=='energy' else theta,
                              'hankel_error':obj(base.astype(float)),'file':p.name,'bytes':p.stat().st_size,'sha256':digest(p)})
    (out/'selection.json').write_text(json.dumps({'basis_layout':'rank,value (transpose of archived files)','rows':summaries,'grid':rows,'seconds':time.perf_counter()-t},indent=2)+'\n')
if __name__=='__main__': main()
