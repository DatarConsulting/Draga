"""Prespecified rank8->4 replay. Numeric observations, not a roundoff proof."""
import json, math, os
from pathlib import Path
os.environ['OMP_NUM_THREADS']='1'
os.environ['MKL_NUM_THREADS']='1'
import numpy as np
import torch
from torch.nn import functional as F
import benchmark as b
from certify import load_checkpoint
ROOT=Path(__file__).resolve().parents[1]

@torch.no_grad()
def replay(model,batches,cert):
    eta=torch.tensor([h['eta_float'] for h in cert['heads']],dtype=torch.float64)
    vinf=model.output_vectors.weight.detach().double().abs().amax(0)
    sums=dict(tv=0.,argmax_flip=0.,numeric_bound_violations=0,nonvacuous_TV_bounds=0,
              sum_bound=0.,max_tv=0.,max_logit_error=0.,max_bound=0.,max_error_bound_ratio=0.)
    n=0
    for x in batches:
        z,main,res,h,u,routes=model(x,return_parts=True)
        zz,mm,rr,hh,uu,routes2=model(x,project_rank=4,return_parts=True)
        assert torch.equal(main,mm) and torch.equal(u,uu) and torch.equal(routes,routes2)
        energy=torch.sqrt(torch.cumsum(u.double().square(),dim=1))
        delta=(energy*eta[None,None,:]*vinf[None,None,:]).sum(-1)
        tv_bound=torch.clamp(delta/2,max=1.)
        error=(z.double()-zz.double()).abs().amax(-1)
        tv=(torch.softmax(z.double(),-1)-torch.softmax(zz.double(),-1)).abs().sum(-1)/2
        flip=(z.argmax(-1)!=zz.argmax(-1)).double()
        n+=tv.numel()
        sums['tv']+=tv.sum().item();sums['argmax_flip']+=flip.sum().item()
        sums['sum_bound']+=tv_bound.sum().item()
        sums['nonvacuous_TV_bounds']+=int((tv_bound<1.).sum())
        # Merely a numeric consistency check; tolerance does not establish rounding safety.
        sums['numeric_bound_violations']+=int((error>delta+1e-5).sum())
        sums['max_tv']=max(sums['max_tv'],tv.max().item())
        sums['max_logit_error']=max(sums['max_logit_error'],error.max().item())
        sums['max_bound']=max(sums['max_bound'],delta.max().item())
        sums['max_error_bound_ratio']=max(sums['max_error_bound_ratio'],(error/(delta+1e-30)).max().item())
    return dict(positions=n,tv_mean=sums['tv']/n,tv_max=sums['max_tv'],
                argmax_flip_fraction=sums['argmax_flip']/n,max_logit_error=sums['max_logit_error'],
                computed_TV_upper_mean=sums['sum_bound']/n,
                fraction_computed_TV_upper_below_one=sums['nonvacuous_TV_bounds']/n,
                max_computed_logit_upper=sums['max_bound'],
                max_error_to_computed_bound_ratio=sums['max_error_bound_ratio'],
                numeric_bound_violations_with_1e_minus5_tolerance=sums['numeric_bound_violations'])

def main():
    torch.set_num_threads(1);torch.use_deterministic_algorithms(True)
    data,_=b.load_text();out=[]
    paths=sorted((ROOT/'checkpoints').glob('*_multi_cert_seed*.pt'))
    if len(paths)!=10: raise RuntimeError('Expected all10 multi_cert checkpoints')
    for path in paths:
        model,ck=load_checkpoint(path);task=ck['args']['task'];seed=ck['args']['seed']
        cert=json.loads((ROOT/'certificates'/f'{path.stem}.json').read_text())
        original=json.loads((ROOT/'results'/f'{path.stem}.json').read_text())
        record=dict(task=task,seed=seed,arm='multi_cert',original_state_bytes=256,projected_state_bytes=192)
        if task=='text':
            full=b.evaluate_text(model,data['test'])
            assert abs(full['bpc']-original['test']['bpc'])<1e-8
            projected=b.evaluate_text(model,data['test'],project_rank=4)
            starts=np.arange(0,len(data['test'])-64,64)
            batches=[b.text_windows(data['test'],starts[o:o+64],64)[0] for o in range(0,len(starts),64)]
            record.update(original=full,projected=projected,
                          projected_validation=b.evaluate_text(model,data['validation'],project_rank=4),
                          numeric_replay=replay(model,batches,cert))
        else:
            full={str(L):b.evaluate_order(model,L,18000+L) for L in [32,64,128]}
            assert full['32']['accuracy']==original['test']['32']['accuracy']
            projected={str(L):b.evaluate_order(model,L,18000+L,project_rank=4) for L in [32,64,128]}
            for L,sd in [(32,19001),(96,19002)]:
                full['pairs'+str(L)]=b.evaluate_order(model,L,sd,n=512,paired=True)
                projected['pairs'+str(L)]=b.evaluate_order(model,L,sd,n=512,paired=True,project_rank=4)
            x=np.random.default_rng(18000+128).integers(0,3,size=(1024,128),dtype=np.int64)
            batches=[torch.from_numpy(x[o:o+128]) for o in range(0,len(x),128)]
            record.update(original=full,projected=projected,numeric_replay=replay(model,batches,cert))
        record['uniform_TV_bounds_float']=cert['uniform_TV_bounds_float']
        out.append(record)
        print(json.dumps({'task':task,'seed':seed,'numeric_replay':record['numeric_replay']}),flush=True)
    (ROOT/'results/compression_replay.json').write_text(json.dumps({'runs':out,
        'scope':'Teacher-forced numerical replay of fixed coordinate projection, no refitting or test-selected coordinates. The accompanying real-arithmetic certificate does not certify FP32 execution roundoff.'},indent=2))
if __name__=='__main__':main()
