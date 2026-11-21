"""Fresh-corpus pilot: scalar residual outputs and counting-stable multiletter heads.
Training is numerical. Exact stored-coefficient stability checks are separate.
No globally optimal AAK approximation solver is implemented.
"""
import argparse, copy, hashlib, itertools, json, math, os, platform, time
from pathlib import Path
import numpy as np
import torch
from torch import nn
from torch.nn import functional as F

ROOT = Path(__file__).resolve().parents[1]
ARMS = ['main_only', 'single', 'single_param', 'multi_path', 'multi_cert']
PERMS = np.array(list(itertools.permutations(range(3))), dtype=np.int64)
LOOKUP = {tuple(p):i for i,p in enumerate(PERMS)}
ENCODED_LOOKUP=np.zeros(27,dtype=np.int64)
for i,p in enumerate(PERMS): ENCODED_LOOKUP[p[0]*9+p[1]*3+p[2]]=i

def labels_s3(x):
    # token0 identity; token1 swaps positions0,1; token2 swaps positions1,2.
    p = np.tile(np.arange(3), (len(x),1))
    ys = np.zeros_like(x)
    for t in range(x.shape[1]):
        for a, idx in [(1,[1,0,2]), (2,[0,2,1])]:
            mask = x[:,t] == a
            p[mask] = p[mask][:,idx]
        ys[:,t] = ENCODED_LOOKUP[p[:,0]*9+p[:,1]*3+p[:,2]]
    return ys

@torch.jit.script
def recurrent(main_in: torch.Tensor, inject: torch.Tensor, routes: torch.Tensor,
              decay: torch.Tensor, matrices: torch.Tensor, write: torch.Tensor,
              has_residual: bool):
    batch, steps, width = main_in.shape
    m = torch.zeros((batch,width), dtype=main_in.dtype, device=main_in.device)
    heads, rank = write.shape
    r = torch.zeros((batch,heads,rank), dtype=main_in.dtype, device=main_in.device)
    states = torch.jit.annotate(list[torch.Tensor], [])
    for t in range(steps):
        m = decay*m + main_in[:,t,:]
        if has_residual:
            # hard forward routing; softmax straight-through derivative at training.
            a = torch.einsum('bm,mhij->bhij', routes[:,t,:], matrices)
            r = torch.einsum('bhij,bhj->bhi', a, r) + inject[:,t,:,None]*write[None,:,:]
            states.append(torch.cat([m,r.reshape(batch,heads*rank)],dim=-1))
        else:
            states.append(m)
    return torch.stack(states,dim=1)

class MemoryModel(nn.Module):
    def __init__(self, arm, vocab, outputs, task, hidden=64):
        super().__init__()
        self.arm, self.task = arm, task
        self.has_residual = arm != 'main_only'
        self.main_dim = 32 if self.has_residual else 64
        self.heads, self.rank = 4, 8
        self.modes = (3 if task=='order' else 4) if arm.startswith('multi_') else 1
        self.embedding = nn.Embedding(vocab,32)
        self.main_input = nn.Linear(32,self.main_dim)
        self.decay_raw = nn.Parameter(torch.linspace(-1.0,4.0,self.main_dim))
        if self.has_residual:
            self.injection = nn.Linear(32,self.heads)
            self.write = nn.Parameter(torch.randn(self.heads,self.rank)/math.sqrt(self.rank))
            base = .97*torch.eye(self.rank)[None].repeat(self.heads,1,1)
            base += .015*torch.randn_like(base)/math.sqrt(self.rank)
            self.raw_matrices = nn.Parameter(base[None].repeat(self.modes,1,1,1))
            self.router = nn.Linear(32,self.modes) if self.modes>1 and task=='text' else None
            self.head_read = nn.Parameter(torch.randn(self.heads,self.rank)/math.sqrt(self.rank))
            self.output_vectors = nn.Linear(self.heads,outputs,bias=False)
        self.norm = nn.LayerNorm(self.main_dim)
        self.readout = nn.Sequential(nn.Linear(self.main_dim,hidden),nn.GELU(),nn.Linear(hidden,outputs))

    def matrices(self):
        if hasattr(self,'frozen_matrices'): return self.frozen_matrices
        if self.arm=='multi_cert':
            # Per head: stack ALL mode matrices vertically. This caps sum R_a^T R_a.
            stack=self.raw_matrices.permute(1,0,2,3).reshape(self.heads,self.modes*self.rank,self.rank)
            sig=torch.linalg.matrix_norm(stack,ord=2,dim=(-2,-1))
            scale=torch.clamp(.99/sig,max=1.)
            return self.raw_matrices*scale[None,:,None,None]
        sig=torch.linalg.matrix_norm(self.raw_matrices,ord=2,dim=(-2,-1),keepdim=True)
        return self.raw_matrices*torch.clamp(.995/sig,max=1.)

    def forward(self,x,collapse=False,project_rank=None,return_parts=False):
        emb=self.embedding(x)
        main=torch.tanh(self.main_input(emb))
        decay=.999*torch.sigmoid(self.decay_raw)
        if self.has_residual:
            inject=torch.tanh(self.injection(emb))
            matrices=self.matrices()
            if self.modes==1: routes=torch.ones((*x.shape,1),device=x.device)
            elif self.task=='order': routes=F.one_hot(x,num_classes=self.modes).float()
            else:
                soft=torch.softmax(self.router(emb),dim=-1)
                hard=F.one_hot(soft.argmax(-1),num_classes=self.modes).float()
                routes=hard+soft-soft.detach() if self.training else hard
            if collapse:
                matrices=matrices.mean(0,keepdim=True)
                routes=torch.ones((*x.shape,1),device=x.device)
            rank=self.rank if project_rank is None else project_rank
            matrices=matrices[:,:,:rank,:rank]
            write=self.write[:,:rank]
        else:
            inject=torch.zeros((*x.shape,4),device=x.device)
            routes=torch.ones((*x.shape,1),device=x.device)
            matrices=torch.zeros((1,4,8,8),device=x.device)
            write=torch.zeros((4,8),device=x.device)
            rank=8
        h=recurrent(main,inject,routes,decay,matrices,write,self.has_residual)
        main_logits=self.readout(self.norm(h[:,:,:self.main_dim]))
        if self.has_residual:
            r=h[:,:,self.main_dim:].reshape(*x.shape,self.heads,rank)
            scalar=torch.einsum('bthr,hr->bth',r,self.head_read[:,:rank])
            residual_logits=self.output_vectors(scalar)
            logits=main_logits+residual_logits
        else:
            residual_logits=torch.zeros_like(main_logits)
            logits=main_logits
        if return_parts:
            return logits,main_logits,residual_logits,h,inject,routes
        return logits

def count_params(model): return sum(p.numel() for p in model.parameters())

def make_model(arm,vocab,outputs,task,seed):
    hidden=64
    if arm=='single_param':
        torch.manual_seed(seed)
        target=count_params(MemoryModel('multi_cert',vocab,outputs,task))
        candidates=[]
        for h in range(1,201):
            torch.manual_seed(seed)
            n=count_params(MemoryModel('single_param',vocab,outputs,task,h))
            candidates.append((abs(target-n),h))
        hidden=min(candidates)[1]
    torch.manual_seed(seed)
    return MemoryModel(arm,vocab,outputs,task,hidden)

def load_text():
    texts={k:(ROOT/'data'/f'{k}.txt').read_text() for k in ['train','validation','test']}
    vocab=['<UNK>']+sorted(set(texts['train']))
    mapping={c:i for i,c in enumerate(vocab)}
    encoded={k:np.array([mapping.get(c,0) for c in text],dtype=np.int64) for k,text in texts.items()}
    return encoded,vocab

def text_windows(data,starts,length):
    idx=starts[:,None]+np.arange(length+1)[None,:]
    z=data[idx]
    return torch.from_numpy(z[:,:-1].copy()),torch.from_numpy(z[:,1:].copy())

@torch.no_grad()
def evaluate_text(model,data,collapse=False,project_rank=None):
    model.eval()
    starts=np.arange(0,len(data)-64,64)
    total,n=0.,0
    for off in range(0,len(starts),64):
        x,y=text_windows(data,starts[off:off+64],64)
        z=model(x,collapse,project_rank)
        total += F.cross_entropy(z.reshape(-1,z.shape[-1]),y.reshape(-1),reduction='sum').item()
        n+=y.numel()
    return {'bpc':total/n/math.log(2),'perplexity':math.exp(total/n),'tokens':n}

@torch.no_grad()
def evaluate_order(model,length,seed,n=1024,paired=False,collapse=False,project_rank=None):
    model.eval(); rng=np.random.default_rng(seed)
    x=rng.integers(0,3,size=(n,length),dtype=np.int64)
    if paired:
        # swap immediately after8 shared events; equal long suffix afterwards.
        x[:,8:10]=[1,2]
        other=x.copy(); other[:,8:10]=[2,1]
        x=np.concatenate([x,other])
    y=labels_s3(x)[:,-1]
    pred=[]; losses=[]
    for off in range(0,len(x),128):
        z=model(torch.from_numpy(x[off:off+128]),collapse,project_rank)[:,-1,:]
        losses.append(F.cross_entropy(z,torch.from_numpy(y[off:off+128]),reduction='sum').item())
        pred.extend(z.argmax(-1).tolist())
    pred=np.array(pred); correct=pred==y
    out={'accuracy':float(correct.mean()),'loss':sum(losses)/len(x),'n':len(x),
         'label_counts':np.bincount(y,minlength=6).tolist()}
    if paired:
        assert np.all(y[:n]!=y[n:])
        out.update(pair_both_correct=float((correct[:n]&correct[n:]).mean()),
                   pair_predictions_different=float((pred[:n]!=pred[n:]).mean()))
    return out

@torch.no_grad()
def diagnostics(model,vocab):
    model.eval(); out={}
    if model.has_residual:
        a=model.matrices().double().numpy()
        rhos=[]; comm=[]
        for h in range(model.heads):
            cp=sum(np.kron(a[j,h],a[j,h]) for j in range(model.modes))
            rhos.append(float(max(abs(np.linalg.eigvals(cp)))))
            comm.extend(float(np.linalg.norm(a[j,h]@a[k,h]-a[k,h]@a[j,h]))
                        for j in range(model.modes) for k in range(j))
        out={'counting_cp_spectral_radius_per_head':rhos,
             'counting_stability_sufficient_condition':all(v<1 for v in rhos),
             'mean_commutator_frobenius':float(np.mean(comm)) if comm else 0.}
        if model.router is not None:
            x=torch.arange(len(vocab)); routes=model.router(model.embedding(x)).argmax(-1).tolist()
            out['character_routes']={str(vocab[i]):routes[i] for i in range(len(vocab))}
    # Batch1 total forward timing at fixed64 context, not a fused decoding benchmark.
    x=torch.zeros((1,64),dtype=torch.long)
    for _ in range(3): model(x)
    times=[]
    for _ in range(7):
        s=time.perf_counter(); model(x); times.append((time.perf_counter()-s)/64)
    out['cpu_us_per_token_batch1_window64']=float(np.median(times)*1e6)
    return out

def train(args):
    torch.set_num_threads(1); torch.set_num_interop_threads(1)
    torch.use_deterministic_algorithms(True)
    np.random.seed(args.seed)
    data,vocab=load_text() if args.task=='text' else (None,['identity','swap01','swap12'])
    outputs=len(vocab) if args.task=='text' else 6
    model=make_model(args.arm,len(vocab),outputs,args.task,args.seed)
    opt=torch.optim.AdamW(model.parameters(),lr=.003,weight_decay=.01)
    rng=np.random.default_rng(10000+args.seed)
    stem=f'{args.task}_{args.arm}_seed{args.seed}'
    result_path=ROOT/'results'/f'{stem}.json'
    trace=[]; best=float('inf'); best_state=None; best_step=0
    start=time.perf_counter()
    for step in range(args.steps+1):
        if step%100==0 or step==args.steps:
            val=evaluate_text(model,data['validation']) if args.task=='text' else evaluate_order(model,32,17001,n=512)
            score=val['bpc'] if args.task=='text' else val['loss']
            trace.append({'step':step,'validation':val,'elapsed_seconds':time.perf_counter()-start})
            if score<best:
                best=score; best_step=step; best_state=copy.deepcopy(model.state_dict())
            print(json.dumps({'run':stem,**trace[-1]}),flush=True)
        if step==args.steps: break
        model.train()
        if args.task=='text':
            starts=rng.integers(0,len(data['train'])-64,size=32)
            x,y=text_windows(data['train'],starts,64)
        else:
            raw=rng.integers(0,3,size=(64,32),dtype=np.int64)
            x,y=torch.from_numpy(raw),torch.from_numpy(labels_s3(raw))
        for group in opt.param_groups:
            group['lr']=.003 if step<.8*args.steps else .0006
        opt.zero_grad(set_to_none=True)
        logits=model(x)
        loss=F.cross_entropy(logits.reshape(-1,outputs),y.reshape(-1))
        if not torch.isfinite(loss): raise RuntimeError('nonfinite loss')
        loss.backward(); nn.utils.clip_grad_norm_(model.parameters(),1.0); opt.step()
    train_secs=time.perf_counter()-start
    model.load_state_dict(best_state)
    if model.has_residual: model.frozen_matrices=model.matrices().detach().clone()
    if args.task=='text':
        test=evaluate_text(model,data['test'])
        collapse=evaluate_text(model,data['test'],True) if args.arm.startswith('multi_') else None
    else:
        test={str(L):evaluate_order(model,L,18000+L) for L in [32,64,128]}
        test['pairs32']=evaluate_order(model,32,19001,n=512,paired=True)
        test['pairs96']=evaluate_order(model,96,19002,n=512,paired=True)
        collapse=evaluate_order(model,32,18032,collapse=True) if args.arm.startswith('multi_') else None
    torch.save({'state_dict':best_state,'args':vars(args),'vocab':vocab,'effective_matrices':model.matrices().detach().clone() if model.has_residual else None},ROOT/'checkpoints'/f'{stem}.pt')
    out={'task':args.task,'arm':args.arm,'seed':args.seed,'steps':args.steps,'selected_step':best_step,
         'parameters':count_params(model),'parameter_bytes_fp32':4*count_params(model),
         'state_scalars':64,'state_bytes_fp32':256,'main_state_scalars':model.main_dim,
         'residual_state_scalars':32 if model.has_residual else 0,
         'residual_head_dim':8 if model.has_residual else 0,'transition_modes':model.modes,
         'readout_hidden':model.readout[0].out_features,'residual_readout':'scalar_per_head_additive_logits','validation_best':best,'test':test,
         'collapsed_transition_test':collapse,'training_seconds':train_secs,
         'trace':trace,'diagnostics':diagnostics(model,vocab),
         'torch_version':torch.__version__,'python':platform.python_version()}
    result_path.write_text(json.dumps(out,indent=2,ensure_ascii=False))
    print('DONE '+stem,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('--task',choices=['text','order'],required=True)
    p.add_argument('--arm',choices=ARMS,required=True); p.add_argument('--seed',type=int,default=101)
    p.add_argument('--steps',type=int,default=1800)
    train(p.parse_args())
