"""Positive, normalized finite-belief machines. All acceptance tests use Fraction.
Synthetic four-state mixture reduction is separate from the fresh trained task.
The trained task is character prediction on the unused end of the same book.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import json,math,hashlib,re,unicodedata,ctypes,subprocess
from pathlib import Path
from fractions import Fraction as F
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
def ser(x):
    if isinstance(x,(list,tuple)):return [ser(z) for z in x]
    if isinstance(x,dict):return {k:ser(v) for k,v in x.items()}
    return str(x) if isinstance(x,F) else x
def nearest(n,d):
    q,r=divmod(n,d);return q+int(2*r>d or (2*r==d and q%2))
def quantize(rows,Q=24):
    D=2**Q;out=[]
    for row in rows:
        assert sum(row)==1 and min(row)>0
        z=[int(x*D) for x in row];rem=D-sum(z)
        for j in sorted(range(len(row)),key=lambda j:row[j]*D-z[j],reverse=True)[:rem]:z[j]+=1
        assert sum(z)==D and min(z)>0
        out.append(z)
    return out
def certificate(rows,delta=F(0),tau_source=None):
    ints=quantize(rows);D=2**24;K=(len(rows[0])-1)//2
    qrows=[[F(x,D) for x in r] for r in ints]
    tau=min(r[-1] for r in qrows);tau0=min(r[-1] for r in rows)
    dr=max(sum(abs(a-b) for a,b in zip(x,y))/2 for x,y in zip(rows,qrows))
    q=F(0)
    for a in range(K):
        aa,ab=qrows[0][2*a:2*a+2];ba,bb=qrows[1][2*a:2*a+2]
        q=max(q,abs(aa*bb-ab*ba)/min(aa+ab,ba+bb)**2)
    assert q<1
    events=[[sum(r[2*a:2*a+2]) for a in range(K)]+[r[-1]] for r in qrows]
    C=sum(abs(a-b) for a,b in zip(*events))/2
    for f in range(2,23):
        E=F(1,2**(f+1))/(1-q);b=delta+dr/tau0+C*E/tau
        if b<=F(1,100):break
    assert b<=F(1,100)
    # Conservative KL after an explicit geometric mixture; logarithms here numeric.
    t=min(tau0,tau_source or tau0);r=1-float(t);d=float(b)
    N=max(0,math.ceil(math.log(d)/math.log(r))) if r else 0
    psi=d*N+(r**(N+1))/float(t)
    kl=-math.log1p(-d)+d*math.log(2)+math.log(2*K)*psi
    return {'rows_exact':ser(rows),'coefficient_integers':ints,'Q':24,'alphabet_size':K,
      'contraction':str(q),'event_lipschitz':str(C),'stopping_min':str(tau),
      'coefficient_row_TV':str(dr),'coefficient_word_TV':str(dr/tau0),
      'reduction_TV':str(delta),'belief_fractional_bits':f,'state_bits':f+1,
      'state_bytes':math.ceil((f+1)/8),'belief_error':str(E),'machine_word_TV_bound':str(C*E/tau),
      'total_word_TV_bound':str(b),'KL_geometric_mixture_upper_numeric':kl,
      'KL_note':'Natural logarithms; applies only after qplus=(1-d)q+d*u, u(w)=1/2*(1/(2K))^len(w). Mixture implementation not supplied.',
      'tail_min_stopping':str(t),'tail_cutoff':N}
def pstep(j,a,ints,G):
    n=j*ints[0][2*a]+(G-j)*ints[1][2*a]
    d=n+j*ints[0][2*a+1]+(G-j)*ints[1][2*a+1]
    return nearest(n*G,d)
def machine_checks(cert,seed,exhaustive=True):
    ints=cert['coefficient_integers'];G=2**cert['belief_fractional_bits'];K=cert['alphabet_size']
    A=np.array(ints,np.uint32);ptr=A.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32))
    lib=ctypes.CDLL(str(ROOT/'code/belief.so'));fn=lib.belief_step
    fn.argtypes=[ctypes.c_uint32,ctypes.c_int,ctypes.c_uint32,ctypes.c_int,ctypes.POINTER(ctypes.c_uint32)];fn.restype=ctypes.c_uint32
    packed=lib.belief_step_packed
    packed.argtypes=[ctypes.POINTER(ctypes.c_uint8),ctypes.c_int,ctypes.c_int,ctypes.c_uint32,ctypes.c_int,ctypes.POINTER(ctypes.c_uint32)]
    def actual(j,a):
        raw=bytearray(j.to_bytes(cert['state_bytes'],'little'));buf=(ctypes.c_uint8*len(raw)).from_buffer(raw)
        assert packed(buf,len(raw),int(a),G,K,ptr)==0
        return int.from_bytes(raw,'little')
    states={G//2};todo=[G//2];steps=0
    while todo:
        j=todo.pop()
        # Emitted next-character and termination probabilities sum to exactly one.
        assert sum(j*ints[0][v]+(G-j)*ints[1][v] for v in range(2*K+1))==G*2**24
        for a in range(K):
            nxt=pstep(j,a,ints,G);assert nxt==actual(j,a);steps+=1
            if nxt not in states:states.add(nxt);todo.append(nxt)
    rng=np.random.default_rng(seed);j=G//2;x=.5;maxerr=0;Afloat=np.array(ints)/2**24
    for a in rng.integers(0,K,100000):
        j2=actual(j,int(a));assert j2==pstep(j,int(a),ints,G);j=j2;steps+=1
        v=x*Afloat[0,2*a:2*a+2]+(1-x)*Afloat[1,2*a:2*a+2];x=v[0]/sum(v)
        maxerr=max(maxerr,abs(x-j/G))
    assert maxerr<=float(F(cert['belief_error']))+1e-12
    return {'C_python_steps':steps,'reachable_states_complete':len(states),'max_observed_belief_error':maxerr,
      'normalization_all_reachable_states':True,'precision':'Packed persistent byte array of state_bytes, uint32 decoded temporary and uint128 arithmetic temporaries'}
def rows_from(T,E,tau):
    rows=[]
    for i in range(2):
        row=[]
        for a in range(E.shape[1]):
            for j in range(2):row.append((1-tau)*E[i,a]*T[i,j])
        row.append(tau);rows.append(row)
    return rows
def synthetic():
    out=[]
    for tau in [F(1,20),F(1,10),F(1,5)]:
      for seed in range(4):
        T=np.array([[F(9+seed,20),F(11-seed,20)],[F(8+seed,20),F(12-seed,20)]],object)
        E=np.array([[F(3,5),F(2,5)],[F(2,5),F(3,5)]],object)
        rows=rows_from(T,E,tau);c=certificate(rows,F(1,1000),tau)
        # Source is exactly (999/1000)*core + (1/1000)*other.
        To=np.array([[F(3,4),F(1,4)],[F(1,3),F(2,3)]],object)
        Eo=np.array([[F(4,5),F(1,5)],[F(1,4),F(3,4)]],object)
        c['other_rows']=ser(rows_from(To,Eo,tau));c['initial_belief']='1/2'
        c['id']=f'pfa_tau{tau.denominator}_s{seed}';c['checks']=machine_checks(c,seed)
        out.append(c);print(c['id'],c['total_word_TV_bound'],flush=True)
    return out
def fresh_data():
    path=ROOT/'prior_evidence/certificate_pilot/data/pg15058_raw.txt';raw=path.read_bytes()
    s=raw.decode('utf-8-sig').replace('\r\n','\n').replace('\r','\n');start=s.index('Ma bonne, ma bonne, dit un jour Sophie')
    body=unicodedata.normalize('NFC',s[start:s.index('*** END OF THE PROJECT GUTENBERG EBOOK',start)])
    paras=[re.sub(r'\s+',' ',p).strip() for p in re.split(r'\n\s*\n',body) if p.strip()]
    full=unicodedata.normalize('NFC','\n\n'.join(paras).lower());rest=full[180000:]
    train=rest[:8000];val=rest[8000:10000];test=rest[10000:]
    vocab=sorted(set(train));idx={a:i for i,a in enumerate(vocab)};unk=len(vocab)
    encode=lambda t:np.array([idx.get(a,unk) for a in t])
    old={full[i:i+64] for i in range(180000-63)};tr={train[i:i+64] for i in range(len(train)-63)}
    meta={'source_sha256':hashlib.sha256(raw).hexdigest(),'full_chars':len(full),'new_region':[180000,len(full)],
      'vocab_from_train_only':vocab+['<UNK>'],'train_chars':len(train),'validation_chars':len(val),'test_chars':len(test),
      'overlap64_test_prior':sum(test[i:i+64] in old for i in range(len(test)-63)),
      'overlap64_test_train':sum(test[i:i+64] in tr for i in range(len(test)-63)),
      'test_sha256':hashlib.sha256(test.encode()).hexdigest(),'freshness':'Unused suffix of the same book in the supplied pilot. No claim of independent author or universal historical freshness.',
      'protocol':'5 seeds, 30 projected Baum-Welch iterations, windows 64; all settings fixed before test evaluation; validation not used for model selection.'}
    return encode(train),encode(val),encode(test),meta
def train_hmm(seq,K,seed):
    rng=np.random.default_rng(seed);X=seq[:len(seq)//64*64].reshape(-1,64);B,L=X.shape
    base=(np.bincount(seq,minlength=K)+1)/(len(seq)+K)
    E=np.stack([base*(1+rng.uniform(-.1,.1,K)),base*(1+rng.uniform(-.1,.1,K))]);E/=E.sum(1)[:,None]
    T=np.array([[.55,.45],[.45,.55]])
    for epoch in range(30):
        alpha=np.empty((B,L,2));sc=np.empty((B,L));v=np.ones((B,2))/2
        for t in range(L):
            u=v*E[:,X[:,t]].T;sc[:,t]=u.sum(1);alpha[:,t]=u/sc[:,t,None];v=alpha[:,t]@T
        beta=np.ones((B,L,2))
        for t in range(L-2,-1,-1):beta[:,t]=(beta[:,t+1]*E[:,X[:,t+1]].T)@T.T/sc[:,t+1,None]
        gamma=alpha*beta;gamma/=gamma.sum(2)[:,:,None]
        counts=np.ones((2,K))*.1
        for i in range(2):np.add.at(counts[i],X.ravel(),gamma[:,:,i].ravel())
        candidate=counts/counts.sum(1)[:,None];ratio=np.clip(candidate/base,.75,1.25)
        E=ratio*base;E/=E.sum(1)[:,None]
        tc=np.ones((2,2))*.1
        for t in range(L-1):
            xi=alpha[:,t,:,None]*T[None,:,:]*(E[:,X[:,t+1]].T*beta[:,t+1])[:,None,:]
            xi/=xi.sum((1,2))[:,None,None];tc+=xi.sum(0)
        first=np.clip(tc[:,0]/tc.sum(1),.4,.6);T=np.stack([first,1-first],axis=1)
    # Rationalize normalized rows exactly before defining the reference model.
    rational=lambda row: [F(str(float(x))) for x in row]
    Tr=[rational(row) for row in T];Er=[rational(row) for row in E]
    Tr=[[x/sum(row) for x in row] for row in Tr];Er=[[x/sum(row) for x in row] for row in Er]
    return np.array(Tr,object),np.array(Er,object),base
def evaluate(seq,T,E,cert,base):
    # Conditional on continuation; stochastic stopping is external with fixed tau=.1.
    T=np.array(T,float);E=np.array(E,float);v=np.array([.5,.5]);j=2**(cert['belief_fractional_bits']-1)
    G=2**cert['belief_fractional_bits'];ints=cert['coefficient_integers'];D=2**24;logref=logq=logiid=0.;mx=0
    for a in seq:
        p=v@E;logref-=math.log2(p[a]);u=v*E[:,a];v=(u/u.sum())@T
        nums=np.array([j*ints[0][2*b]+(G-j)*ints[1][2*b]+j*ints[0][2*b+1]+(G-j)*ints[1][2*b+1] for b in range(E.shape[1])],float)
        pq=nums/nums.sum();logq-=math.log2(pq[a]);mx=max(mx,.5*abs(p-pq).sum());j=pstep(j,int(a),ints,G)
        logiid-=math.log2(base[a])
    return {'reference_BPC':logref/len(seq),'quantized_BPC':logq/len(seq),'iid_BPC':logiid/len(seq),'observed_max_conditional_TV':mx,'characters':len(seq)}
def run():
    subprocess.run(['gcc','-std=c11','-O2','-shared','-fPIC',str(ROOT/'code/belief.c'),'-o',str(ROOT/'code/belief.so')],check=True)
    syn=synthetic();tr,va,te,meta=fresh_data();trained=[];K=len(meta['vocab_from_train_only'])
    for seed in [101,102,103,104,105]:
        T,E,base=train_hmm(tr,K,seed);c=certificate(rows_from(T,E,F(1,10)))
        c['id']=f'fresh_char_seed{seed}';c['checks']=machine_checks(c,seed)
        c['test']=evaluate(te,T,E,c,base);c['validation']=evaluate(va,T,E,c,base)
        c['reference_T']=ser(T.tolist());c['reference_E']=ser(E.tolist());trained.append(c)
        print(c['id'],c['test'],flush=True)
    out={'synthetic':syn,'trained':trained,'data':meta,
      'scope':'Probabilistic branch separate from signed Hankel rank-four pilot. Word-law TV proof uses proper stopping; it is not a uniform next-token TV theorem. Training is a new small char-HMM experiment, not a rerun of the earlier neural model.'}
    (ROOT/'results/probability.json').write_text(json.dumps(out,indent=2,ensure_ascii=False))
if __name__=='__main__':run()
