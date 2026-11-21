"""Constructive rank-two finite-bit pilot. No training, no claim of global optimum.
All acceptance bounds in this file are rational; floating point is diagnostic.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
import json, math, time, itertools, hashlib, platform
from pathlib import Path
from fractions import Fraction as F
import numpy as np
import scipy.linalg as la
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results'; MOD=ROOT/'models'
OUT.mkdir(exist_ok=True);MOD.mkdir(exist_ok=True)
Q=24
TOL=F(1,100)

def sm(x): return sp.Matrix(x)
def frac(x): return F(int(sp.numer(x)),int(sp.denom(x)))
def arr(M): return np.array(M,dtype=float)
def ser(M): return [[str(x) for x in row] for row in np.array(M,dtype=object).tolist()]
def scalar(x): return str(x)
def rnd(x):
    x=frac(x); q,r=divmod(x.numerator,x.denominator)
    return q+int(2*r>x.denominator or (2*r==x.denominator and q%2))
def dyad(M,q=Q): return M.applyfunc(lambda x:sp.Rational(rnd(x*2**q),2**q))
def psd(M,strict=False):
    M=sp.Matrix(M)
    if M!=M.T:return False
    while M.rows:
        p=M[0,0]
        if p<0 or (strict and p==0):return False
        if p==0:
            if any(M[0,j]!=0 for j in range(M.cols)):return False
            M=M[1:,1:]
        else:M=M[1:,1:]-M[1:,:1]*M[:1,1:]/p
    return True
def gram(As,B):
    n=As[0].rows
    L=sp.eye(n*n)-sum((sp.kronecker_product(a,a) for a in As),sp.zeros(n*n))
    v=L.inv(method='DM')*sp.Matrix(list(B))
    return sp.Matrix(n,n,list(v))
def hsq(As,alpha,beta):
    P=gram(As,beta*beta.T)
    Qg=gram([a.T for a in As],alpha.T*alpha)
    assert psd(P) and psd(Qg)
    return sp.trace(P*Qg),P,Qg
def sqrt_up(x,bits=48):
    x=frac(x)
    n=x.numerator*2**(2*bits)//x.denominator
    r=math.isqrt(n)
    while F(r*r,2**(2*bits))<x:r+=1
    return F(r,2**bits)
def hankel_lower_diagnostic(As,alpha,beta,depth=4):
    words=[()]
    for l in range(1,depth+1):words.extend(itertools.product(range(len(As)),repeat=l))
    rows=[];cols=[]
    for w in words:
        m=np.eye(As[0].rows)
        for a in w:m=m@arr(As[a])
        rows.append((arr(alpha)@m).ravel());cols.append((m@arr(beta)).ravel())
    H=np.array(rows)@np.array(cols).T
    return float(np.linalg.svd(H,compute_uv=False)[2])
def rank_certificate(As,alpha,beta,n):
    words=[()]+[w for l in range(1,n) for w in itertools.product(range(2),repeat=l)]
    rows=[];cols=[]
    for w in words:
        m=sp.eye(n)
        for a in w:m=m*As[a]
        rows.append(list(alpha*m));cols.append(list(m*beta))
    R=sp.Matrix(rows); C=sp.Matrix(cols).T
    ri=list(R.T.rref()[1])[:n];ci=list(C.rref()[1])[:n]
    assert len(ri)==len(ci)==n
    dr=R[ri,:].det();dc=C[:,ci].det()
    assert dr and dc
    return {'row_words':[words[i] for i in ri],'column_words':[words[i] for i in ci],
            'row_determinant':str(dr),'column_determinant':str(dc)}
def instance(family,margin,seed,k=2):
    rng=np.random.default_rng(seed)
    q=sp.Rational([35,45,55,67][margin],100)
    # Two noncommuting substochastic rational transition matrices.
    eps=sp.Rational([1,2,4][seed%3],128)
    a0=sm([[sp.Rational(4,5),sp.Rational(1,10)], [eps,sp.Rational(3,5)]])*q
    a1=sm([[sp.Rational(1,2),sp.Rational(1,5)], [2*eps,sp.Rational(4,5)]])*q
    if family=='boundary':
        a0=sm([[sp.Rational(94,100),sp.Rational(1,100)],[eps/4,sp.Rational(93,100)]])*q
        a1=sm([[sp.Rational(92,100),sp.Rational(2,100)],[eps/2,sp.Rational(94,100)]])*q
    A=[a0,a1];alpha=sm([[sp.Rational(1,2),sp.Rational(1,4)]])
    beta=sm([sp.Rational(1,2),sp.Rational([1,4,16][seed%3],256)])
    if family=='nonnormal':
        S=sm([[1,sp.Rational(7,4)],[0,1]])
        A=[S.inv()*a*S for a in A];alpha=alpha*S;beta=S.inv()*beta
    # Weak, reachable and observable extra block gives a genuine rank-four source.
    R=[sm([[sp.Rational(1,5),sp.Rational(1,11)],[0,sp.Rational(2,7)]]),
       sm([[sp.Rational(1,7),0],[sp.Rational(1,13),sp.Rational(1,4)]])]
    ar=sm([[sp.Rational(1,2),sp.Rational(1,3)]])
    br=sm([sp.Rational(seed%3+1,10000),sp.Rational(1,13000)])
    As=[sp.diag(A[i],R[i]) for i in range(2)]
    af=alpha.row_join(ar);bf=beta.col_join(br)
    # Normalization stipulated by the frozen protocol, prior to reduction/quantization.
    sq,_,_=hsq(As,af,bf)
    while sq>1:
        beta=beta/2;br=br/2;bf=bf/2;sq=sq/4
    mix=sm([[1,1,0,0],[0,1,1,0],[0,0,1,1],[0,0,0,1]])
    dense=[mix.inv()*a*mix for a in As];ad=af*mix;bd=mix.inv()*bf
    return A,alpha,beta,R,ar,br,As,af,bf,dense,ad,bd
def balanced_basis(A,alpha,beta):
    P=np.array(gram(A,beta*beta.T),float);O=np.array(gram([a.T for a in A],alpha.T*alpha),float)
    # Standard balanced similarity from reachability and observability.
    L=la.cholesky(P,lower=True)
    vals,U=la.eigh(L.T@O@L)
    vals=np.maximum(vals,1e-30)
    S=L@U@np.diag(vals**(-.25))
    S=sp.Matrix([[sp.Rational(int(round(v*65536)),65536) for v in row] for row in S])
    if S.det()==0:raise ValueError('rounded singular balancing transform')
    return S
def enc_frame(A,alpha,beta,S):
    Ai=[S.inv()*a*S for a in A];al=alpha*S;be=S.inv()*beta
    Aq=[dyad(a) for a in Ai];aq=dyad(al);bq=dyad(be)
    k=be.rows
    # Common nonnegative majorant and an exact resolvent witness.
    M=sm(k,k,lambda i,j:max(abs(a[i,j]) for a in Aq)) if False else sp.Matrix(k,k,lambda i,j:max(abs(a[i,j]) for a in Aq))
    K=sp.eye(k)-M.T
    if K.det()==0:return None
    inv=K.inv(method='DM')
    if any(v<0 for v in inv):return None
    R0=inv*aq.T.applyfunc(abs)
    if any(v<=0 for v in R0):return None
    exp=[int(math.ceil(math.log2(2*float(v)+1e-5))) for v in R0]
    J=[sp.diag(Ai[i],Aq[i]) for i in range(2)]
    ja=al.row_join(-aq);jb=be.col_join(bq)
    sq,P,O=hsq(J,ja,jb);coef=sqrt_up(sq)
    X=gram(Aq,sp.eye(k))
    assert psd(X,True)
    return {'A':Aq,'alpha':aq,'beta':bq,'M':M,'inv':inv,'exp':exp,
            'coef':coef,'coef_sq':sq,'coef_P':P,'coef_Q':O,'stable_X':X,'S':S,
            'unrounded_A':Ai,'unrounded_alpha':al,'unrounded_beta':be}
def allocation(frame,bits):
    k=len(bits);ex=frame['exp'];inv=frame['inv']
    steps=[sp.Rational(2)**(ex[i]-bits[i]+1) for i in range(k)]
    half=sm([h/2 for h in steps]);aq=frame['alpha']
    init=[rnd(aq[i]/steps[i]) for i in range(k)]
    actual=sm([init[i]*steps[i] for i in range(k)])
    R=inv*(actual.applyfunc(abs)+half)
    limits=[(2**(bits[i]-1)-1)*steps[i] for i in range(k)]
    if any(R[i]>limits[i] or init[i]<-2**(bits[i]-1) or init[i]>=2**(bits[i]-1) for i in range(k)):return None
    E=inv*half
    assert all(abs(actual[i]-aq[i])<=half[i] for i in range(k))
    err=(frame['beta'].T.applyfunc(abs)*E)[0]
    # Bound integer accumulator width for the exact implementation.
    shifts=[ex[i]-bits[i]+1 for i in range(k)]
    need=0
    for j in range(k):
        pows=[shifts[i]-Q for i in range(k)]
        minpow=min(pows)
        bound=sum(abs(int(frame['A'][a][i,j]*2**Q))*2**(bits[i]-1)*2**(pows[i]-minpow)
                  for i in range(k) for a in [0])
        bound=max(bound,sum(abs(int(frame['A'][1][i,j]*2**Q))*2**(bits[i]-1)*2**(pows[i]-minpow) for i in range(k)))
        need=max(need,int(bound).bit_length()+1)
    if need>126:return None
    return {'bits':list(bits),'state_bits':sum(bits),'state_bytes':math.ceil(sum(bits)/8),
            'steps':[str(s) for s in steps],'shifts':shifts,'initial_integers':init,
            'invariant_R':ser(R),'error_E':ser(E),'execution_bound':str(err),
            'accumulator_signed_bits_bound':need}
def run():
    start=time.time();models=[];records=[]
    for family in ['near_rank2','nonnormal','boundary']:
      for margin in range(4):
       for seed in [101,102,103]:
        name=f'{family}_m{margin}_s{seed}'
        if (MOD/(name+'.json')).exists():
            model=json.loads((MOD/(name+'.json')).read_text());models.append(model)
            u,m=model['best']['uniform'],model['best']['mixed']
            records.append({'id':name,'family':family,'model_bound':float(F(model['model_error_bound'])),
              'uniform_bits':u['state_bits'] if u else None,'mixed_bits':m['state_bits'] if m else None,
              'uniform_bytes':u['state_bytes'] if u else None,'mixed_bytes':m['state_bytes'] if m else None,
              'byte_saving_percent':100*(1-m['state_bytes']/u['state_bytes']) if u and m else None,
              'mixed_bound':float(F(m['total_bound'])) if m else None,
              'uniform_basis':u['basis'] if u else None,'mixed_basis':m['basis'] if m else None})
            continue
        A,alpha,beta,R,ar,br,As,af,bf,dense,ad,bd=instance(family,margin,seed)
        rank=rank_certificate(dense,ad,bd,4)
        sqr,Pr,Qr=hsq(R,ar,br);model_bound=sqrt_up(sqr)
        sqf,Pf,Qf=hsq(As,af,bf);assert sqf<=1
        source_X=gram(As,sp.eye(4));assert psd(source_X,True)
        frames={};candidates={}
        for label,S in [('native',sp.eye(2)),('balanced',balanced_basis(A,alpha,beta))]:
            fr=enc_frame(A,alpha,beta,S)
            if fr is None:
                frames[label]={'accepted':False,'reason':'entrywise comparison resolvent not nonnegative'};continue
            frameout={'accepted':True,'A':[ser(a) for a in fr['A']],'alpha':ser(fr['alpha']),'beta':ser(fr['beta']),
                'M':ser(fr['M']),'inverse_comparison':ser(fr['inv']),'range_exponents':fr['exp'],
                'S':ser(S),'coefficient_fractional_bits':Q,'coefficient_error_bound':str(fr['coef']),
                'coefficient_error_HS_squared':str(fr['coef_sq']),
                'coefficient_error_P':ser(fr['coef_P']),'coefficient_error_Q':ser(fr['coef_Q']),
                'stability_X':ser(fr['stable_X']),
                'unrounded_A':[ser(a) for a in fr['unrounded_A']],
                'unrounded_alpha':ser(fr['unrounded_alpha']),'unrounded_beta':ser(fr['unrounded_beta'])}
            allc=[]
            for p0 in range(3,33):
              for p1 in range(3,33):
                if p0+p1>64:continue
                c=allocation(fr,[p0,p1])
                if c is None:continue
                total=model_bound+fr['coef']+F(c['execution_bound'])
                c['total_bound']=str(total)
                c['passes_split']=model_bound<=F(5,1000) and fr['coef']<=F(2,1000) and F(c['execution_bound'])<=F(3,1000)
                c['passes_total']=total<=TOL
                c['basis']=label
                allc.append(c)
            frameout['allocation_count']=len(allc)
            for policy in ['uniform','mixed']:
                valid=[c for c in allc if c['passes_split'] and (policy=='mixed' or c['bits'][0]==c['bits'][1])]
                best=min(valid,key=lambda c:(c['state_bytes'],c['state_bits'],F(c['execution_bound']))) if valid else None
                frameout['best_'+policy]=best
            frameout['budget_frontier']={str(b):min((c for c in allc if c['state_bits']<=b),key=lambda c:F(c['total_bound']),default=None) for b in [8,12,16,24,32,48,64]}
            frames[label]=frameout;candidates[label]=allc
        best={}
        for policy in ['uniform','mixed']:
            opts=[fr['best_'+policy] for fr in frames.values() if fr['accepted'] and fr['best_'+policy]]
            best[policy]=min(opts,key=lambda c:(c['state_bytes'],c['state_bits'],F(c['total_bound']))) if opts else None
        model={'id':name,'family':family,'margin_index':margin,'seed':seed,'source_n':4,'target_k':2,
            'source_A':[ser(a) for a in dense],'source_alpha':ser(ad),'source_beta':ser(bd),
            'source_hidden_A':[ser(a) for a in As],'source_hidden_alpha':ser(af),'source_hidden_beta':ser(bf),
            'source_similarity':ser(sm([[1,1,0,0],[0,1,1,0],[0,0,1,1],[0,0,0,1]])),
            'source_P':ser(Pf),'source_Q':ser(Qf),'source_HS_squared':str(sqf),'source_stability_X':ser(source_X),
            'rank_certificate':rank,'reduced_A':[ser(a) for a in A],'reduced_alpha':ser(alpha),'reduced_beta':ser(beta),
            'residual_A':[ser(a) for a in R],'residual_alpha':ser(ar),'residual_beta':ser(br),
            'residual_P':ser(Pr),'residual_Q':ser(Qr),'model_error_HS_squared':str(sqr),'model_error_bound':str(model_bound),
            'sigma3_finite_section_diagnostic':hankel_lower_diagnostic(As,af,bf),'frames':frames,'best':best,
            'global_optimality':'NOT_CERTIFIED; explicit block reduction, finite-section SVD is numeric only'}
        (MOD/(name+'.json')).write_text(json.dumps(model,indent=2))
        models.append(model)
        u,m=best['uniform'],best['mixed']
        rec={'id':name,'family':family,'model_bound':float(model_bound),'uniform_bits':u['state_bits'] if u else None,
             'mixed_bits':m['state_bits'] if m else None,'uniform_bytes':u['state_bytes'] if u else None,'mixed_bytes':m['state_bytes'] if m else None,
             'byte_saving_percent':100*(1-m['state_bytes']/u['state_bytes']) if u and m else None,
             'mixed_bound':float(F(m['total_bound'])) if m else None,
             'uniform_basis':u['basis'] if u else None,'mixed_basis':m['basis'] if m else None}
        records.append(rec);print(json.dumps(rec),flush=True)
        (OUT/'pilot_progress.json').write_text(json.dumps(records,indent=2))
    summary={'instances':len(models),'records':records,'elapsed_seconds':time.time()-start,
             'scope':'constructed rational sources, known core plus weak block; no learned-model compression or optimality claim',
             'gate_25_percent_by_family':{fam:all(r['byte_saving_percent'] is not None and r['byte_saving_percent']>=25 for r in records if r['family']==fam) for fam in ['near_rank2','nonnormal','boundary']}}
    (OUT/'pilot_summary.json').write_text(json.dumps(summary,indent=2))
    print('COMPLETE',time.time()-start,flush=True)

if __name__=='__main__':run()
