"""Standalone rational certificate checker. Does not import the producer."""
import json,copy,math,time
from pathlib import Path
from fractions import Fraction as F
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
def mat(v):return np.array([[F(x) for x in row] for row in v],dtype=object)
def eye(n):return np.eye(n,dtype=object)*F(1)
def zero(n):return np.full((n,n),F(0),dtype=object)
def sympos(M,strict=False):
    assert np.array_equal(M,M.T)
    C=M.copy()
    while len(C):
        p=C[0,0];assert p>=0 and (not strict or p>0)
        if p==0:assert all(x==0 for x in C[0]);C=C[1:,1:]
        else:C=C[1:,1:]-C[1:,:1]@C[:1,1:]/p
def det(M):
    A=M.copy();out=F(1)
    for i in range(len(A)):
        p=next((j for j in range(i,len(A)) if A[j,i]),None)
        if p is None:return F(0)
        if p!=i:A[[i,p]]=A[[p,i]];out=-out
        z=A[i,i];out*=z
        for j in range(i+1,len(A)):A[j,i:]-=A[j,i]/z*A[i,i:]
    return out
def block(A,B):
    C=zero(len(A)+len(B));C[:len(A),:len(A)]=A;C[len(A):,len(A):]=B;return C
def stein(A,P,forcing,transpose=False):
    S=sum((a.T@P@a if transpose else a@P@a.T for a in A),zero(len(P)))
    assert np.array_equal(P-S,forcing)
def checkgrams(A,alpha,beta,P,Q,square):
    sympos(P);sympos(Q)
    stein(A,P,beta@beta.T);stein(A,Q,alpha.T@alpha,True)
    assert sum((P@Q)[i,i] for i in range(len(P)))==square
def check(m):
    As=[mat(a) for a in m['source_hidden_A']];alpha=mat(m['source_hidden_alpha']);beta=mat(m['source_hidden_beta'])
    dense=[mat(a) for a in m['source_A']];ad=mat(m['source_alpha']);bd=mat(m['source_beta']);S=mat(m['source_similarity'])
    assert np.array_equal(ad,alpha@S) and np.array_equal(S@bd,beta)
    for a,d in zip(As,dense):assert np.array_equal(S@d,a@S)
    P,Q=mat(m['source_P']),mat(m['source_Q']);sq=F(m['source_HS_squared'])
    checkgrams(As,alpha,beta,P,Q,sq);assert sq<=1
    X=mat(m['source_stability_X']);sympos(X,True);stein(As,X,eye(4))
    rows=[];cols=[]
    for w in m['rank_certificate']['row_words']:
        v=ad.copy()
        for a in w:v=v@dense[a]
        rows.append(v[0])
    for w in m['rank_certificate']['column_words']:
        V=eye(4)
        for a in w:V=V@dense[a]
        cols.append((V@bd)[:,0])
    assert det(np.array(rows,dtype=object))==F(m['rank_certificate']['row_determinant'])!=0
    assert det(np.array(cols,dtype=object).T)==F(m['rank_certificate']['column_determinant'])!=0
    R=[mat(a) for a in m['residual_A']];ar=mat(m['residual_alpha']);br=mat(m['residual_beta'])
    A=[mat(a) for a in m['reduced_A']];al=mat(m['reduced_alpha']);be=mat(m['reduced_beta'])
    assert np.array_equal(np.concatenate([al,ar],axis=1),alpha)
    assert np.array_equal(np.concatenate([be,br],axis=0),beta)
    for a,r,s in zip(A,R,As):assert np.array_equal(block(a,r),s)
    sq=F(m['model_error_HS_squared']);bound=F(m['model_error_bound'])
    checkgrams(R,ar,br,mat(m['residual_P']),mat(m['residual_Q']),sq)
    assert bound>=0 and bound**2>=sq and bound<=F(5,1000)
    checked=0
    for label,f in m['frames'].items():
        if not f['accepted']:continue
        aa=[mat(a) for a in f['A']];ab=mat(f['alpha']);bb=mat(f['beta'])
        realA=[mat(a) for a in f['unrounded_A']];realalpha=mat(f['unrounded_alpha']);realbeta=mat(f['unrounded_beta'])
        S=mat(f['S']);assert det(S)!=0
        assert np.array_equal(realalpha,al@S) and np.array_equal(S@realbeta,be)
        for a,t in zip(A,realA):assert np.array_equal(S@t,a@S)
        for x in [*aa,ab,bb]:assert all((v*2**24).denominator==1 for v in x.flat)
        J=[block(a,b) for a,b in zip(realA,aa)]
        ja=np.concatenate([realalpha,-ab],axis=1);jb=np.concatenate([realbeta,bb],axis=0)
        sqc=F(f['coefficient_error_HS_squared']);bc=F(f['coefficient_error_bound'])
        checkgrams(J,ja,jb,mat(f['coefficient_error_P']),mat(f['coefficient_error_Q']),sqc)
        assert bc>=0 and bc**2>=sqc
        X=mat(f['stability_X']);sympos(X,True);stein(aa,X,eye(2))
        M=mat(f['M']);inv=mat(f['inverse_comparison'])
        assert all(x>=0 for x in inv.flat)
        assert np.array_equal((eye(2)-M.T)@inv,eye(2))
        for a in aa:assert all(abs(a[i,j])<=M[i,j] for i in range(2) for j in range(2))
        # Positive resolvent applied to ones is a strict supersolution.
        v=inv@np.full((2,1),F(1));assert all(x>0 for x in v.flat)
        assert np.array_equal(v-M.T@v,np.full((2,1),F(1)))
        for pol in ['uniform','mixed']:
            c=f['best_'+pol]
            if c is None:continue
            steps=np.array([[F(x)] for x in c['steps']],dtype=object);eta=steps/2
            init=np.array(c['initial_integers'],dtype=object)[:,None]*steps
            E=mat(c['error_E']);Rr=mat(c['invariant_R'])
            assert np.array_equal(E-M.T@E,eta)
            assert all(abs(init[i,0]-ab[0,i])<=eta[i,0] for i in range(2))
            assert np.array_equal(Rr-M.T@Rr,abs(init)+eta)
            assert all(Rr[i,0] <= (2**(c['bits'][i]-1)-1)*steps[i,0] for i in range(2))
            err=(abs(bb).T@E)[0,0]
            assert err==F(c['execution_bound'])
            assert bound+bc+err==F(c['total_bound'])<=F(1,100)
            assert bound<=F(5,1000) and bc<=F(2,1000) and err<=F(3,1000)
            # Independently bound signed integer products using full bit ranges.
            shifts=c['shifts'];base=min(shifts)-24
            maxbound=0
            for a in aa:
              for j in range(2):
                z=sum(abs(int(a[i,j]*2**24))*2**(c['bits'][i]-1)*2**(shifts[i]-24-base) for i in range(2))
                maxbound=max(maxbound,z)
            assert maxbound<2**126
            assert int(maxbound).bit_length()+1==c['accumulator_signed_bits_bound']
            assert c['state_bytes']==math.ceil(sum(c['bits'])/8)
            checked+=1
    return {'id':m['id'],'source_rank_exactly_four':True,'source_Hankel_norm_at_most_one':True,
            'exact_gramians':True,'all_horizons_execution_certificate':True,'checked_allocations':checked,
            'additive_optimality_0_001_from_trivial_lower_bound':bound<=F(1,1000)}
def main():
    start=time.time();out=[]
    for p in sorted((ROOT/'models').glob('*.json')):
        m=json.loads(p.read_text())
        if 'source_hidden_A' not in m:continue
        r=check(m);out.append(r);print(json.dumps(r),flush=True)
    assert len(out)==36
    first=json.loads(next((ROOT/'models').glob('*.json')).read_text())
    rejects=[]
    for target in ['bound','matrix','invariant']:
        bad=copy.deepcopy(first)
        if target=='bound':bad['model_error_bound']='0'
        elif target=='matrix':bad['source_A'][0][0][0]='999'
        else:
            f=next(f for f in bad['frames'].values() if f['accepted']);f['best_mixed']['error_E'][0][0]='0'
        try:check(bad)
        except (AssertionError,ValueError):rejects.append(target)
    assert len(rejects)==3
    report={'all_passed':True,'instances':len(out),'records':out,'negative_controls_rejected':rejects,'elapsed_seconds':time.time()-start,
            'scope':'ordinary exact rational proof checker; compiler not formally verified'}
    (ROOT/'results/exact_verification.json').write_text(json.dumps(report,indent=2));print('PASS',len(out),flush=True)
if __name__=='__main__':main()
