"""Exact finite-state separation lower bounds and k=3,4 examples."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import json,itertools,math,ctypes
from pathlib import Path
from fractions import Fraction as F
import numpy as np
import sympy as sp
import experiment as ex
import verify as ve
from replay import pack,unpack,py_step
ROOT=Path(__file__).resolve().parents[1]
def lower_bounds():
    out=[]
    suffixes=[()]+[w for d in range(1,5) for w in itertools.product(range(2),repeat=d)]
    for path in sorted((ROOT/'models').glob('*.json')):
        m=json.loads(path.read_text());A=[ve.mat(a) for a in m['source_A']];alpha=ve.mat(m['source_alpha']);beta=ve.mat(m['source_beta'])
        Af=[np.array(a,float) for a in A];cols=[]
        for w in suffixes:
            v=np.eye(4)
            for a in w:v=v@Af[a]
            cols.append((v@np.array(beta,float))[:,0])
        cols=np.array(cols).T;words=[()];states=[np.array(alpha,float)[0]];lev=[0]
        for d in range(10):
            nxt=[]
            for i in lev:
                for a in range(2):nxt.append(len(words));words.append(words[i]+(a,));states.append(states[i]@Af[a])
            lev=nxt
        features=np.array(states)@cols;chosen=[]
        for i in range(len(words)):
            if not chosen or np.all(np.max(abs(features[chosen]-features[i]),axis=1)>.02000001):chosen.append(i)
        witnesses=[]
        for ii,i in enumerate(chosen):
          for j in chosen[:ii]:
            t=int(np.argmax(abs(features[i]-features[j])));u=words[i];v=words[j];w=suffixes[t]
            def exact(word):
                s=alpha.copy()
                for a in word:s=s@A[a]
                return (s@beta)[0,0]
            gap=abs(exact(u+w)-exact(v+w));assert gap>F(1,50)
            witnesses.append({'u':u,'v':v,'suffix':w,'gap':str(gap)})
        out.append({'id':m['id'],'separated_prefixes':len(chosen),'lower_bound_bits':math.ceil(math.log2(len(chosen))),
          'constructive_bits':m['best']['mixed']['state_bits'],'witnesses':witnesses})
    (ROOT/'results/lower_bounds.json').write_text(json.dumps({'epsilon':'1/100','records':out,
      'scope':'Any deterministic time-homogeneous machine with bounded total persistent state and pointwise error <=epsilon on every word. Search finite, witnesses exact; no tightness claim.'},indent=2))
def highrank():
    records=[]
    lib=ctypes.CDLL(str(ROOT/'code/kernel.so'));fn=lib.step
    fn.argtypes=[ctypes.POINTER(ctypes.c_uint8),ctypes.c_int,ctypes.c_int,ctypes.POINTER(ctypes.c_int32),ctypes.POINTER(ctypes.c_int32),ctypes.POINTER(ctypes.c_int32),ctypes.c_int,ctypes.c_int]
    for k in [3,4]:
      for scale in [sp.Rational(3,4),sp.Rational(1)]:
        A0=sp.diag(*[sp.Rational(42+5*i,128) for i in range(k)])
        A1=sp.diag(*[sp.Rational(65-4*i,128) for i in range(k)])
        for i in range(k):A0[i,(i+1)%k]=sp.Rational(1,64);A1[(i+1)%k,i]=sp.Rational(1,128)
        A=[A0*scale,A1*scale];al=sp.ones(1,k)/4;be=sp.Matrix([sp.Rational(1,2**(2*i+2)) for i in range(k)])
        R=[sp.diag(sp.Rational(1,8),sp.Rational(1,16)),sp.Matrix([[sp.Rational(1,16),sp.Rational(1,64)],[sp.Rational(1,128),sp.Rational(1,8)]])]
        ar=sp.ones(1,2)/4;br=sp.ones(2,1)/2**16
        As=[sp.diag(a,r) for a,r in zip(A,R)];af=al.row_join(ar);bf=be.col_join(br)
        sq,P,O=ex.hsq(As,af,bf);assert sq<=1
        model_sq,Pr,Or=ex.hsq(R,ar,br);model_bound=ex.sqrt_up(model_sq)
        M=sp.Matrix(k,k,lambda i,j:max(abs(a[i,j]) for a in A));inv=(sp.eye(k)-M.T).inv()
        R0=inv*al.T;exp=[math.ceil(math.log2(2*float(v)+1e-5)) for v in R0]
        frame={'A':A,'alpha':al,'beta':be,'M':M,'inv':inv,'exp':exp}
        uniform=next(c for p in range(3,33) if (c:=ex.allocation(frame,[p]*k)) and F(c['execution_bound'])<=F(3,1000))
        weights=be.T*inv
        # Dynamic program for a separable error surrogate, then exact feasibility check.
        dp={0:(F(0),[])}
        for i in range(k):
            nd={}
            for used,(cost,bits) in dp.items():
              for p in range(3,33):
                b=used+p
                if b>64:continue
                v=cost+ex.frac(weights[i])*F(2)**(exp[i]-p)
                if b not in nd or v<nd[b][0]:nd[b]=(v,bits+[p])
            dp=nd
        mixed=None
        for B in sorted(dp):
            c=ex.allocation(frame,dp[B][1])
            if c and F(c['execution_bound'])<=F(3,1000):mixed=c;break
        assert mixed
        # Independent exact residual checks (checker uses Fraction, not producer).
        aa=[ve.mat(ex.ser(a)) for a in A];mm=ve.mat(ex.ser(M));iv=ve.mat(ex.ser(inv))
        assert np.array_equal((ve.eye(k)-mm.T)@iv,ve.eye(k));assert all(x>=0 for x in iv.flat)
        ve.checkgrams([ve.mat(ex.ser(a)) for a in As],ve.mat(ex.ser(af)),ve.mat(ex.ser(bf)),ve.mat(ex.ser(P)),ve.mat(ex.ser(O)),F(str(sq)))
        X=ex.gram(As,sp.eye(k+2));ve.sympos(ve.mat(ex.ser(X)),True);ve.stein([ve.mat(ex.ser(a)) for a in As],ve.mat(ex.ser(X)),ve.eye(k+2))
        rank=ex.rank_certificate(As,af,bf,k+2)
        trials=[]
        for pol,c in [('uniform',uniform),('mixed',mixed)]:
            eta=np.array([[F(v)/2] for v in c['steps']],object);Er=ve.mat(c['error_E']);Rr=ve.mat(c['invariant_R'])
            ini=np.array(c['initial_integers'],object)[:,None]*eta*2
            assert np.array_equal(Er-mm.T@Er,eta) and np.array_equal(Rr-mm.T@Rr,abs(ini)+eta)
            assert all(Rr[i,0]<=(2**(c['bits'][i]-1)-1)*eta[i,0]*2 for i in range(k))
            bits=np.array(c['bits'],np.int32);sh=np.array(c['shifts'],np.int32)
            ai=np.array([[[int(v*2**24) for v in row] for row in np.array(a,object)] for a in A],np.int32)
            buf=bytearray(pack(c['initial_integers'],bits));rng=np.random.default_rng(k);steps=100000;source=np.array(af,float)[0];maxerr=0
            for a in rng.integers(0,2,steps):
                prev=unpack(buf,bits);expected=py_step(prev,ai,bits.tolist(),sh.tolist(),int(a),24)
                ret=fn((ctypes.c_uint8*len(buf)).from_buffer(buf),len(buf),k,ai.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),bits.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),sh.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),int(a),24)
                assert ret==0 and unpack(buf,bits)==expected
                source=source@np.array(As[a],float);y=float(np.array(expected)*2.**sh@np.array(be,float)[:,0]);maxerr=max(maxerr,abs(float(source@np.array(bf,float)[:,0])-y))
            trials.append({'policy':pol,'allocation':c,'C_python_steps':steps,'observed_error':maxerr,'total_bound':str(model_bound+F(c['execution_bound']))})
        records.append({'k':k,'n':k+2,'scale':str(scale),'source_rank_certificate':rank,'source_HS_squared':str(sq),
          'source_A':[ex.ser(a) for a in As],'source_alpha':ex.ser(af),'source_beta':ex.ser(bf),
          'source_P':ex.ser(P),'source_Q':ex.ser(O),'source_X':ex.ser(X),
          'M':ex.ser(M),'inverse':ex.ser(inv),'model_bound':str(model_bound),'trials':trials})
        print('higher-rank',k,scale,uniform['state_bits'],mixed['state_bits'],flush=True)
    (ROOT/'results/higher_rank.json').write_text(json.dumps({'all_passed':True,'records':records,
      'scope':'Dyadic positive examples; known weak block removed; general entrywise invariant proof, no universal compression or new rank theorem.'},indent=2))
if __name__=='__main__':lower_bounds();highrank()
