#!/usr/bin/env python3
"""Independent rational replay, byte codec, and forward (all-controls) bounds."""
from fractions import Fraction as F
import json
from pathlib import Path
import numpy as np
import experiment as ex


def matvec(A,x):return [sum(a*b for a,b in zip(row,x)) for row in A]
def add(x,y):return [a+b for a,b in zip(x,y)]
def sub(x,y):return [a-b for a,b in zip(x,y)]
def absvec(x):return list(map(abs,x))
def fractions(A):return [[F(str(v)) for v in row] for row in A.tolist()]
def col(A):return [F(str(v)) for v in A]


def encode_decode(pre,bits,R):
    codes=[];word=0;shift=0
    for i,b in enumerate(bits):
        mid=2**(b-1)-1
        code=round(pre[i]/(R[i]/mid))+mid
        assert 0<=code<=2*mid
        word |= code<<shift;shift+=b;codes.append(code)
    payload=word.to_bytes(sum(bits)//8,'little')
    assert len(payload)==sum(bits)//8
    word=int.from_bytes(payload,'little');q=[];shift=0
    for i,b in enumerate(bits):
        mid=2**(b-1)-1;code=(word>>shift)&(2**b-1)
        assert code==codes[i]
        q.append((code-mid)*R[i]/mid);shift+=b
    return q


def main():
    out=ex.OUT
    rng=np.random.default_rng(917)
    tokens=rng.integers(0,2,size=96)
    nums=rng.integers(-1024,1025,size=96)
    traces=0;pointwise=0;future=0;floats_match=0
    sample=[]
    for ep in ('0.000001','0.001','0.01'):
        setup=ex.rational_setup(ep)
        a0,a1,km,bv,rv,invm,t1=setup['rational']
        A=[fractions(a0),fractions(a1)];K=fractions(km);B=col(bv);R=col(rv);inv=fractions(invm)
        for triangular in (False,True):
            T=[A[0],fractions(t1) if triangular else A[1]]
            E=[[[A[a][i][j]-T[a][i][j] for j in range(4)] for i in range(4)] for a in range(2)]
            D=[[max(abs(E[a][i][j]) for a in range(2)) for j in range(4)] for i in range(4)]
            for bits in ((6,4,4,2),(8,6,6,4)):
                delta=[R[i]/(2**bits[i]-2) for i in range(4)]
                source=add(matvec(D,R),delta)
                global_v=[min(z,2*r) for z,r in zip(matvec(inv,source),R)]
                x=[F(0)]*4;q=x.copy();v=x.copy();history=[]
                for t,(a,num) in enumerate(zip(tokens,nums)):
                    a=int(a);u=F(int(num),1024)
                    pre=add(matvec(T[a],q),[z*u for z in B])
                    assert all(abs(z)<=r for z,r in zip(pre,R))
                    qnew=encode_decode(pre,bits,R)
                    eta=sub(qnew,pre)
                    assert all(abs(z)<=d for z,d in zip(eta,delta))
                    local=sub(matvec(E[a],q),eta)
                    v=add(matvec([[abs(z) for z in row] for row in A[a]],v),absvec(local))
                    x=add(matvec(A[a],x),[z*u for z in B])
                    error=absvec(sub(x,qnew));q=qnew
                    assert all(e<=bb and e<=gg for e,bb,gg in zip(error,v,global_v))
                    history.append((v.copy(),error.copy()))
                    pointwise+=4
                for t in (0,8,16,32,48):
                    vf=history[t][0].copy()
                    for h in range(1,33):
                        vf=add(matvec(K,vf),source)
                        if h in (1,8,32):
                            assert all(e<=bb for e,bb in zip(history[t+h][1],vf))
                            future+=4
                sample.append({'epsilon':ep,'transition':'triangular' if triangular else 'full','bits':bits,
                               'last_error':[str(v) for v in history[-1][1]],
                               'last_bound':[str(v) for v in history[-1][0]],
                               'global_bound':[str(v) for v in global_v]})
                traces+=1
    # Necessary-term negative controls, independently and exactly.
    setup=ex.rational_setup('0.000001');a0,a1,km,bv,rv,invm,t1=setup['rational']
    q=bv
    structural_error=(a1-t1)*q
    assert any(v!=0 for v in structural_error)  # zero bound if coupling omitted fails
    R=col(rv);pre=col(bv);quant=encode_decode(pre,(4,4,4,4),R)
    assert any(a!=b for a,b in zip(pre,quant))  # zero bound if quantization omitted fails
    # A one-step local rounding counterexample: a tiny discarded coupling
    # crosses a grid boundary. The admissible grid state is not claimed to be
    # reachable from zero under the experiment's common input stream.
    rounding=[]
    for ep in ('0.000001','0.000000001'):
        ss=ex.rational_setup(ep)
        aa0,aa1,kk,bb,rr,ii,tt=ss['rational']
        r=col(rr);b=col(bb);s=[x/7 for x in r]
        q=[s[0],F(0),F(0),F(0)]
        epsilon=F(ep)
        drive=(s[1]/2-epsilon*q[0]/2)/b[1]
        assert abs(drive)<=1
        full_pre=add(matvec(fractions(aa1),q),[v*drive for v in b])
        tri_pre=add(matvec(fractions(tt),q),[v*drive for v in b])
        full_q=encode_decode(full_pre,(4,4,4,4),r)
        tri_q=encode_decode(tri_pre,(4,4,4,4),r)
        assert full_q[1]-tri_q[1]==s[1]
        rounding.append({'epsilon':ep,'local_input':str(drive),
                         'pre_quantization_difference':str(full_pre[1]-tri_pre[1]),
                         'post_quantization_difference':str(full_q[1]-tri_q[1]),
                         'post_quantization_difference_float':float(full_q[1]-tri_q[1]),
                         'scope':'one-step admissible grid state; reachability from zero not claimed'})
    ex.save_json(out/'rounding_boundary_control.json',rounding)
    rows=json.loads((out/'results.json').read_text());raw=np.load(out/'per_sequence_evidence.npz')
    for row in rows:
        prefix=f"{row['epsilon']}_{row['budget_bits']}_{row['transition']}"
        bits=raw[prefix+'_bits'];g2=np.square(row['gains'])
        if row['policy']=='calibrated':
            chosen=int(np.argmin(np.mean(raw[prefix+'_calibration_mse']@g2/4,axis=1)))
            assert bits[chosen].tolist()==row['bits']
        if row['policy']=='global_certificate':
            chosen=int(np.argmin(raw[prefix+'_test_global_v']**2@g2/4))
            assert bits[chosen].tolist()==row['bits']
        idx=int(np.flatnonzero(np.all(bits==row['bits'],axis=1))[0])
        mse=float(np.mean(raw[prefix+'_test_mse'][idx]@g2/4))
        assert abs(mse-row['test_mse'])<1e-12
        assert sum(row['bits'])==8*row['persistent_state_bytes']
    ex.save_json(out/'exact_replays.json',sample)
    result={'status':'PASS','rational_replay_traces':traces,'steps_per_trace':96,
            'exact_pointwise_coordinate_checks':pointwise,
            'exact_future_coordinate_checks_horizons_1_8_32':future,
            'negative_control_omitting_coupling':'detected failure',
            'negative_control_omitting_quantization':'detected failure',
            'local_quantization_boundary_control':'one-bin jump at epsilon 1e-6 and 1e-9, exact rational calculation',
            'published_rows_recomputed':len(rows),'calibration_and_analytic_selection_verified':True,
            'byte_budget_verified':True,
            'scope':'Exact tests concern the rational recurrence and rational nearest-grid codec, not all floating-point executions.'}
    ex.save_json(out/'independent_checks.json',result)
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
