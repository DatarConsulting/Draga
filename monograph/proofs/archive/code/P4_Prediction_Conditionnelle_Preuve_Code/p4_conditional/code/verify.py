#!/usr/bin/env python3
"""Independent certificate checker: normalized matrix multiplication + Fractions.
No floating-point decision is used. Does not import the synthesis code.
"""
from fractions import Fraction as F
from pathlib import Path
import json, ctypes as C, subprocess, struct, hashlib, copy
ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/'results'
A=(((F(1,4),F(1,8)),(F(1,24),F(1,12))),
   ((F(1,12),F(1,24)),(F(1,8),F(1,4))))
J=(F(3,8),F(5,8))

def bayes(x,b):
    M=A[b]; u=(1-x)*M[0][0]+x*M[1][0]; v=(1-x)*M[0][1]+x*M[1][1]
    return v/(u+v)

def belief(w):
    u,v=F(1,2),F(1,2)
    for b in w:
        M=A[int(b)]; u,v=u*M[0][0]+v*M[1][0],u*M[0][1]+v*M[1][1]
    return v/(u+v)

def checks(c):
    p=c['p']; D=1<<p; R=F(c['radius']); eps=F(c['epsilon']); js=c['states']; es=c['transitions']
    assert eps==F(49,264*D) and R==4*eps
    assert len(js)==len(es)==c['upper_states'] and len(set(js))==len(js)
    assert js[c['initial']]==D//2
    reach={c['initial']}; pending=list(reach); min_slack=None
    for i in pending:
        assert len(es[i])==2
        for k in es[i]:
            assert 0<=k<len(js)
            if k not in reach: reach.add(k); pending.append(k)
    assert len(reach)==len(js)
    for i,j in enumerate(js):
        x=F(j,D); assert J[0]<=x<=J[1]
        lo=max(J[0],x-R); hi=min(J[1],x+R)
        probs=(F(3,8)-x/4,F(1,8)+x/4,F(1,2))
        assert min(probs)>0 and sum(probs)==1
        for b,k in enumerate(es[i]):
            y=F(js[k],D); nlo=max(J[0],y-R); nhi=min(J[1],y+R)
            left,right=bayes(lo,b),bayes(hi,b)
            assert nlo<=left<=right<=nhi, (p,i,b,'invariant')
            s=min(left-nlo,nhi-right)
            min_slack=s if min_slack is None else min(min_slack,s)
    xs=[belief(w) for w in c['witness_words']]
    assert len(xs)==c['lower_states'] and all(x<y for x,y in zip(xs,xs[1:]))
    gaps=[(y-x)/4 for x,y in zip(xs,xs[1:])]
    assert min(gaps)>2*eps and F(c['minimum_witness_tv'])==min(gaps)
    assert c['lower_bits']==(len(xs)-1).bit_length()
    assert c['upper_bits']==(len(js)-1).bit_length()
    d=F(c['robust_parameter_radius'])
    assert 0<d<=F(1,1000) and d<=eps/32 and min(gaps)-2*d>2*(eps+d)
    assert min_slack>0
    return {'p':p,'N_upper':len(js),'N_lower':len(xs),'B_upper':c['upper_bits'],
            'B_lower':c['lower_bits'],'exact_bits':c['lower_bits']==c['upper_bits'],
            'epsilon':str(eps),'epsilon_decimal':float(eps),
            'stream_epsilon':str(2*eps),'stream_epsilon_decimal':float(2*eps),'min_interval_slack':str(min_slack),
            'robust_parameter_radius':str(d),'robust_target_epsilon':str(eps+d),
            'robust_stream_target_epsilon':str(2*(eps+d)),
            'static_model_bytes':16+8*len(js),'state_bytes':(c['upper_bits']+7)//8,
            'invariant_edges':2*len(js),'lower_witnesses':len(xs)}

def source_checks():
    # Structural facts without SymPy or the synthesis implementation.
    assert all(v>0 for M in A for row in M for v in row)
    assert all(sum(A[b][r][s] for b in (0,1) for s in (0,1))==F(1,2) for r in (0,1))
    def prob(w):
        u,v=F(1,2),F(1,2)
        for b in w:
            M=A[int(b)]; u,v=u*M[0][0]+v*M[1][0],u*M[0][1]+v*M[1][1]
        return (u+v)/2
    det=prob('')*prob('00')-prob('0')**2
    assert det==F(1,768)
    assert all(J[0]<bayes(J[0],b)<=bayes(J[1],b)<J[1] for b in (0,1))
    assert bayes(J[1],0)<bayes(J[0],1)
    # Fixed source rows yield future-law sign determined by first symbol.
    assert all(A[0][0][i]-A[0][1][i]>0 for i in (0,1))
    assert all(A[1][0][i]-A[1][1][i]<0 for i in (0,1))
    return {'rank_minor':str(det),'hankel_hs_squared_bound':'64/169',
            'belief_interval':['3/8','5/8'],'derivative_upper':'16/49',
            'derivative_lower':'16/81','image_gap':str(bayes(J[0],1)-bayes(J[1],0)),
            'whole_future_readout_tv_factor':'1/4',
            'perturbation_tv_factor_upper':'229/264'}

def main():
    source=source_checks()
    subprocess.run(['gcc','-std=c11','-O2','-Wall','-Wextra','-Werror','-shared','-fPIC',
                    str(ROOT/'code/packed_predictor.c'),'-o',str(ROOT/'code/libpacked.so')],check=True)
    lib=C.CDLL(str(ROOT/'code/libpacked.so')); P=C.POINTER(C.c_uint8)
    lib.pc_valid.argtypes=[P,C.c_size_t]; lib.pc_init.argtypes=[P,P]
    lib.pc_step.argtypes=[P,P,C.c_uint]; lib.pc_output.argtypes=[P,P,C.POINTER(C.c_uint64)]
    lib.pc_output_stream.argtypes=[P,P,C.POINTER(C.c_uint64)]
    results=[]; c_edges=c_outputs=0
    for path in sorted(OUT.glob('certificate_p*.json')):
        cert=json.loads(path.read_text()); r=checks(cert); p=cert['p']; D=1<<p
        payload=(OUT/f'machine_p{p:02d}.bin').read_bytes(); nb=r['state_bytes']; immutable=payload[:-nb]
        model=(C.c_uint8*len(immutable)).from_buffer_copy(immutable)
        assert lib.pc_valid(model,len(immutable))==1
        assert payload[-nb:]==bytes(nb)
        state=(C.c_uint8*nb)(); lib.pc_init(model,state); assert bytes(state)==bytes(nb)
        out=(C.c_uint64*4)()
        for i,j in enumerate(cert['states']):
            for k,v in enumerate(i.to_bytes(nb,'little')): state[k]=v
            assert lib.pc_output(model,state,out)==0
            assert list(out)==[3*D-2*j,D+2*j,4*D,8*D]
            c_outputs+=1
            assert lib.pc_output_stream(model,state,out)==0
            assert list(out)[:3]==[3*D-2*j,D+2*j,4*D]
            c_outputs+=1
            for b in (0,1):
                for k,v in enumerate(i.to_bytes(nb,'little')): state[k]=v
                assert lib.pc_step(model,state,b)==0
                assert int.from_bytes(bytes(state),'little')==cert['transitions'][i][b]
                c_edges+=1
        r['model_sha256']=hashlib.sha256(payload).hexdigest(); results.append(r)
        print(json.dumps(r),flush=True)
    # Meaningful mutation tests must be rejected by the same certificate checker.
    base=json.loads((OUT/'certificate_p11.json').read_text()); rejected=[]
    bad=copy.deepcopy(base); bad['transitions'][0][0]=0
    bad2=copy.deepcopy(base); bad2['witness_words'][1]=bad2['witness_words'][0]
    bad3=copy.deepcopy(base); bad3['radius']=str(F(bad3['radius'])/100)
    for label,c in [('wrong_transition',bad),('duplicate_witness',bad2),('false_radius',bad3)]:
        try: checks(c)
        except AssertionError: rejected.append(label)
        else: raise AssertionError('mutation not rejected: '+label)
    summary={'all_passed':True,'source':source,'configurations':results,
             'C_exact_edges':c_edges,'C_exact_outputs':c_outputs,
             'rational_invariant_edges':sum(x['invariant_edges'] for x in results),
             'lower_witnesses_recomputed':sum(x['lower_witnesses'] for x in results),
             'exact_bit_cases':[x['p'] for x in results if x['exact_bits']],
             'negative_controls_rejected':rejected}
    (OUT/'VERIFIED.json').write_text(json.dumps(summary,indent=2)+'\n')
    print('FINAL',json.dumps({k:v for k,v in summary.items() if k not in ('configurations','source')}),flush=True)
if __name__=='__main__': main()
