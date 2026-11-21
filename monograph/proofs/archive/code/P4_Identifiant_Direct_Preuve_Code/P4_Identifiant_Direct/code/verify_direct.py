#!/usr/bin/env python3
"""Exact exhaustive certificate; enumeration occurs only in this verifier."""
from pathlib import Path
from fractions import Fraction as F
import ctypes as C
import json, subprocess, sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'reference'))
from oracle import get_nodes, bounds, value, bayes

def encode(p,word):
    # Independent string rule, not the integer formula used by C.
    if p in (19,21):
        B=(p+1)//2
        if len(word)==B:
            assert all(word[i]!=word[i+1] for i in range(B-2))
            return (1<<(B-1))+2*int(word[0])+int(word[-1])
    return int('1'+word,2)

def state(index,nb): return (C.c_uint8*nb).from_buffer_copy(index.to_bytes(nb,'little'))
def fraction(out): return F(out[0]+(out[1]<<64),out[2]+(out[3]<<64))

def main():
    subprocess.run(['gcc','-std=c11','-O2','-Wall','-Wextra','-Werror',
                    '-Wno-misleading-indentation','-shared','-fPIC',
                    str(ROOT/'code/direct_filter.c'),'-o',str(ROOT/'code/libdirect.so')],check=True)
    lib=C.CDLL(str(ROOT/'code/libdirect.so'))
    U8=C.POINTER(C.c_uint8);U64=C.POINTER(C.c_uint64)
    lib.df_step.argtypes=[C.c_uint,U8,C.c_uint,C.c_uint]
    lib.df_output.argtypes=[C.c_uint,U8,C.c_uint,U64]
    lib.df_step_output.argtypes=[C.c_uint,U8,C.c_uint,C.c_uint,U64]
    lib.df_advance.argtypes=[C.c_uint,U8,C.c_uint,C.c_uint,U64]
    rows=[];transitions=outputs=invalid=0
    for p in range(4,23):
        nodes=get_nodes(p); byword={r['word']:r for r in nodes}
        n=len(nodes);B=(n-1).bit_length();nb=(B+7)//8
        codes={r['word']:encode(p,r['word']) for r in nodes}
        assert len(set(codes.values()))==n
        assert max(codes.values())<1<<B
        expected={}
        for r in nodes:
            center=(r['lo']+r['hi'])/2 if r['leaf'] else value(r['word'],F(1,2))
            expected[codes[r['word']]]=F(3,4)-center/2
        out=(C.c_uint64*4)()
        for index in range(1<<B):
            st=state(index,nb)
            rc=lib.df_output(p,st,nb,out)
            if index in expected:
                assert rc==0 and fraction(out)==expected[index];outputs+=1
            else:
                assert rc==-2;invalid+=1
        eps=F(49,132*(1<<p))
        for r in nodes:
            w=r['word'];index=codes[w]
            for b in (0,1):
                raw=str(b)+w;new=raw
                for k in range(1,len(raw)+1):
                    lo,hi=bounds(raw[:k])
                    if hi-lo<=4*eps: new=raw[:k];break
                assert new in byword
                nr=byword[new]
                if r['leaf']:
                    assert nr['leaf']
                    assert nr['lo']<=bayes(r['lo'],b)<=bayes(r['hi'],b)<=nr['hi']
                else:
                    x=bayes(value(w,F(1,2)),b)
                    assert nr['lo']<=x<=nr['hi']
                    if not nr['leaf']: assert x==value(new,F(1,2))
                for mode in ('separate','fused','advance'):
                    st=state(index,nb)
                    if mode=='fused':rc=lib.df_step_output(p,st,nb,b,out)
                    elif mode=='advance':rc=lib.df_advance(p,st,nb,b,out)
                    else:
                        rc=lib.df_step(p,st,nb,b)
                        assert lib.df_output(p,st,nb,out)==0
                    assert rc==0
                    assert int.from_bytes(bytes(st),'little')==codes[new]
                    assert fraction(out)==expected[codes[new]]
                    transitions+=1;outputs+=1
        previous=json.loads((ROOT/'reference'/f'certificate_p{p:02d}.json').read_text())
        xs=[value(w[::-1],F(1,2)) for w in previous['witness_words']]
        assert all((b-a)/2>2*eps for a,b in zip(xs,xs[1:]))
        LB=(len(xs)-1).bit_length()
        excess=[r for r in nodes if int('1'+r['word'],2)>=1<<B]
        special=[]
        if excess:
            assert p in (19,21) and len(excess)==4
            H=1<<(B-1)
            for h in range(H,H+4):assert bin(h)[3:] not in byword
            for r in excess:
                special.append({'word':r['word'],'heap_code':int('1'+r['word'],2),'direct_id':codes[r['word']]})
        row={'p':p,'epsilon':str(eps),'states':n,'bits':B,'lower_witnesses':len(xs),
             'lower_bits':LB,'optimal_bits_certified':LB==B,'depth':max(len(r['word']) for r in nodes),
             'state_bytes':nb,'remapped':special}
        rows.append(row);print(json.dumps(row),flush=True)
        (ROOT/'results'/f'direct_p{p:02d}.bin').write_bytes(bytes([p,nb])+int(1).to_bytes(nb,'little'))
    # Finite-width bounds from the existing mathematical proof (p<=30).
    q=F(16,49)
    assert q**20/F(4)<=F(49,33*(1<<30))
    assert 7**20<2**64
    assert max(49*(16*7**20)**2,528*9**20*(1<<30))<2**128
    # Falsify extension of the four-exception rule to the next odd precision.
    ns=get_nodes(23);B=(len(ns)-1).bit_length()
    over=[r['word'] for r in ns if int('1'+r['word'],2)>=1<<B]
    assert len(over)==44
    result={'passed':True,'configs':rows,'transition_checks':transitions,'exact_output_checks':outputs,
            'unused_id_rejections':invalid,'optimal_p':[r['p'] for r in rows if r['optimal_bits_certified']],
            'runtime_transition_tables':0,'runtime_exception_tables':0,
            'immutable_parameter_bytes':2,'finite_width_proof_passed':True,
            'scope_counterexample':{'p':23,'states':len(ns),'bits':B,'excess_paths':len(over)},
            'complexity':'O(depth) integer operations, O(depth) bits of transient state; no state-tree enumeration',
            'qualification':'C API p=4..22; arbitrary precision optimal direct coding not proved'}
    (ROOT/'results/VERIFIED.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASSED',transitions,outputs,invalid,flush=True)

if __name__=='__main__':main()
