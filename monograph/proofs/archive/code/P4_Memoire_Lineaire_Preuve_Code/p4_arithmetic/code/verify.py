"""Exhaustive comparison with independent rational belief-coordinate oracle."""
import sys,ctypes as C,json,time
from fractions import Fraction as F
from pathlib import Path
from fast import encode,bits,append,branch,leaf
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'reference'))
from oracle import nodes,bayes
lib=C.CDLL(str(ROOT/'code/libarithmetic.so'))
U=C.c_uint;U64=C.c_uint64;PTR=C.POINTER
lib.ac_branch.argtypes=[U64,U64];lib.ac_branch.restype=U
lib.ac_encode.argtypes=[U,U64,PTR(U64)]
lib.ac_decode.argtypes=[U,U64,PTR(U64)]
lib.ac_output.argtypes=[U,PTR(C.c_uint8),U,PTR(U64)]
lib.ac_advance.argtypes=[U,PTR(C.c_uint8),U,U,PTR(U64)]
def packed(x,n):return (C.c_uint8*n).from_buffer_copy(x.to_bytes(n,'little'))
def fraction(o):return F(o[0]+(o[1]<<64),o[2]+(o[3]<<64))
def main():
    transitions=outputs=roundtrips=invalid=0;rows=[]
    for p in [19,21,23,25,27,29,30]:
        start=time.perf_counter();ns=nodes(p);byword={z['word']:z for z in ns};B=bits(p);nb=(B+7)//8
        ids={};min_width=1<<B;out=(U64*4)();rid=U64();path=U64()
        for z in ns:
            w=z['word'];i,N=encode(p,w);min_width=min(min_width,N);ids[w]=i
            assert lib.ac_encode(p,int('1'+w,2),C.byref(rid))==0 and rid.value==i
            assert lib.ac_decode(p,i,C.byref(path))==0 and path.value==int('1'+w,2)
            cur=(1,0,0,1)
            for b in map(int,w):cur=append(cur,b)
            S,T=cur[0]+cur[2],cur[1]+cur[3]
            assert lib.ac_branch(S,T)==branch(S,T)
            assert leaf(cur,p)==z['leaf'];roundtrips+=1
        assert len(set(ids.values()))==len(ns) and min_width>=2
        with (ROOT/f'results/graph_{p}.txt').open('w') as stream:
            stream.write(f'{p} {B} {nb} {len(ns)}\n')
            for z in ns:
                w=z['word'];i=ids[w]
                assert lib.ac_output(p,packed(i,nb),nb,out)==0 and fraction(out)==z['out'];outputs+=1
                expected=[]
                for bit in [0,1]:
                    raw=str(bit)+w;new=raw
                    for k in range(1,len(raw)+1):
                        if byword[raw[:k]]['leaf']:new=raw[:k];break
                    nz=byword[new];expected.append(ids[new])
                    if z['leaf']:
                        assert nz['leaf'] and nz['lo']<=bayes(z['lo'],bit)<=bayes(z['hi'],bit)<=nz['hi']
                        assert (nz['hi']-nz['lo'])/4<=F(49,132*(1<<p))
                    else:
                        x=bayes(z['center'],bit);assert nz['lo']<=x<=nz['hi']
                        if not nz['leaf']:assert x==nz['center']
                    st=packed(i,nb);assert lib.ac_advance(p,st,nb,bit,out)==0
                    assert int.from_bytes(bytes(st),'little')==ids[new] and fraction(out)==nz['out']
                    transitions+=1;outputs+=1
                stream.write(f'{i} {expected[0]} {expected[1]}\n')
        # Exhaustive rejection including all unused IDs, not just sampled faults.
        valid=set(ids.values());bad=0
        for i in range(1<<B):
            rc=lib.ac_decode(p,i,C.byref(path));assert (rc==0)==(i in valid)
            bad+=bool(rc)
        invalid+=bad
        row={'p':p,'actual_states':len(ns),'minimum_index_bits':(len(ns)-1).bit_length(),
             'arithmetic_bits':B,'state_bytes':nb,'minimum_interval_width':min_width,
             'unused_ids_rejected':bad,'max_depth':max(map(len,ids)),
             'seconds_local':time.perf_counter()-start}
        rows.append(row);print(json.dumps(row),flush=True)
    result={'passed':True,'graphs':rows,'roundtrips':roundtrips,'transitions':transitions,
            'exact_outputs':outputs,'unused_ids_rejected':invalid}
    (ROOT/'results/VERIFIED.json').write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':main()
