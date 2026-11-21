from pathlib import Path
from fractions import Fraction as F
from collections import Counter
import ctypes as C
import json,subprocess,time
from uniform import Codec,compile_lib,Context,count_rolling
ROOT=Path(__file__).resolve().parents[1]

# Independent projective belief coordinates, derived from original T_b.
T=[((F(1,2),F(1,4)),(F(1,12),F(1,6))),((F(1,6),F(1,12)),(F(1,4),F(1,2)))]
RB=[]
for ((a,b),(c,d)) in T:
    vals=[12*(d-b),12*b,12*(c+d-a-b),12*(a+b)]
    assert all(x.denominator==1 for x in vals)
    RB.append(tuple(int(x) for x in vals))
assert RB==[(-1,3,-6,9),(5,1,6,3)]
def bayes(x,b):
    ((a,c),(d,e))=T[b]
    left=(1-x)*a+x*d;right=(1-x)*c+x*e
    return right/(left+right)
def nodes(p):
    eps=F(49,132*(1<<p));out=[]
    def visit(w,a,b,c,d):
        lo=F(3*a+8*b,3*c+8*d);hi=F(5*a+8*b,5*c+8*d)
        leaf=hi-lo<=4*eps
        center=(lo+hi)/2 if leaf else F(a+2*b,c+2*d)
        out.append({'word':w,'lo':lo,'hi':hi,'leaf':leaf,'center':center,'out':F(3,4)-center/2})
        if not leaf:
            for bit,(A,B,C,D) in enumerate(RB):
                visit(w+str(bit),a*A+b*C,a*B+b*D,c*A+d*C,c*B+d*D)
    visit('',1,0,0,1);return out
def fold23(w):
    if len(w)!=12:return int('1'+w,2)
    places=[i for i in range(1,11) if w[i]==w[i-1]]
    assert len(places)<=1
    defect=places[0] if places else 0
    k=2*defect+int(w[-1])
    return 4095-k if w[0]=='1' else 2048+k
def fraction(o):return F(o[0]+(o[1]<<64),o[2]+(o[3]<<64))
def pack(i,nb):return (C.c_uint8*nb).from_buffer_copy(i.to_bytes(nb,'little'))
def main():
    lib=compile_lib()
    subprocess.run(['gcc','-std=c11','-O2','-shared','-fPIC',str(ROOT/'code/fold23.c'),'-o',str(ROOT/'code/libfold23.so')],check=True)
    fold=C.CDLL(str(ROOT/'code/libfold23.so'))
    fold.df_output.argtypes=[C.c_uint,C.POINTER(C.c_uint8),C.c_uint,C.POINTER(C.c_uint64)]
    fold.df_advance.argtypes=[C.c_uint,C.POINTER(C.c_uint8),C.c_uint,C.c_uint,C.POINTER(C.c_uint64)]
    # Independent Python DP cross-check and rolling-space counter.
    for p in [4,6,8]:
        a=Codec(p,'python');b=Codec(p,'c',lib)
        assert list(b.counts)==a.counts
        assert count_rolling(p,a.initial,a.R)==a.size
    configs=[];transitions=outputs=roundtrips=0;foldchecks=0;structural=[]
    for p in [19,21,23,25,27,29,30]:
        start=time.perf_counter();co=Codec(p,'c',lib);build=time.perf_counter()-start
        ns=nodes(p);byword={z['word']:z for z in ns}
        truebits=(len(ns)-1).bit_length()
        ex=[z['word'] for z in ns if len(z['word'])>=truebits]
        parents={w[:-1] for w in ex}
        defects=Counter(sum(a==b for a,b in zip(w,w[1:])) for w in parents)
        structural.append({'p':p,'states':len(ns),'construction_bits':truebits,'exception_paths':len(ex),'parent_defects':dict(sorted(defects.items()))})
        if p==23:
            assert len(ex)==44 and dict(defects)=={0:2,1:20}
            assert all(bin(q)[3:] not in byword for q in list(range(2048,2070))+list(range(4074,4096)))
            fids={w:fold23(w) for w in byword}
            assert len(set(fids.values()))==len(ns) and max(fids.values())<4096
        # Brute language enumeration only in this verifier, to check all ranks.
        allwords=[]
        def visit(w,j,r):
            allwords.append(w)
            for b in (0,1):
                k=co.cost[2*j+b]
                if r>=k:visit(w+str(b),co.next[2*j+b],r-k)
        visit('',co.initial,co.R)
        assert len(allwords)==co.size
        assert co.bits<=F(50501,100000)*(p+1)+8 # stricter finite check using lower bounds on d and L
        for i,w in enumerate(allwords):
            assert co.encode(w)==i and co.decode(i)==w
            outid=C.c_uint64();path=C.c_uint64()
            assert co.lib.uc_encode(C.byref(co.ctx),int('1'+w,2),C.byref(outid))==0 and outid.value==i
            assert co.lib.uc_decode(C.byref(co.ctx),i,C.byref(path))==0 and path.value==int('1'+w,2)
            roundtrips+=1
        ids={w:co.encode(w) for w in byword}
        assert len(set(ids.values()))==len(ns)
        out=(C.c_uint64*4)();eps=F(49,132*(1<<p))
        for z in ns:
            w=z['word'];st=pack(ids[w],co.bytes)
            assert co.lib.uc_output(C.byref(co.ctx),st,co.bytes,out)==0 and fraction(out)==z['out'];outputs+=1
            for bit in (0,1):
                raw=str(bit)+w;new=raw
                for k in range(1,len(raw)+1):
                    if byword[raw[:k]]['leaf']:new=raw[:k];break
                nz=byword[new]
                if z['leaf']:
                    assert nz['leaf'] and nz['lo']<=bayes(z['lo'],bit)<=bayes(z['hi'],bit)<=nz['hi']
                    assert (nz['hi']-nz['lo'])/4<=eps
                else:
                    x=bayes(z['center'],bit)
                    assert nz['lo']<=x<=nz['hi']
                    if not nz['leaf']:assert x==nz['center']
                st=pack(ids[w],co.bytes)
                assert co.lib.uc_advance(C.byref(co.ctx),st,co.bytes,bit,out)==0
                assert int.from_bytes(bytes(st),'little')==ids[new] and fraction(out)==nz['out']
                transitions+=1;outputs+=1
                if p==23:
                    sf=pack(fids[w],2)
                    assert fold.df_advance(23,sf,2,bit,out)==0
                    assert int.from_bytes(bytes(sf),'little')==fids[new] and fraction(out)==nz['out'];foldchecks+=1
        row={'p':p,'states':len(ns),'construction_bits':truebits,'coding_words':co.size,'coding_bits':co.bits,
             'grid_states':co.J,'budget':co.R,'DP_entries':len(co.counts),'DP_bytes_uint64':len(co.counts)*8,
             'transition_data_bytes_uint32':16*co.J,'context_bytes':C.sizeof(Context),
             'build_seconds_local':build,'max_actual_word_cost':max(co.word_cost(w) for w in byword)}
        configs.append(row);print(json.dumps(row),flush=True)
    # All unused 12-bit IDs are rejected by the p23 checked output API.
    valid=set(fids.values());rejected=0
    for i in range(4096):
        rc=fold.df_output(23,pack(i,2),2,out)
        assert (rc==0)==(i in valid)
        if rc:rejected+=1
    assert rejected==4096-3707
    result={'passed':True,'configs':configs,'uniform_transitions':transitions,'uniform_exact_outputs':outputs,
            'all_coding_language_roundtrips':roundtrips,'fold23_transitions':foldchecks,'fold23_unused_ids_rejected':rejected,
            'structural':structural,'full_python_C_DP_comparisons':[4,6,8],'rolling_counter_comparisons':[4,6,8]}
    (ROOT/'results/VERIFIED.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASSED',transitions,outputs,roundtrips,foldchecks,flush=True)
if __name__=='__main__':main()
