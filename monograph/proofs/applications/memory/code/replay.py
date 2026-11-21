"""Execute the actual packed C recurrence; independently replay using integers."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import json,ctypes,subprocess,time,itertools,math,hashlib
from pathlib import Path
from fractions import Fraction as F
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
def read_matrix(x):return np.array([[float(F(v)) for v in row] for row in x])
def iround(n,d):
    q,r=divmod(n,d)
    return q+int(2*r>d or (2*r==d and q%2))
def pack(vals,bits):
    n=0;o=0
    for x,p in zip(vals,bits):
        x=int(x);p=int(p);n|=(x%2**p)<<o;o+=p
    return n.to_bytes(math.ceil(o/8),'little')
def unpack(raw,bits):
    n=int.from_bytes(raw,'little');o=0;out=[]
    for p in bits:
        p=int(p)
        a=(n>>o)&(2**p-1);out.append(a if a<2**(p-1) else a-2**p);o+=p
    return out
def py_step(s,A,bits,shift,a,Q):
    base=min(shift)-Q;out=[]
    for j in range(len(bits)):
        n=sum(int(s[i])*int(A[a,i,j])*2**(shift[i]-Q-base) for i in range(len(bits)))
        sh=shift[j]-base
        out.append(iround(n,2**sh) if sh>=0 else n*2**(-sh))
    return out
def run():
    so=ROOT/'code/kernel.so'
    subprocess.run(['gcc','-std=c11','-O2','-shared','-fPIC','-fstack-usage',str(ROOT/'code/kernel.c'),'-o',str(so)],check=True)
    lib=ctypes.CDLL(str(so));fn=lib.step
    fn.argtypes=[ctypes.POINTER(ctypes.c_uint8),ctypes.c_int,ctypes.c_int,ctypes.POINTER(ctypes.c_int32),ctypes.POINTER(ctypes.c_int32),ctypes.POINTER(ctypes.c_int32),ctypes.c_int,ctypes.c_int]
    def cstep(vals,bits,shift,A,tok):
        buf=(ctypes.c_uint8*len(vals)).from_buffer(vals)
        ret=fn(buf,len(vals),len(bits),A.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),bits.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),shift.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),tok,24)
        assert ret==0,ret
    rows=[];total_steps=0
    shard=int(os.environ.get('SHARD','0'));nshards=int(os.environ.get('NSHARDS','1'))
    for idx,path in enumerate(sorted((ROOT/'models').glob('*.json'))):
        if idx%nshards!=shard:continue
        m=json.loads(path.read_text())
        if 'best' not in m:continue
        for policy in ['uniform','mixed']:
            receipt=ROOT/'results'/f"replay_{m['id']}_{policy}.json"
            if receipt.exists():rows.append(json.loads(receipt.read_text()));continue
            alloc=m['best'][policy]
            if not alloc:continue
            f=m['frames'][alloc['basis']]
            bits=np.array(alloc['bits'],np.int32);shift=np.array(alloc['shifts'],np.int32)
            A=np.array([[[int(F(v)*2**24) for v in row] for row in a] for a in f['A']],np.int32)
            sourceA=np.array([read_matrix(a) for a in m['source_A']]);source0=read_matrix(m['source_alpha']).ravel();sourceb=read_matrix(m['source_beta']).ravel()
            beta=read_matrix(f['beta']).ravel();scales=2.**shift
            init=alloc['initial_integers'];initial=pack(init,bits)
            besterr=abs(float(source0@sourceb)-float(np.array(init)*scales@beta));comparisons=0;states={initial}
            # Every binary word <=16 through a prefix tree, each prefix once.
            level=[(initial,source0)]
            for length in range(1,17):
                nextlevel=[]
                for raw,source in level:
                    pystate=unpack(raw,bits)
                    for a in [0,1]:
                        buf=bytearray(raw);cstep(buf,bits,shift,A,a)
                        exact=py_step(pystate,A,bits.tolist(),shift.tolist(),a,24)
                        assert bytes(buf)==pack(exact,bits)
                        newsource=source@sourceA[a]
                        error=abs(float(newsource@sourceb)-float(np.array(exact)*scales@beta))
                        besterr=max(besterr,error);comparisons+=1;states.add(bytes(buf))
                        nextlevel.append((bytes(buf),newsource))
                level=nextlevel
            # Deterministic long sequences: 100k steps each, including common suffix pair.
            rng=np.random.default_rng(991)
            seqs=[np.zeros(100000,dtype=np.int8),np.ones(100000,dtype=np.int8),np.arange(100000,dtype=np.int32)%2,rng.integers(0,2,100000,dtype=np.int8)]
            common=rng.integers(0,2,99998,dtype=np.int8)
            seqs.extend([np.r_[0,1,common],np.r_[1,0,common]])
            for seq in seqs:
                raw=bytearray(initial);s=list(init);source=source0.copy()
                for a in seq:
                    cstep(raw,bits,shift,A,int(a));s=py_step(s,A,bits.tolist(),shift.tolist(),int(a),24)
                    assert bytes(raw)==pack(s,bits)
                    source=source@sourceA[a]
                    err=abs(float(source@sourceb)-float(np.array(s)*scales@beta));besterr=max(besterr,err);comparisons+=1
            assert besterr<=float(F(alloc['total_bound']))+1e-11
            # Exact reachability closure of the finite machine, no sequence-length cutoff.
            seen={initial};front=[initial]
            while front:
                raw=front.pop()
                for a in [0,1]:
                    buf=bytearray(raw);cstep(buf,bits,shift,A,a);b=bytes(buf)
                    if b not in seen:
                        seen.add(b);front.append(b)
                if len(seen)>1000000:raise RuntimeError('reachable-state cap')
            total_steps+=comparisons
            # Serialized deployed coefficient metadata and initial state.
            binary=b'R2Q1'+bytes([2,24])+bytes(bits.tolist())+shift.astype('<i4').tobytes()+A.astype('<i4').tobytes()
            binary+=np.array([int(F(v)*2**24) for row in f['alpha'] for v in row],dtype='<i4').tobytes()
            binary+=np.array([int(F(v)*2**24) for row in f['beta'] for v in row],dtype='<i4').tobytes()
            packet=ROOT/'models'/f"{m['id']}_{policy}.bin";packet.write_bytes(binary+initial)
            row={'id':m['id'],'policy':policy,'state_bytes':len(initial),'static_bytes':len(binary),'serialized_bytes':len(binary)+len(initial),
                 'exact_C_python_steps':comparisons,'reachable_machine_states':len(seen),'max_observed_abs_error':besterr,
                 'certified_bound':float(F(alloc['total_bound'])),'source_reference':'float64 diagnostic; exact proof is separate',
                 'packet_sha256':hashlib.sha256(packet.read_bytes()).hexdigest()}
            rows.append(row);print(json.dumps(row),flush=True)
            receipt.write_text(json.dumps(row,indent=2))
    (ROOT/f'results/replay_summary_{shard}.json').write_text(json.dumps({'records':rows,'total_exact_step_comparisons':sum(r['exact_C_python_steps'] for r in rows),'overflow_events':0,'mismatches':0},indent=2))
if __name__=='__main__':run()
