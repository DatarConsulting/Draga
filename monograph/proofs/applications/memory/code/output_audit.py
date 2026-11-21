"""Exact C output check over complete reachable graphs, plus table-state baseline."""
import json,ctypes,subprocess,struct,hashlib,math
from pathlib import Path
from fractions import Fraction as F
import numpy as np
from replay import pack,unpack,py_step
ROOT=Path(__file__).resolve().parents[1]
def run():
    subprocess.run(['gcc','-std=c11','-O2','-shared','-fPIC','-fstack-usage',str(ROOT/'code/kernel.c'),'-o',str(ROOT/'code/kernel.so')],check=True)
    lib=ctypes.CDLL(str(ROOT/'code/kernel.so'));fn=lib.output_value
    ptr=ctypes.POINTER(ctypes.c_int32);bptr=ctypes.POINTER(ctypes.c_uint8)
    fn.argtypes=[bptr,ctypes.c_int,ctypes.c_int,ptr,ptr,ptr,ctypes.c_int,ctypes.POINTER(ctypes.c_int64),ctypes.POINTER(ctypes.c_int)]
    ts=lib.table_step;ts.argtypes=[bptr,ctypes.c_int,ctypes.c_int,bptr]
    rows=[];dest=ROOT/'models/tables';dest.mkdir(exist_ok=True)
    for path in sorted((ROOT/'models').glob('*.json')):
      m=json.loads(path.read_text())
      for policy in ['uniform','mixed']:
        c=m['best'][policy];f=m['frames'][c['basis']];bits=np.array(c['bits'],np.int32);sh=np.array(c['shifts'],np.int32)
        A=np.array([[[int(F(v)*2**24) for v in row] for row in a] for a in f['A']],np.int32)
        beta=np.array([int(F(v[0])*2**24) for v in f['beta']],np.int32)
        bound=sum(2**(int(bits[i])-1)*abs(int(beta[i]))*2**(int(sh[i])-int(min(sh))) for i in range(2));assert bound<2**63
        start=tuple(c['initial_integers']);states=[start];ids={start:0};nxt=[];outputs=[];e=None
        for s in states:
            trans=[]
            for a in range(2):
                t=tuple(py_step(s,A,bits.tolist(),sh.tolist(),a,24))
                if t not in ids:ids[t]=len(states);states.append(t)
                trans.append(ids[t])
            nxt.append(trans)
            raw=bytearray(pack(s,bits));buf=(ctypes.c_uint8*len(raw)).from_buffer(raw);num=ctypes.c_int64();exp=ctypes.c_int()
            assert fn(buf,len(raw),2,beta.ctypes.data_as(ptr),bits.ctypes.data_as(ptr),sh.ctypes.data_as(ptr),24,ctypes.byref(num),ctypes.byref(exp))==0
            expected=sum(F(s[i])*F(c['steps'][i])*F(f['beta'][i][0]) for i in range(2))
            assert F(num.value)*F(2)**exp.value==expected
            e=exp.value;outputs.append(num.value)
        # Moore minimization preserves exact dyadic output on every continuation.
        def classes(signatures):
            d={};out=[]
            for sig in signatures:
                if sig not in d:d[sig]=len(d)
                out.append(d[sig])
            return out
        cls=classes(outputs)
        while True:
            nc=classes([(outputs[i],cls[nxt[i][0]],cls[nxt[i][1]]) for i in range(len(states))])
            if nc==cls:break
            cls=nc
        N=max(cls)+1;assert N<=256 and cls[0]==0
        reps=[cls.index(i) for i in range(N)];tn=[[cls[nxt[r][a]] for a in range(2)] for r in reps];to=[outputs[r] for r in reps]
        table=(ctypes.c_uint8*(2*N))(*[x for row in tn for x in row])
        for i in range(len(states)):
            assert to[cls[i]]==outputs[i]
            for a in range(2):
                st=ctypes.c_uint8(cls[i]);assert ts(ctypes.byref(st),a,N,table)==0 and st.value==cls[nxt[i][a]]
        packet=struct.pack('<4sHh',b'TBL1',N,e)+bytes(table)+b''.join(struct.pack('<q',v) for v in to)
        (dest/f"{m['id']}_{policy}.tbin").write_bytes(packet+bytes([0]))
        rows.append({'id':m['id'],'policy':policy,'reachable_states':len(states),'minimal_equivalent_Moore_states':N,
          'information_bits_for_exact_reencoding':math.ceil(math.log2(N)),'table_persistent_bytes':1,
          'table_static_bytes':len(packet),'table_serialized_bytes':len(packet)+1,'original_state_bytes':c['state_bytes'],
          'original_serialized_bytes':64+c['state_bytes'],'exact_C_outputs_checked':len(states),'exact_C_table_edges_checked':2*len(states),
          'output_signed_accumulator_bits_bound':bound.bit_length()+1,'table_sha256':hashlib.sha256(packet+bytes([0])).hexdigest()})
    (ROOT/'results/output_and_table_audit.json').write_text(json.dumps({'all_passed':True,'records':rows,
      'scope':'Complete reachable graphs; exact integer output. Tables are behaviorally identical to the quantized machines, not to the unquantized source. One-byte table IDs remove the mixed-vs-uniform persistent-byte advantage at increased static storage.'},indent=2))
    print('PASS output + table',len(rows),flush=True)
if __name__=='__main__':run()
