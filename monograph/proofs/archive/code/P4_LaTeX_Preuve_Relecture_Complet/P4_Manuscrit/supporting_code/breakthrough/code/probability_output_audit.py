import json,ctypes,subprocess,struct,hashlib
from pathlib import Path
from probability import pstep
ROOT=Path(__file__).resolve().parents[1]
def run():
    subprocess.run(['gcc','-std=c11','-O2','-shared','-fPIC','-fstack-usage',str(ROOT/'code/belief.c'),'-o',str(ROOT/'code/belief.so')],check=True)
    lib=ctypes.CDLL(str(ROOT/'code/belief.so'));fn=lib.belief_probabilities
    fn.argtypes=[ctypes.POINTER(ctypes.c_uint8),ctypes.c_int,ctypes.c_uint32,ctypes.c_int,ctypes.POINTER(ctypes.c_uint32),ctypes.POINTER(ctypes.c_uint64)];fn.restype=ctypes.c_uint64
    d=json.loads((ROOT/'results/probability.json').read_text());records=[];dest=ROOT/'models/probability';dest.mkdir(exist_ok=True)
    for c in d['synthetic']+d['trained']:
        K=c['alphabet_size'];G=2**c['belief_fractional_bits'];A=c['coefficient_integers'];ints=[v for row in A for v in row];arr=(ctypes.c_uint32*len(ints))(*ints)
        states={G//2};todo=[G//2];checks=0
        while todo:
            j=todo.pop();raw=bytearray(j.to_bytes(c['state_bytes'],'little'));buf=(ctypes.c_uint8*len(raw)).from_buffer(raw);nums=(ctypes.c_uint64*(K+1))()
            den=fn(buf,len(raw),G,K,arr,nums);assert den==G*2**24 and sum(nums)==den
            for a in range(K):
                assert nums[a]==j*sum(A[0][2*a:2*a+2])+(G-j)*sum(A[1][2*a:2*a+2])
                nxt=pstep(j,a,A,G)
                if nxt not in states:states.add(nxt);todo.append(nxt)
            assert nums[K]==j*A[0][-1]+(G-j)*A[1][-1];checks+=K+1
        static=struct.pack('<4sHBB',b'PFA1',K,c['belief_fractional_bits'],24)+struct.pack('<'+'I'*len(ints),*ints)
        packet=static+(G//2).to_bytes(c['state_bytes'],'little');(dest/(c['id']+'.bin')).write_bytes(packet)
        records.append({'id':c['id'],'exact_C_output_probabilities_checked':checks,'persistent_bytes':c['state_bytes'],'static_bytes':len(static),'packet_bytes':len(packet),'packet_sha256':hashlib.sha256(packet).hexdigest(),'reachable_states':len(states)})
    (ROOT/'results/probability_output_audit.json').write_text(json.dumps({'all_passed':True,'records':records},indent=2));print('PASS probability output',len(records))
if __name__=='__main__':run()
