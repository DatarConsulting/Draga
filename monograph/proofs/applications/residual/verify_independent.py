#!/usr/bin/env python3
"""Independent Python packet parser with a separate scalar FMA evaluator."""
import ctypes,json,struct,subprocess
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parent
CODE=r'''#include <math.h>
#include <stddef.h>
void eval(const float *u,const float *v,int rank,float *out){
 for(int i=0;i<128*128;i++)out[i]=0;
 for(int a=0;a<rank;a++)for(int i=0;i<128;i++)for(int j=0;j<128;j++)
  out[i*128+j]=fmaf(u[a*128+i],v[a*128+j],out[i*128+j]);
}
'''
class Bits:
    def __init__(self,b):self.b=b;self.pos=0
    def read(self,n):
        if not n:return 0
        p=self.pos//8;shift=self.pos%8;length=(n+shift+7)//8
        assert p+length<=len(self.b)
        v=int.from_bytes(self.b[p:p+length],'little')>>shift
        self.pos+=n;return v&((1<<n)-1)
def main():
    (ROOT/'src/independent_eval.c').write_text(CODE)
    subprocess.run(['gcc','-O2','-march=native','-ffp-contract=off','-fno-fast-math','-shared','-fPIC',str(ROOT/'src/independent_eval.c'),'-lm','-o',str(ROOT/'independent_eval.so')],check=True)
    lib=ctypes.CDLL(str(ROOT/'independent_eval.so'));ptr=np.ctypeslib.ndpointer(dtype=np.float32,flags='C_CONTIGUOUS')
    lib.eval.argtypes=[ptr,ptr,ctypes.c_int,ptr]
    def factors(blob):
        e=struct.unpack_from('<h',blob)[0];c=np.frombuffer(blob[2:258],dtype='<i2').astype(np.float32)
        assert -126<=e<=100
        return np.ldexp(c,e)
    reference=np.fromfile(ROOT/'inputs/atelier/native_state.f32',dtype='<u4').reshape(16,16384)
    tested=[]
    for name in ['energy_96','hankel_96','cross_32','row_ref']:
        data=(ROOT/'packets'/f'{name}_atelier_native.frame').read_bytes();assert data[:8]==b'HANKEL01'
        kind,rank,heads,basis_bytes=struct.unpack_from('<4I',data,8);assert heads==16
        basis=np.frombuffer(data[128:128+basis_bytes],dtype='<f4').copy()
        position=128+basis_bytes
        for h in range(16):
            length=struct.unpack_from('<I',data,32+4*h)[0];p=data[position:position+length];position+=length
            assert p[0]==0xa7 and p[1]==kind and struct.unpack_from('<I',p,12)[0]==length
            r,b,flag=p[2:5];ne=struct.unpack_from('<H',p,6)[0];nm=struct.unpack_from('<I',p,8)[0]
            meta=p[16:16+nm];pred=np.zeros(16384,dtype=np.float32)
            if kind==0:
                assert nm==512
                pred=np.repeat(np.frombuffer(meta,dtype='<f4'),128).copy()
            else:
                U=[];V=[]
                for a in range(r):
                    if kind==1:U.append(basis.reshape(rank,128)[a]);V.append(factors(meta[a*258:(a+1)*258]))
                    else:U.append(factors(meta[a*516:a*516+258]));V.append(factors(meta[a*516+258:(a+1)*516]))
                u=np.asarray(U,dtype=np.float32).reshape(-1);v=np.asarray(V,dtype=np.float32).reshape(-1)
                lib.eval(u,v,r,pred)
            start=16+nm;signs=p[start:start+2048];start+=2048
            lows=Bits(p[start:start+2048*b]);start+=2048*b
            bitmap=p[start:start+2048] if flag else b'';start+=2048 if flag else 0
            high=Bits(p[start:]);words=pred.view(np.uint32);back=np.empty_like(words);count=0
            for i in range(16384):
                z=lows.read(b)
                if flag and (bitmap[i//8]>>(i%8))&1:z|=high.read(32-b)<<b;count+=1
                d=-(z//2)-1 if z&1 else z//2
                a=int(words[i]&0x7fffffff)+d;assert 0<=a<=0x7fffffff
                back[i]=a|(((signs[i//8]>>(i%8))&1)<<31)
            assert count==ne and np.array_equal(back,reference[h])
        assert position==len(data)
        tested.append({'method':name,'heads':16,'words':262144,'exact':True,'frame_bytes':len(data)})
    result={'passed':True,'independent_parser':True,'shared_predictor_evaluation':'Separate C loop with scalar fmaf; no archived/new codec header included','frames':tested,'total_words':sum(x['words'] for x in tested)}
    (ROOT/'results/independent_checks.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
if __name__=='__main__':main()
