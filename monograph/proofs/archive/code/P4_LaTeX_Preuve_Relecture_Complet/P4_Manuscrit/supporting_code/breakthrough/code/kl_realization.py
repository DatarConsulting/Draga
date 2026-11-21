"""Same-state local smoothing, exact C output, rational TV and KL upper bounds.
KL theorem: D(p||q_s)<=-ln(1-d)-ln(1-s)/tau
  +d ln(2/s)+Psi_p(d) ln(2m/s), if TV(p,q)<=d.
Proof supplied in the report; log upper bounds are rational with certified tails.
"""
import json,math,ctypes,subprocess,hashlib
from pathlib import Path
from fractions import Fraction as F
from probability import pstep
ROOT=Path(__file__).resolve().parents[1]
def logup(z):
    assert z>=1
    k=0
    while z>=2:z/=2;k+=1
    def series(x):
        u=(x-1)/(x+1);N=24
        return 2*sum((u**(2*j+1)/F(2*j+1) for j in range(N)),F(0))+2*u**(2*N+1)/(F(2*N+1)*(1-u*u))
    return k*series(F(2))+series(z)
def round_up(x,bits=40):return F((x.numerator*2**bits+x.denominator-1)//x.denominator,2**bits)
def run():
    subprocess.run(['gcc','-std=c11','-O2','-shared','-fPIC','-fstack-usage',str(ROOT/'code/belief.c'),'-o',str(ROOT/'code/belief.so')],check=True)
    lib=ctypes.CDLL(str(ROOT/'code/belief.so'));fn=lib.belief_probabilities_smoothed
    fn.argtypes=[ctypes.POINTER(ctypes.c_uint8),ctypes.c_int,ctypes.c_uint32,ctypes.c_int,ctypes.POINTER(ctypes.c_uint32),ctypes.c_int,ctypes.POINTER(ctypes.c_uint64)];fn.restype=ctypes.c_uint64
    data=json.loads((ROOT/'results/probability.json').read_text());records=[]
    for c in data['synthetic']+data['trained']:
        tv=F(c['total_word_TV_bound']);tau=F(c['stopping_min']);tau0=F(c['tail_min_stopping']);d=round_up(tv,16)
        s=1
        while F(1,2**s)>(F(1,100)-tv)*tau/2:s+=1
        lam=F(1,2**s);totalTV=tv+lam/tau;assert totalTV<=F(1,100)
        r=1-tau0;N=max(0,math.ceil(math.log(float(d))/math.log(float(r))));psi=d*N+r**(N+1)/tau0
        m=c['alphabet_size'];kl=logup(1/(1-d))+logup(1/(1-lam))/tau0+d*logup(2/lam)+psi*logup(2*m/lam)
        KL=round_up(kl);assert KL>=kl
        A=c['coefficient_integers'];flat=[v for row in A for v in row];arr=(ctypes.c_uint32*len(flat))(*flat);G=2**c['belief_fractional_bits']
        maxden=G*2**24*2**s*(2*m);assert maxden<2**64
        seen={G//2};todo=[G//2];checks=0
        while todo:
            j=todo.pop();raw=bytearray(j.to_bytes(c['state_bytes'],'little'));buf=(ctypes.c_uint8*len(raw)).from_buffer(raw);nums=(ctypes.c_uint64*(m+1))()
            den=fn(buf,len(raw),G,m,arr,s,nums);assert den==maxden and sum(nums)==den and min(nums)>0
            for a in range(m+1):
                old=j*(sum(A[0][2*a:2*a+2]) if a<m else A[0][-1])+(G-j)*(sum(A[1][2*a:2*a+2]) if a<m else A[1][-1])
                assert F(nums[a],den)==(1-lam)*F(old,G*2**24)+lam*(F(1,2*m) if a<m else F(1,2));checks+=1
                if a<m:
                    nxt=pstep(j,a,A,G)
                    if nxt not in seen:seen.add(nxt);todo.append(nxt)
        original=ROOT/'models/probability'/f"{c['id']}.bin";oldpacket=original.read_bytes()
        # The last byte(s) are initial state; insert the smoothing exponent into static metadata.
        n=c['state_bytes'];packet=b'PFS1'+oldpacket[4:-n]+bytes([s])+oldpacket[-n:]
        (ROOT/'models/probability'/f"{c['id']}_KL.bin").write_bytes(packet)
        records.append({'id':c['id'],'lambda':str(lam),'lambda_exponent':s,'base_TV':str(tv),'smoothed_TV_bound':str(totalTV),
          'd_for_KL':str(d),'tail_stopping':str(tau0),'tail_cutoff':N,'tail_first_moment_bound':str(psi),
          'KL_rational_upper':str(KL),'KL_nats_upper_display':float(KL),'exact_C_probabilities_checked':checks,
          'maximum_unsigned_arithmetic_bits':maxden.bit_length(),'state_bytes':n,'packet_bytes':len(packet),
          'packet_sha256':hashlib.sha256(packet).hexdigest()})
    (ROOT/'results/KL_realization.json').write_text(json.dumps({'all_passed':True,'records':records,
      'scope':'Locally smoothed machine, not the global geometric mixture. Same persistent state; one extra static byte; output integer arithmetic. KL against the reference word law, not language data distribution.'},indent=2))
    print('PASS KL realization',len(records),'range',min(x['KL_nats_upper_display'] for x in records),max(x['KL_nats_upper_display'] for x in records))
if __name__=='__main__':run()
