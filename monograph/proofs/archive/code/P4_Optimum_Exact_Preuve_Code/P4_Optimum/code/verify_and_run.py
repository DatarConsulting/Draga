"""Exact tests of the rank-two memory theorem and full C transition graphs.
Proof is in theorem.tex; the finite checks do not replace the quantified proof.
No imported producer certificate or floating-point acceptance test.
"""
import json,ctypes,subprocess,math,hashlib,itertools,time,platform
from pathlib import Path
from fractions import Fraction as F
import sympy as sp
ROOT=Path(__file__).resolve().parents[1]
def spec(a,m):
    L=(a+1)*m+2;T=(L-2+a)//a
    widths=[min(t,L-1-a*t) for t in range(T)]
    offsets=[0]
    for n in widths:offsets.append(offsets[-1]+2**n)
    return L,T,widths,offsets,offsets[-1]+1
def decode(i,widths,offsets):
    if i==offsets[-1]:return len(widths),0
    for t,n in enumerate(widths):
        if i<offsets[t+1]:return t,i-offsets[t]
    raise AssertionError('invalid state')
def expected_step(t,b,token,T,widths,offsets):
    if t>=T-1:return offsets[-1]
    n=widths[t];nn=widths[t+1]
    # Independent reference: prepend the newest bit to a textual suffix.
    old=format(b,f'0{n}b') if n else ''
    new=(str(token)+old)[:nn]
    return offsets[t+1]+int(new,2)
def interval(a,t,n,b):
    lo=F(b,2**(a*t+n))
    hi=lo+F(2**(t-n)-1,2**((a+1)*t))
    return lo,hi
def lower_small(a,m):
    rho=F(1,2**a);eps=F(1,2**((a+1)*m+2));pts=[]
    for t in range(m+1):
        for b in range(2**t):pts.append((rho**t,F(b,2**((a+1)*t))))
    t=m+1
    for b in [0,2**t-1]:pts.append((rho**t,F(b,2**((a+1)*t))))
    mingap=None;pairs=0
    for i,(z,x) in enumerate(pts):
        for zz,xx in pts[:i]:
            d=max(abs(x-xx),rho*abs(x+z-xx-zz)/2)
            assert d>2*eps
            mingap=d if mingap is None else min(mingap,d);pairs+=1
    return {'points':len(pts),'pairs':pairs,'minimum_gap':str(mingap),'threshold':str(2*eps)}
def rank_stability(a):
    rho=sp.Rational(1,2**a);A0=sp.diag(rho,rho/2);A1=sp.Matrix([[rho,rho/2],[0,rho/2]])
    al=sp.Matrix([[1,0]]);be=sp.Matrix([0,1]);As=[A0,A1]
    K=sp.eye(4)-sum((sp.kronecker_product(A,A) for A in As),sp.zeros(4))
    P=sp.Matrix(2,2,list(K.inv()*sp.Matrix(list(be*be.T))))
    Kt=sp.eye(4)-sum((sp.kronecker_product(A.T,A.T) for A in As),sp.zeros(4))
    Q=sp.Matrix(2,2,list(Kt.inv()*sp.Matrix(list(al.T*al))))
    X=sp.Matrix(2,2,list(K.inv()*sp.Matrix(list(sp.eye(2)))))
    for M in [P,Q,X]:assert M==M.T and M[0,0]>0 and M.det()>0
    assert X-sum((A*X*A.T for A in As),sp.zeros(2))==sp.eye(2)
    assert P-sum((A*P*A.T for A in As),sp.zeros(2))==be*be.T
    assert Q-sum((A.T*Q*A for A in As),sp.zeros(2))==al.T*al
    H=sp.Matrix([[0,rho/2],[rho/2,3*rho*rho/4]])
    assert H.det()==-rho*rho/4!=0 and A0*A1!=A1*A0
    hs=sp.trace(P*Q);assert hs<=sp.Rational(292,441)<1
    return {'a':a,'rho':str(rho),'H_minor_det':str(H.det()),'HS_squared':str(hs),
      'P':[[str(v) for v in row] for row in P.tolist()], 'Q':[[str(v) for v in row] for row in Q.tolist()],
      'X':[[str(v) for v in row] for row in X.tolist()]}
def run():
    start=time.time();root=ROOT/'code';subprocess.run(['gcc','-O2','-std=c11','-shared','-fPIC','-fstack-usage',str(root/'optimal_memory.c'),'-o',str(root/'optimal_memory.so')],check=True)
    lib=ctypes.CDLL(str(root/'optimal_memory.so'));bs=ctypes.POINTER(ctypes.c_uint8)
    step=lib.p4_step;step.argtypes=[bs,ctypes.c_int,ctypes.c_int,ctypes.c_int,ctypes.c_int]
    output=lib.p4_output;output.argtypes=[bs,ctypes.c_int,ctypes.c_int,ctypes.c_int,ctypes.POINTER(ctypes.c_uint64),ctypes.POINTER(ctypes.c_int)]
    records=[];ranks=[];small=[]
    for a in range(1,7):
      ranks.append(rank_stability(a))
      for m in range(a+2,13):
        if (a+1)*m+2>62:continue
        L,T,widths,offsets,N=spec(a,m);eps=F(1,2**L);rho=F(1,2**a);B=m+2;nb=(B+7)//8
        lower=2**(m+1)+1
        assert 2**(B-1)<lower<=N<2**B
        # Uniform analytic separation margins, covering all pairs in the witness.
        within=F(1,2**((a+1)*m));assert within==4*eps
        cross=rho/2*(rho**m*(1-rho)-2*eps);assert cross>2*eps
        extreme=F(2**(m+1)-1,2**((a+1)*(m+1)));assert extreme>2*eps
        if m<=5:small.append({'a':a,'m':m,**lower_small(a,m)})
        edges=outputs=0;err=F(0)
        for i in range(N):
            t,b=decode(i,widths,offsets);raw=bytearray(i.to_bytes(nb,'little'));buf=(ctypes.c_uint8*nb).from_buffer(raw)
            numerator=ctypes.c_uint64();power=ctypes.c_int();assert output(buf,nb,a,m,ctypes.byref(numerator),ctypes.byref(power))==0
            assert power.value==L;y=F(numerator.value,2**L);outputs+=1
            if t<T:
                lo,hi=interval(a,t,widths[t],b);err=max(err,abs(y-lo),abs(y-hi))
                assert abs(y-lo)<=eps and abs(y-hi)<=eps
            else:
                lo=F(0);hi=rho**T;assert hi<=2*eps and y==eps
            for token in [0,1]:
                raw=bytearray(i.to_bytes(nb,'little'));buf=(ctypes.c_uint8*nb).from_buffer(raw)
                assert step(buf,nb,a,m,token)==0
                dest=int.from_bytes(raw,'little');assert dest==expected_step(t,b,token,T,widths,offsets);edges+=1
                tt,bb=decode(dest,widths,offsets)
                if t<T:
                    image_lo=rho*(lo+token*rho**t)/2;image_hi=rho*(hi+token*rho**t)/2
                    if tt<T:
                        dlo,dhi=interval(a,tt,widths[tt],bb);assert dlo<=image_lo<=image_hi<=dhi
                    else:assert 0<=image_lo<=image_hi<=rho**T
                else:assert dest==N-1
        assert err==eps
        # Explicit serialization: eight static bytes plus packed state, no clock field.
        packet=b'P4B1'+bytes([a,m,B,nb])+bytes(nb)
        (ROOT/'models'/f'a{a}_m{m}.bin').write_bytes(packet)
        row={'a':a,'m':m,'rho':str(rho),'epsilon':str(eps),'epsilon_power':L,
          'lower_witness_states':lower,'constructed_states':N,'optimal_persistent_bits':B,
          'persistent_bytes':nb,'static_parameter_bytes':8,'packet_bytes':len(packet),'horizon_to_sink':T,
          'exact_C_edges_checked':edges,'exact_C_outputs_checked':outputs,'worst_interval_error':str(err),
          'cross_depth_separation_margin':str(cross-2*eps),'packet_sha256':hashlib.sha256(packet).hexdigest()}
        records.append(row);print(json.dumps(row),flush=True)
    # Two deliberate violations: false bit budget and output perturbation at a tight cell.
    controls={'one_less_bit_budget_rejected':all(r['lower_witness_states']>2**(r['optimal_persistent_bits']-1) for r in records),
      'output_above_tight_midpoint_rejected':F(1,256)+F(1,512)>F(1,256),
      'zero_output_tail_too_early_rejected':F(31,1024)>F(1,256)}
    assert all(controls.values())
    result={'all_passed':True,'records':records,'rank_certificates':ranks,'small_pairwise_witnesses':small,
      'negative_controls':controls,'total_C_edges':sum(r['exact_C_edges_checked'] for r in records),
      'total_C_outputs':sum(r['exact_C_outputs_checked'] for r in records),'elapsed_seconds':time.time()-start,
      'contract':'Deterministic time-homogeneous Moore machine, all persistent history and clock charged; arbitrary immutable code/tables allowed; worst-case absolute output error on all binary words.',
      'scope':'Exact bit optimum for the specified rank-two family and tolerances, not global optimization of the earlier 36 sources; compiler not formally verified.'}
    (ROOT/'results/exact_results.json').write_text(json.dumps(result,indent=2))
    print('PASS',len(records),result['total_C_edges'],result['total_C_outputs'],flush=True)
if __name__=='__main__':run()
