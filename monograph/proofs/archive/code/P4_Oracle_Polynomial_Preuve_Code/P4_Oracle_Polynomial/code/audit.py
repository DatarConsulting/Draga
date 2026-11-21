"""Exact references and enumeration occur only in this test/benchmark file."""
from spectral_codec import SpectralCodec
from spectral_oracle import coordinates,interval_matrix,verify_polynomial
from reduced_codec import ReducedCodec
from certificate import certify,hmm,point
from model import Model
from dyadic import Dyadic
from fractions import Fraction as F
from pathlib import Path
from math import comb
import json,time,random,statistics,sys

ROOT=Path(__file__).resolve().parents[1]
CASES=[('S',hmm('1/4','1/4','6/25','19/25')),
       ('A',hmm('6/25','13/50','6/25','3/4'))]

def save(name,x):
    (ROOT/'results'/name).write_text(json.dumps(x,indent=2)+'\n')

def make(ms,p):return Model.from_rationals(p,certify(ms,J=(F(1,4),F(3,4))))
def mul(x,y):
    a,b,c,d=x;A,B,C,D=y
    return a*A+b*C,a*B+b*D,c*A+d*C,c*B+d*D
def out(ms,x):
    a,b,c,d=ms[0]
    return (a+c)*x+(b+d)*(1-x)

def local_checks():
    checks=0
    for _,ms in CASES:
        mod=make(ms,16);ctx=Dyadic(128);m=12
        M=interval_matrix(mod,F(1),m,ctx)
        maps=coordinates(mod)[2]
        for k in range(m+1):
            for j in range(m+1):
                exact=F(0)
                for a,b,xi,kap in maps:
                    exact+=kap*sum(F(comb(j,l))*a**(j-l)*b**l *
                                   ((-1)**(k-l))*comb(j+k-l+1,k-l)*xi**(k-l)
                                   for l in range(min(j,k)+1))
                assert F(M[k][j][0],ctx.S)<=exact<=F(M[k][j][1],ctx.S)
                checks+=1
    # Actual nominal S map: a point inside the proposed disk already leaves D_1.
    z=-F(5,4);alpha=-F(52,75);beta=F(1,2);xi=-F(39,200)
    image=(alpha+beta*z)/(1+xi*z)
    H=certify(CASES[0][1],J=(F(1,4),F(3,4)))
    assert z*z*H.taylor_radius<1 and abs(image)>1
    # Exact non-commutation at m=0, n=2, sigma=1.
    kap=F(228,625)
    V1=2*kap
    G0=lambda x:(alpha+beta*x)/(1+xi*x)
    G1=lambda x:(-alpha+beta*x)/(1-xi*x)
    g0=lambda x:kap/(1+xi*x)**2
    g1=lambda x:kap/(1-xi*x)**2
    V2=kap*(g0(G0(0))+g1(G0(0))+g0(G1(0))+g1(G1(0)))
    assert V2!=V1**2
    N=3;k=3;Q=4
    left=(2*(N-1)*k+Q)//(2*Q)
    assert N-1-left==0
    data=dict(passed=True,exact_matrix_entries=checks,
        disk_counterexample=dict(z=str(z),image=str(image),r=str(H.taylor_radius)),
        repeated_truncation_counterexample=dict(V2_center=str(V2),matrix_square_center=str(V1**2),difference=str(V2-V1**2)),
        capacity_equality_counterexample=dict(N=N,q='3/4',left=left,right=N-1-left))
    save('LOCAL_CHECKS.json',data)
    return data

def graph(name,ms,p):
    model=make(ms,p);c=SpectralCodec(model);H=model.cert;eps=F(1,1<<p)
    I=(F(1),F(0),F(0),F(1));rows=[];cache={};trueprob=c.probability
    stored=tuple(tuple(F(x,model.H) for x in M) for M in model.mats)
    bmin=F(c.capacity_certificate['rounded_branch_min_lower'])
    loss=F(c.capacity_certificate['worst_capacity_loss'])
    def visit(w,z,actual,quant,lo,N,mu):
        leaf=model.leaf(z);xs=H.J if leaf else (F(1,2),)
        tv=max(abs(F(model.output(z),model.Q)-out(ms,point(actual,x))) for x in xs)
        assert tv<=eps
        zl,zr=model.enclosure(z);el,er=point(quant,H.J[0]),point(quant,H.J[1])
        assert F(zl,model.Q)<=el<=er<=F(zr,model.Q)
        assert max(el-F(zl,model.Q),F(zr,model.Q)-er)<=model.point_error_bound
        assert N>=F(1<<c.B)*mu-loss and N>=1
        rows.append(dict(w=w,z=z,actual=actual,lo=lo,N=N,leaf=leaf,tv=tv))
        if not leaf:
            assert F(1<<c.B)*mu>2/bmin
            k,Q=trueprob(z);cache[z]=(k,Q);q=F(k,Q)
            assert bmin<=q<=1-bmin
            left=c.split(N,k,Q)
            for bit,size,shift,prob in [(0,left,1,q),(1,N-1-left,1+left,1-q)]:
                visit(w+str(bit),model.step(z,bit),mul(actual,ms[bit]),mul(quant,stored[bit]),lo+shift,size,mu*prob)
    visit('',model.root(),I,I,0,1<<c.B,F(1))
    # Test-only cache: it avoids timing the same pure branch evaluation again.
    c.probability=lambda z:cache[z]
    words={r['w']:r for r in rows};ids={r['lo'] for r in rows}
    assert len(ids)==len(rows)
    transitions=0
    for r in rows:
        assert c.encode(r['w'])==(r['lo'],r['N'],r['z'])
        assert c.decode(r['lo'])[0]==r['w']
        for bit in (0,1):
            raw=str(bit)+r['w'];dest=raw
            for j in range(len(raw)+1):
                if raw[:j] in words and words[raw[:j]]['leaf']:
                    dest=raw[:j];break
            target=words[dest]
            assert c.advance(r['lo'],bit)==target['lo']
            if r['leaf']:assert target['leaf']
            for x in H.J if r['leaf'] else (F(1,2),):
                y=point(ms[bit],point(r['actual'],x))
                assert point(target['actual'],H.J[0])<=y<=point(target['actual'],H.J[1])
                assert abs(F(model.output(target['z']),model.Q)-out(ms,y))<=eps
            transitions+=1
    return dict(model=name,p=p,states=len(rows),transitions=transitions,bits=c.B,
                offset=c.capacity_certificate['offset_bits'],
                max_TV_over_epsilon=float(max(r['tv'] for r in rows)/eps),
                test_only_state_enumeration=True)

def configurations():
    import oracle as old
    original={name:getattr(old,name) for name in ['transfer','transfer_coefficients']}
    def forbidden(*a,**kw):raise AssertionError('Forbidden enumerating oracle called')
    for name in original:setattr(old,name,forbidden)
    rows=[]
    try:
        for name,ms in CASES:
            for p in (4,8,12,16,32,64,128):
                model=make(ms,p);start=time.perf_counter();c=SpectralCodec(model);elapsed=time.perf_counter()-start
                # Bad candidates must not inherit a certificate from a good one.
                bad=list(c.oracle.coeffs);bad[1]+=10
                assert verify_polynomial(model,c.sigma,bad) is None
                state=0;x=F(1,2);peak=F(0);eps=F(1,1<<p);rng=random.Random(3000+p)
                for step in range(64):
                    b=(step%2) if step<16 else rng.randrange(2)
                    state=c.advance(state,b);x=point(ms[b],x)
                    error=abs(F(c.output(state),model.Q)-out(ms,x))
                    assert error<=eps;peak=max(peak,error)
                row=dict(model=name,p=p,bits=c.B,offset=c.capacity_certificate['offset_bits'],
                         compilation_seconds=elapsed,max_TV_over_epsilon=float(peak/eps),
                         stream_steps=64,oracle=c.oracle.json(),capacity=c.capacity_certificate)
                rows.append(row)
                print('config',name,p,'bits',c.B,'seconds',round(elapsed,4),flush=True)
                save('CONFIGURATIONS.json',dict(passed=len(rows)==14,rows=rows,
                    enumerating_oracles_patched_to_raise=True,production_state_table=False,
                    universal_R6_closed=False))
    finally:
        for name,f in original.items():setattr(old,name,f)

def graphs():
    rows=[]
    for name,ms in CASES:
        for p in (4,8,12):
            row=graph(name,ms,p);rows.append(row);print('graph',row,flush=True)
            save('GRAPHS.json',dict(passed=len(rows)==6,rows=rows,
                 states=sum(x['states'] for x in rows),transitions=sum(x['transitions'] for x in rows)))

def benchmarks():
    rows=[]
    for name,ms in CASES:
        for p in (16,64,128):
            model=make(ms,p);rng=random.Random(8128);bits=[rng.randrange(2) for _ in range(32)]
            row=dict(model=name,p=p,steps_per_trial=len(bits),trials=3)
            for label,cls in [('polynomial',SpectralCodec),('enumeration',ReducedCodec)]:
                start=time.perf_counter();codec=cls(model);ct=time.perf_counter()-start
                times=[]
                for repeat in range(3):
                    state=0;start=time.perf_counter()
                    for b in bits:state=codec.advance(state,b);codec.output(state)
                    times.append(time.perf_counter()-start)
                row[label]=dict(compile_seconds=ct,step_output_ms=1000*statistics.median(times)/len(bits),
                    bits=codec.B,horizon=codec.n,all_trial_seconds=times)
            row['observed_speedup']=row['enumeration']['step_output_ms']/row['polynomial']['step_output_ms']
            rows.append(row);print('bench',row,flush=True)
            save('BENCHMARKS.json',dict(passed=len(rows)==6,rows=rows,
                 timing_scope='CPython, uncached production advance+output; exact verifier excluded',
                 asymptotic_complexity_validated_by_benchmarks=False))

if __name__=='__main__':
    mode=sys.argv[1] if len(sys.argv)>1 else 'all'
    if mode in ('local','all'):print(local_checks(),flush=True)
    if mode in ('configs','all'):configurations()
    if mode in ('graphs','all'):graphs()
    if mode in ('bench','all'):benchmarks()
