"""Exact rational verifier. Full-tree enumeration and memoization are TEST ONLY."""
from reduced_codec import ReducedCodec,positive_power_lower
from compile_comparison import ROOT,CASES
from certificate import certify,point
from model import Model
from fractions import Fraction as F
from dyadic import Dyadic
from oracle import power
import json,random,time

def save(name,data):
    (ROOT/'results'/name).write_text(json.dumps(data,indent=2)+'\n')
def mul(x,y):
    a,b,c,d=x;A,B,C,D=y
    return a*A+b*C,a*B+b*D,c*A+d*C,c*B+d*D
def out(ms,x):
    a,b,c,d=ms[0];return (a+c)*x+(b+d)*(1-x)

def graph(name,ms,p):
    start=time.perf_counter();H=certify(ms,J=(F(1,4),F(3,4)))
    m=Model.from_rationals(p,H);codec=ReducedCodec(m);realprob=codec.probability
    eps=F(1,1<<p);I=(F(1),F(0),F(0),F(1));rows=[];cache={}
    stored=tuple(tuple(F(x,m.H) for x in mat) for mat in m.mats)
    cert=codec.capacity_certificate
    bmin=F(cert['rounded_branch_min_lower']);mass=F(cert['internal_mass_coefficient_lower'])
    loss=F(cert['worst_capacity_loss'])
    ctx=Dyadic(192);epspower=power(ctx,ctx.rat(eps),codec.sigma)
    mass_bound=mass*F(epspower[1],ctx.S)
    def visit(w,z,true,compiled,lo,N,mu):
        leaf=m.leaf(z);tv=F(0)
        endpoints=H.J if leaf else (F(1,2),)
        for x in endpoints:tv=max(tv,abs(F(m.output(z),m.Q)-out(ms,point(true,x))))
        assert tv<=eps
        zl,zr=m.enclosure(z)
        exact_l,exact_r=point(compiled,H.J[0]),point(compiled,H.J[1])
        assert F(zl,m.Q)<=exact_l<=exact_r<=F(zr,m.Q)
        assert max(exact_l-F(zl,m.Q),F(zr,m.Q)-exact_r)<=m.point_error_bound
        assert N>=1 and N>=F(1<<codec.B)*mu-loss
        row=dict(w=w,z=z,true=true,lo=lo,N=N,leaf=leaf,tv=tv,output=m.output(z))
        rows.append(row)
        if not leaf:
            assert mu>mass_bound
            assert F(1<<codec.B)*mu>2/bmin
            k,Q=realprob(z);cache[z]=(k,Q);q=F(k,Q)
            assert bmin<=q<=1-bmin
            # Independently compute nearest integer and both child sizes.
            left=((N-1)*k)//Q
            if 2*((N-1)*k-left*Q)>=Q:left+=1
            assert left==codec.split(N,k,Q) and 0<left<N-1
            for bit,size,shift,prob in [(0,left,1,q),(1,N-1-left,1+left,1-q)]:
                visit(w+str(bit),m.step(z,bit),mul(true,ms[bit]),
                      mul(compiled,stored[bit]),lo+shift,size,mu*prob)
    visit('',m.root(),I,I,0,1<<codec.B,F(1))
    # Memoization below caches only already-evaluated pure oracle calls;
    # the shipped encoder/decoder/transition logic is exercised unchanged.
    codec.probability=lambda z:cache[z]
    words={r['w']:r for r in rows};ids={r['lo'] for r in rows}
    assert len(ids)==len(rows)
    transitions=0;identities=0
    for r in rows:
        i,N,z=codec.encode(r['w'])
        assert (i,N,z)==(r['lo'],r['N'],r['z'])
        assert codec.decode(i)[0]==r['w'] and codec.output(i)==r['output']
        identities+=1
        for bit in (0,1):
            raw=str(bit)+r['w'];dest=raw
            for j in range(len(raw)+1):
                if raw[:j] in words and words[raw[:j]]['leaf']:
                    dest=raw[:j];break
            target=words[dest]
            assert codec.advance(i,bit)==target['lo']
            assert codec.output(target['lo'])==target['output']
            if r['leaf']:assert target['leaf']
            xs=H.J if r['leaf'] else (F(1,2),)
            for x in xs:
                y=point(ms[bit],point(r['true'],x))
                assert point(target['true'],H.J[0])<=y<=point(target['true'],H.J[1])
                assert abs(F(target['output'],m.Q)-out(ms,y))<=eps
            transitions+=1
    invalid=0
    for i in range(1<<codec.B):
        try:
            w,z=codec.decode(i);assert i in ids and words[w]['lo']==i
        except ValueError:
            assert i not in ids;invalid+=1
    # Re-evaluate several cached entries with the production oracle.
    sample=list(cache.items())[::max(1,len(cache)//7)]
    for z,v in sample:assert realprob(z)==v
    return dict(model=name,p=p,bits=codec.B,old_bits=codec.old_B,
        states=len(rows),transitions=transitions,unused_ids_exhaustively_rejected=invalid,
        oracle_calls=len(cache),roundtrips=identities,minimum_capacity=min(r['N'] for r in rows),
        minimum_internal_capacity=min(r['N'] for r in rows if not r['leaf']),
        max_TV_over_epsilon=float(max(r['tv'] for r in rows)/eps),
        seconds=time.perf_counter()-start,
        verifier_cache_only=True,production_tree_table=False)

def graphs():
    rows=[]
    for name,ms in CASES:
        for p in (4,8,12):
            r=graph(name,ms,p);rows.append(r);print(json.dumps(r),flush=True)
            save('GRAPHS.json',dict(passed=len(rows)==6,rows=rows,
                 total_states=sum(r['states'] for r in rows),
                 total_transitions=sum(r['transitions'] for r in rows),
                 total_unused_ids=sum(r['unused_ids_exhaustively_rejected'] for r in rows)))

def streams():
    rows=[]
    for name,ms in CASES:
        for p in (16,24):
            H=certify(ms,J=(F(1,4),F(3,4)));m=Model.from_rationals(p,H);codec=ReducedCodec(m)
            rng=random.Random(1700+p);x=F(1,2);state=0;peak=F(0);eps=F(1,1<<p)
            count=128 if (name,p) in [('S',16),('A',24)] else 64
            start=time.perf_counter()
            for t in range(count):
                b=(t&1) if t<32 else rng.randrange(2)
                state=codec.advance(state,b);x=point(ms[b],x)
                error=abs(F(codec.output(state),m.Q)-out(ms,x))
                assert error<=eps;peak=max(peak,error)
            row=dict(model=name,p=p,bits=codec.B,old_bits=codec.old_B,steps=count,
                     max_TV_over_epsilon=float(peak/eps),seconds=time.perf_counter()-start,
                     timing_includes_exact_reference=True,production_oracle_uncached=True)
            rows.append(row);print(json.dumps(row),flush=True)
            save('STREAMS.json',dict(passed=len(rows)==4,rows=rows))

def controls():
    checks=0;excluded=0
    for N in range(2,65):
        for Q in (2,4,8,16,32,64):
            for k in range(1,Q):
                q=F(k,Q);left=(2*(N-1)*k+Q)//(2*Q)
                assert left>=q*N-q-F(1,2)
                assert N-1-left>=(1-q)*N-(1-q)-F(1,2)
                if (N-1)*min(q,1-q)>F(1,2):
                    assert 0<left<N-1
                else:excluded+=1
                checks+=1
    try:ReducedCodec.split(3,3,4);raise AssertionError('boundary should be excluded')
    except AssertionError as e:
        assert str(e)=='certified capacity violated'
    ms=CASES[0][1];H=certify(ms,J=(F(1,4),F(3,4)));m=Model.from_rationals(16,H)
    x=F(3,1<<4096);lower,q=positive_power_lower(m,x,F(3,4))
    tolerance=x**2/4096
    assert lower>0 and lower**4<=x**3<=(lower+tolerance)**4
    # Compiling and using the wrong smaller margin must not be called certified.
    c=ReducedCodec(m);cc=c.capacity_certificate
    factor=F(cc['required_root_factor'])
    assert (1<<cc['offset_bits'])>factor
    assert (1<<(cc['offset_bits']-1))<=factor
    save('CONTROLS.json',dict(passed=True,exact_split_checks=checks,
        boundary_cases_outside_sufficient_condition=excluded,
        empty_child_boundary_rejected=True,tiny_base_numerator=3,tiny_base_denominator_bits=4097,
        tiny_base_power_precision=q,smaller_uncertified_offset_rejected=True))

if __name__=='__main__':
    import sys
    mode=sys.argv[1] if len(sys.argv)>1 else 'all'
    if mode in ('all','controls'):controls()
    if mode in ('all','graphs'):graphs()
    if mode in ('all','streams'):streams()
