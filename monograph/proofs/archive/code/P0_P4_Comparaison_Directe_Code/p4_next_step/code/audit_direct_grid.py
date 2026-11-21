"""Independent exact-rational checks for the direct-grid baseline and selector."""
from pathlib import Path
from fractions import Fraction as F
import json, random, sys, time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'reference_code'))
if not (ROOT/'reference_code').exists():
    sys.path.insert(0, str(ROOT.parent/'p4_cone_audit/code'))
from certificate import certify, hmm, point
from model import Model
from direct_grid import DirectGrid, HybridCodec

CASES = [
    ('S', hmm('1/4','1/4','6/25','19/25')),
    ('A', hmm('6/25','13/50','6/25','3/4')),
]

def out(ms, x):
    a,b,c,d=ms[0]
    return (a+c)*x+(b+d)*(1-x)

def save(name, data):
    (ROOT/'results'/name).write_text(json.dumps(data, indent=2)+'\n')

def audit():
    graphs=[]; streams=[]; comparisons=[]; conjugacies=0
    for name, ms in CASES:
        cert=certify(ms, J=(F(1,4), F(3,4)))
        for p in (4,8,12,16,24,64,128,256):
            model=Model.from_rationals(p, cert)
            start=time.perf_counter(); hybrid=HybridCodec(model)
            prep=time.perf_counter()-start
            grid=DirectGrid(model)
            comparisons.append(dict(model=name,p=p,grid_bits=grid.B,
                fractal_bits=hybrid.fractal_bits,selected=hybrid.mode,
                selected_bits=hybrid.B,grid_states=grid.N+1,
                preparation_seconds=prep))
            print(json.dumps(comparisons[-1]), flush=True)
            # Test the identity on 65 rational states and both branches.
            for j in range(65):
                x=F(1,4)+F(j,128); q=out(grid.matrices,x)
                for b in (0,1):
                    assert grid.continuous_update(q,b)==out(grid.matrices, point(grid.matrices[b],x))
                    conjugacies+=1
            if p<=12:
                peak=F(0)
                for i in range(grid.N+1):
                    q=grid.output(i)
                    for b in (0,1):
                        qnext=grid.continuous_update(q,b)
                        assert grid.lo<=qnext<=grid.hi
                        target=grid.advance(i,b)
                        assert 0<=target<=grid.N
                        # The local rounding condition plus the proved contraction
                        # establishes the invariant for arbitrarily long histories.
                        err=abs(grid.output(target)-qnext)
                        assert err<=grid.spacing/2
                        peak=max(peak,err)
                graphs.append(dict(model=name,p=p,states=grid.N+1,
                    transitions=2*(grid.N+1),bits=grid.B,
                    max_rounding_over_epsilon=float(peak/grid.epsilon),
                    proved_total_bound_over_epsilon=float(grid.error_total/grid.epsilon)))
            if p in (16,24):
                rng=random.Random(7100+p); i=grid.initial; x=F(1,2)
                ih=hybrid.initial; peak=abs(grid.output(i)-out(ms,x))
                for t in range(128):
                    b=t%2 if t<32 else rng.randrange(2)
                    i=grid.advance(i,b); ih=hybrid.advance(ih,b)
                    assert hybrid.mode=='grid' and i==ih
                    assert hybrid.output(ih)==grid.output(i)
                    x=point(ms[b],x)
                    err=abs(grid.output(i)-out(ms,x))
                    assert err<=grid.error_total
                    peak=max(peak,err)
                streams.append(dict(model=name,p=p,steps=128,
                    bits=grid.B,max_TV_over_epsilon=float(peak/grid.epsilon),
                    bound_over_epsilon=float(grid.error_total/grid.epsilon)))
        save('DIRECT_GRID.json',dict(passed=False,graphs=graphs,streams=streams,
             comparison=comparisons,exact_conjugacy_checks=conjugacies))
    # Integer half-way boundaries (including the asymmetric sign of alpha)
    # are tested separately, so exact ties do not rely on random coverage.
    m=Model.from_rationals(8,certify(CASES[0][1],J=(F(1,4),F(3,4))))
    g=DirectGrid(m); ties=0
    for i in range(g.N):
        q=g.lo+(F(i)+F(1,2))*g.spacing
        assert g.quantize(q)==i+1; ties+=1
    data=dict(passed=True,graphs=graphs,streams=streams,comparison=comparisons,
              exact_conjugacy_checks=conjugacies,exact_tie_checks=ties,
              total_graph_states=sum(x['states'] for x in graphs),
              total_graph_transitions=sum(x['transitions'] for x in graphs),
              scope='Full local graph checks only at p=4,8,12; 128-step source streams at p=16,24. No full graph claim at larger p. No new test of the fractal codec transitions.')
    save('DIRECT_GRID.json',data)
    print(json.dumps({k:v for k,v in data.items() if k not in ('graphs','streams','comparison')}),flush=True)

if __name__=='__main__':
    audit()
