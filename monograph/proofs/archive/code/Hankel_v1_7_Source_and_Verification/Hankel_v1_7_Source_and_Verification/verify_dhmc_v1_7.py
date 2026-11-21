"""Exact finite checks supporting v1.7; not a proof of its general claims.

Python 3.12 standard library. Run from any directory. Results are written
beside this script. No network, model API, or benchmark data is used.
"""
from pathlib import Path
from fractions import Fraction as F
from collections import deque, defaultdict
import importlib.util
import hashlib
import itertools
import json
import random

ROOT = Path(__file__).resolve().parent
result = {'version': '1.7', 'scope': 'Exact finite checks, not proof-assistant verification, semantic validation or an LLM experiment.', 'checks': []}


def record(name, **kw):
    result['checks'].append({'name': name, 'passed': True, **kw})


def words(alphabet, depth):
    for n in range(depth + 1):
        yield from itertools.product(alphabet, repeat=n)


def distinguish(states, p, q, alphabet):
    queue = deque([(p, q, ())]); seen = {(p, q)}
    while queue:
        u, v, w = queue.popleft()
        if states[u]['output'] != states[v]['output']:
            return w
        for a in alphabet:
            pair = states[u]['transitions'][a], states[v]['transitions'][a]
            if pair not in seen:
                seen.add(pair); queue.append((*pair, w + (a,)))
    return None


def check_quotients():
    path = ROOT / 'implementation_snapshot/compiler.py'
    spec = importlib.util.spec_from_file_location('frozen_dhmc', path)
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    rng = random.Random(1701); pairs = 0; witnesses = 0; cell_checks = 0
    for n in range(2, 11):
        for _ in range(12):
            names = [f's{i}' for i in range(n)]; alphabet = ['a', 'b', 'r']
            states = {s: {'output': str(rng.randrange(3)),
                         'transitions': {a: rng.choice(names) for a in alphabet}} for s in names}
            machine = {'id': 'synthetic', 'title': 'Synthetic', 'scope': 'test',
                       'trigger': 'test', 'unknown_action': 'suspend',
                       'events': {a:a for a in alphabet}, 'states': states, 'initial': names[0]}
            quotient, witness = mod.minimize(machine)
            mapping = witness['state_map']; reachable = sorted(mapping)
            for p, q in itertools.combinations(reachable, 2):
                w = distinguish(states, p, q, alphabet)
                assert (w is None) == (mapping[p] == mapping[q])
                if w is not None:
                    assert len(w) <= len(reachable) - 1
                    witnesses += 1
                pairs += 1
            cells = defaultdict(set)
            for q, z in mapping.items(): cells[z].add(q)
            for z, cell in cells.items():
                for q in cell:
                    assert states[q]['output'] == quotient['states'][z]['output']
                    for a in alphabet:
                        assert states[q]['transitions'][a] in cells[quotient['states'][z]['transitions'][a]]
                        cell_checks += 1
            # Full-state residual reconstruction, independent of the quotient.
            fibres = defaultdict(list)
            for q in names: fibres[states[q]['output']].append(q)
            enc = {q: (c,j) for c,qs in fibres.items() for j,q in enumerate(qs)}
            dec = {v:q for q,v in enc.items()}
            for q in names:
                assert dec[enc[q]] == q
                for a in alphabet:
                    assert dec[enc[states[dec[enc[q]]]['transitions'][a]]] == states[q]['transitions'][a]
    record('quotient_witnesses_invariant_cells_and_residual', machines=108,
           state_pairs=pairs, distinguishing_suffixes=witnesses, transition_inclusions=cell_checks)


def check_guards_and_persistence():
    radius, delta, numerical, true = F(4,256),F(2,256),F(5,256),F(3,256)
    assert abs(numerical-true) <= delta
    assert abs(numerical)>radius and not abs(numerical)>radius+delta
    assert true < radius
    R=F(1,4); lipschitz=2*R; centre=R*R
    assert centre-lipschitz*R < 0
    for j in range(-20,21):
        x=R*F(j,20); assert R*R+x*x>0
    # Positive intrinsic margin: every admissible error in this grid passes.
    refinements=0
    for d in [F(1,2**k) for k in range(2,10)]:
        for j in range(-8,9):
            val=F(3,4)+F(j,8)*d
            if 2*d<F(1,2): assert abs(val)>F(1,4)+d
            refinements+=1
    for n in range(30):
        x=F(1)
        for _ in range(n):x=x/(1+x)
        assert x==F(1,n+1)
    p0=((F(1),F(0)),(F(0),F(0)))
    p1=((F(0),F(0)),(F(0),F(1)))
    def matvec(a,x):return tuple(sum(u*v for u,v in zip(row,x)) for row in a)
    for e in [(F(1),F(0)),(F(0),F(1))]:
        assert matvec(p1,matvec(p0,e))==(0,0)
        assert matvec(p0,matvec(p1,e))==(0,0)
    assert matvec(p0,(F(1),F(0)))==(1,0)
    assert matvec(p1,(F(0),F(1)))==(0,1)
    # Chain merging: outputs 0,1,2,3 has radius 3/2, not one.
    assert F(3,2)>1
    record('guards_refinement_and_persistence_counterexamples', refinements=refinements, positive_quadratic_points=41)


def check_hypotheses():
    trans={0:{'a':1,'b':2},1:{'a':1,'b':0},2:{'a':0,'b':2}}
    checked=0
    # An ambiguous observation permits both letters; concrete truth is enumerated.
    for w in words('ab',7):
        truth=0; hyp={0}
        for j,a in enumerate(w):
            allowed={'a','b'} if j%2==0 else {a}
            hyp={trans[q][b] for q in hyp for b in allowed}
            truth=trans[truth][a]
            assert truth in hyp;checked+=1
    record('sound_hypothesis_propagation', concrete_transitions=checked)


def check_base_replay():
    ratios={'a':F(1,2),'b':F(1,3),'r':F(1,4)}
    def mode(m,a):return 0 if a=='r' else (1 if a=='b' else 1-m)
    def step(x,m,a):
        r=ratios[a]
        return r*x+((1-r) if (m+(a=='b'))%2 else 0)
    checks=0; truncations=0; resets=0
    for theta in [F(1,8),F(1,16),F(1,32)]:
        def walk(depth,x,m,v,base,modes):
            nonlocal checks,truncations,resets
            if depth==8:return
            for a in ratios:
                nx=step(x,m,a);nm=mode(m,a)
                w=(a,)+v;prod=F(1);keep=[]
                for b in w:
                    keep.append(b);prod*=ratios[b]
                    if prod<=theta:break
                nv=tuple(keep);drop=w[len(nv):];nb=base
                for b in reversed(drop):nb=mode(nb,b)
                newmodes=modes+(nm,)
                assert nb==newmodes[len(newmodes)-1-len(nv)]
                rx=F(1,3);rm=nb
                for b in reversed(nv):rx,rm=step(rx,rm,b),mode(rm,b)
                assert rm==nm and abs(rx-nx)<=prod
                if len(nv)==depth+1:assert rx==nx
                else:
                    assert prod<=theta and abs(rx-nx)<=theta
                    truncations+=1
                checks+=1;resets+=a=='r'
                walk(depth+1,nx,nm,nv,nb,newmodes)
        walk(0,F(1,3),0,(),0,(0,))
    record('base_mode_replay_with_reset_and_unequal_letter_ratios', prefixes=checks,
           truncated_prefixes=truncations, reset_steps=resets, max_depth=8)


def distribution(alpha, kernels):
    law=defaultdict(F)
    def visit(i,w,mass):
        law[w]+=mass*kernels[i]['stop']
        for (a,j),v in kernels[i]['edges'].items():
            if v:visit(j,w+(a,),mass*v)
    for i,p in enumerate(alpha):visit(i,(),p)
    return dict(law)


def tv(p,q):
    return sum((abs(p.get(k,F(0))-q.get(k,F(0))) for k in p.keys()|q.keys()),F(0))/2


def check_rounding_and_flow():
    rng=random.Random(1702);checked=0
    for _ in range(30):
        n=4;raw=[rng.randrange(1,10) for _ in range(n)];alpha=[F(v,sum(raw)) for v in raw]
        kernels=[]
        for i in range(n):
            edges={(a,j):F(rng.randrange(1,5),100) for a in 'ab' for j in range(i+1,n)}
            kernels.append({'edges':edges,'stop':1-sum(edges.values(),F(0))})
        source=distribution(alpha,kernels);assert sum(source.values())==1
        visits=list(alpha)
        for i in range(n):
            for (a,j),v in kernels[i]['edges'].items():visits[j]+=visits[i]*v
        s=min(k['stop'] for k in kernels);assert sum(visits)<=1/s
        for bits in [3,5,8]:
            h=F(1,2**bits)
            ra=[(v//h)*h for v in alpha[:-1]];ra.append(1-sum(ra))
            rounded=[];local=[]
            for k in kernels:
                edges={aj:(v//h)*h for aj,v in k['edges'].items()}
                rounded.append({'edges':edges,'stop':1-sum(edges.values(),F(0))})
                local.append(sum(k['edges'].values())-sum(edges.values()))
            actual=tv(source,distribution(ra,rounded))
            init=sum(abs(u-v) for u,v in zip(alpha,ra))/2
            bound=init+sum(v*d for v,d in zip(visits,local))
            assert actual<=bound<= (n-1)*h+2*n*h/s
            # Occupation closure uses true continuation laws and rounded flows.
            f=[]
            for i in range(n):
                start=[F(0)]*n;start[i]=F(1);f.append(distribution(start,kernels))
            rv=list(ra)
            for i in range(n):
                for (a,j),v in rounded[i]['edges'].items():rv[j]+=rv[i]*v
            defect=F(0)
            for i in range(n):
                row=defaultdict(F,{w:rv[i]*v for w,v in f[i].items()})
                row[()]-=rv[i]*rounded[i]['stop']
                for (a,j),v in rounded[i]['edges'].items():
                    for w,mass in f[j].items():row[(a,)+w]-=rv[i]*v*mass
                defect+=sum(abs(x) for x in row.values())
            assert actual<=init+defect/2
            checked+=1
    record('dyadic_coupling_and_occupation_flow', exact_finite_word_models=checked)


def check_conditioning():
    distributions=[]
    for a in range(6):
        for b in range(6-a):distributions.append({0:F(a,5),1:F(b,5),2:F(5-a-b,5)})
    checked=0
    for p,q in itertools.product(distributions,repeat=2):
        eps=tv(p,q)
        for event in [{0},{1},{2},{0,1},{0,2},{1,2}]:
            mu=sum(p[k] for k in event);xi=sum(q[k] for k in event)
            if mu>eps:
                assert xi>0
                assert tv({k:p[k]/mu for k in event},{k:q[k]/xi for k in event})<=eps/mu
                checked+=1
    record('rare_context_conditional_tv', exact_cases=checked)


def check_endpoint_hulls():
    fs=[lambda t:(1+t)/(2+t),lambda t:(2-t)/3]
    lo=hi=F(1,2);points={lo};checked=0
    for depth in range(1,11):
        vals=[F(1,2)]+[f(v) for f in fs for v in (lo,hi)]
        lo,hi=min(vals),max(vals)
        points|={f(v) for f in fs for v in points}
        assert (lo,hi)==(min(points),max(points));checked+=len(points)
    record('mixed_orientation_endpoint_hulls', depths=10, enumerated_prefix_states=checked)


def check_norm():
    exact=F(81,49)+F(450,169)+F(625,529)
    triangle=F(9,7)+F(25,23)
    assert exact<=triangle**2 and triangle/2<F(9,7)
    partial=sum(F(n+1)*2**n*(F(1,3)**n+F(1,5)**n)**2 for n in range(60))
    assert partial<exact and exact-partial<F(1,10**35)
    discount=F(1,3); alphabet_size=2
    target=1/(1-alphabet_size*discount**2)
    for horizon in range(20):
        total=sum(F(alphabet_size**n)*discount**(2*n) for n in range(horizon+1))
        assert target-total==target*(alphabet_size*discount**2)**(horizon+1)
    record('corrected_rank_four_norm', exact_squared_norm=str(exact), numerical_norm=float(exact)**0.5,
           half_triangle_bound=str(triangle/2), discounted_rank_one_checks=20)


def main():
    for fn in [check_quotients,check_guards_and_persistence,check_hypotheses,
               check_base_replay,check_rounding_and_flow,check_conditioning,
               check_endpoint_hulls,check_norm]:fn()
    result['status']='all checks passed'
    result['script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result['compiler_sha256']=hashlib.sha256((ROOT/'implementation_snapshot/compiler.py').read_bytes()).hexdigest()
    (ROOT/'dhmc_verification_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
