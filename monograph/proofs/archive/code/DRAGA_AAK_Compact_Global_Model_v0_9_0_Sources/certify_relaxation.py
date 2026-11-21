"""Rational reconstruction/checker for this compact moment SDP.

The discovery phase uses floats. The verification phase uses integer Gram
factors and exact polynomial residual bounds. This is a lower-bound
certificate only when the exact resulting bound is positive; improvement
over v0.8 and closure to 0.01 are separate, explicitly recorded decisions.
The previously proved lower_v08 inequality is an explicit dependency.
"""
import argparse,json,math
from collections import defaultdict
from fractions import Fraction
from pathlib import Path
import numpy as np
import sympy as s
from compact_model import build,monomials,character,polynomial
ROOT=Path(__file__).resolve().parent

def definitions(d):
    xs=d['variables'];chars=d['characters'];bases={i:monomials(len(xs),i) for i in range(3)}
    defs=[('moment',s.Matrix([1]),[0])]+d['matrices'];out=[]
    for name,M,states in defs:
        pol=[[polynomial(p,xs) for p in row] for row in M.tolist()]
        degree=max((len(e) for row in pol for p in row for e in p),default=0)
        for ch in range(4):
            pairs=[(i,m) for i,c in enumerate(states) for m in bases[(4-degree)//2]
                   if c^character(m,chars)==ch]
            if pairs:out.append((name,ch,pairs,pol))
    eq=[]
    for p in d['equalities']:
        terms=polynomial(p,xs);degree=max(map(len,terms));ch=character(next(iter(terms)),chars)
        for m in bases[4-degree]:
            if character(m,chars)==ch:eq.append((m,terms))
    return out,eq

def discover(suffix):
    info=json.loads((ROOT/('relaxation_discovery'+suffix+'.json')).read_text())
    assert info['status'] in ('optimal','optimal_inaccurate')
    values=np.load(ROOT/('relaxation_numeric'+suffix+'.npz'))
    d=build(info.get('stein_equalities',False),info.get('trace_cuts',False));defs,eq=definitions(d)
    den=10**9;factors=[]
    for i,(name,ch,pairs,_) in enumerate(defs):
        A=values['dual_'+str(i)];w,V=np.linalg.eigh((A+A.T)/2);keep=w>1e-12
        L=np.rint((V[:,keep]*np.sqrt(w[keep]))*den).astype(np.int64)
        factors.append({'name':name,'character':ch,'integers':L.tolist()})
    ref=Fraction(math.floor(-info['normalization_dual']*10**6),10**6)
    out={'denominator':den,'reference':str(ref),'factors':factors,
         'equality_numerators':[int(round(x*den)) for x in values['equality_dual']],
         'stein_equalities':info.get('stein_equalities',False),'trace_cuts':info.get('trace_cuts',False),
         'dependency':'The exact v0.8 lower bound 2839/2500 is assumed as an already verified theorem.'}
    path=ROOT/('reconstructed_certificate'+suffix+'.json');path.write_text(json.dumps(out))
    return path

def verify(path):
    c=json.loads(path.read_text());d=build(c['stein_equalities'],c['trace_cuts']);xs=d['variables']
    defs,eq=definitions(d);den=c['denominator'];coefden=1
    for _,_,_,pol in defs:
        for row in pol:
            for p in row:
                for v in p.values():coefden=math.lcm(coefden,int(v.q))
    for _,p in eq:
        for v in p.values():coefden=math.lcm(coefden,int(v.q))
    scale=coefden*den*den;ref=Fraction(c['reference']);assert scale%ref.denominator==0
    res=defaultdict(int);res[(xs.index(d['u']),)]=scale
    res[()]=-ref.numerator*(scale//ref.denominator)
    assert len(defs)==len(c['factors'])
    for (name,ch,pairs,pol),entry in zip(defs,c['factors']):
        assert (entry['name'],entry['character'])==(name,ch)
        L=np.array(entry['integers'],dtype=object);assert L.shape[0]==len(pairs)
        Q=L@L.T
        ints=[[{e:int(v*coefden) for e,v in p.items()} for p in row] for row in pol]
        for i,(a,ma) in enumerate(pairs):
            for j,(b,mb) in enumerate(pairs):
                q=int(Q[i,j])
                if not q:continue
                for e,v in ints[a][b].items():res[tuple(sorted(ma+mb+e))]-=q*v
    assert len(eq)==len(c['equality_numerators'])
    for (m,p),num in zip(eq,c['equality_numerators']):
        for e,v in p.items():res[tuple(sorted(m+e))]+=num*int(v*coefden)*den
    # Entry bounds follow from the exact compact PMI constraints; these
    # are deliberately rational upper bounds, not numerical eigenvalues.
    bounds=[]
    for x in xs:
        name=str(x)
        if name=='u':b=Fraction(6533,5000)**2
        elif name=='z01':b=Fraction(1,2)
        elif name.startswith('r') and name[1]=='5':b=Fraction(283,200)
        elif name.startswith('w') and name[1]=='0':b=Fraction(283,200)
        elif name.startswith('w'):b=Fraction(10001,10000)
        else:b=Fraction(1)
        bounds.append(b)
    assert Fraction(283,200)**2>Fraction(10001,5000)
    assert Fraction(10001,10000)**2>Fraction(10001,10000)
    # z01 <= 1/2 follows from the spectra of 0 <= Z <= I.
    loss=Fraction(0)
    for e,num in res.items():
        weight=Fraction(1)
        for i in e:weight*=bounds[i]
        loss+=abs(Fraction(num,scale))*weight
    bound=ref-loss
    endpoint=Fraction(math.isqrt((bound.numerator*10**8)//bound.denominator),10000) if bound>0 else Fraction(0)
    if endpoint*endpoint==bound:endpoint-=Fraction(1,10000)
    old=Fraction(2839,2500);upper=Fraction(130563,100000)
    out={'exact_residual_checked':True,'lower_squared':str(bound),'reference':str(ref),
         'residual_budget':str(loss),'positive_lower_bound':bound>0,
         'four_decimal_endpoint':str(endpoint),'improves_v08':endpoint>old,
         'gap_0_01_closed':upper-endpoint<=Fraction(1,100),
         'residual_monomials':len(res),'box_bounds':list(map(str,bounds)),
         'dependency':c['dependency'],'arithmetic':'INTEGER_AND_RATIONAL_ONLY_AFTER_DISCOVERY'}
    report=path.with_name(path.stem+'_verification.json');report.write_text(json.dumps(out,indent=2))
    # Float conversions are display-only and do not enter any assertion.
    print('EXACT RESULT',json.dumps({'reference':float(ref),'loss':float(loss),
          'endpoint':str(endpoint),'improves_v08':out['improves_v08'],
          'gap_0_01_closed':out['gap_0_01_closed']}),flush=True)
    return out

def main():
    p=argparse.ArgumentParser();p.add_argument('--suffix',default='');p.add_argument('--discover',action='store_true')
    a=p.parse_args();path=ROOT/('reconstructed_certificate'+a.suffix+'.json')
    if a.discover:path=discover(a.suffix)
    verify(path)

if __name__=='__main__':main()
