"""Focused reproducible mathematical checks for the v1.2 patch audit.

These checks corroborate explicit counterexamples and finite instances;
they are not proofs of the general theorems. Python standard library only.
"""
from fractions import Fraction as F
from decimal import Decimal, localcontext
from itertools import product
from pathlib import Path
import hashlib
import json

def dec(x):
    if isinstance(x,F): return Decimal(x.numerator)/Decimal(x.denominator)
    return Decimal(x)

def entropy_bits(probs):
    return -sum(dec(p)*dec(p).ln() for p in probs if p)/Decimal(2).ln()

def tv(p,q):
    return sum(abs(a-b) for a,b in zip(p,q))/2

def mm(a,b):
    return [[sum(a[i][r]*b[r][j] for r in range(len(b)))
             for j in range(len(b[0]))] for i in range(len(a))]

def transpose(a): return [list(x) for x in zip(*a)]

def run():
    out={}
    delta=F(1,1000);t=F(1,10**6)
    wa=(1-delta)*t+delta/3;wb=delta/3
    omega=wa/(wa+wb)
    p=[1-t,t,F(0)]
    qa=[(1-delta)*(1-t)+delta/3,(1-delta)*t+delta/3,delta/3]
    qb=[(1-delta)/2+delta/3,delta/3,(1-delta)/2+delta/3]
    q=[omega*a+(1-omega)*b for a,b in zip(qa,qb)]
    assert sum(p)==sum(q)==1
    kl=sum(dec(a)*(dec(a)/dec(b)).ln() for a,b in zip(p,q) if a)
    rejected_bound=dec(delta)*(Decimal(3)/dec(delta)).ln()
    assert kl>Decimal('0.2877') and rejected_bound<Decimal('0.0081')
    visits=F(1,2)/(1-t)+1
    assert visits<2
    out['conditional_KL_counterexample']={
        'p_equals_q':True,'delta':str(delta),'prefix_mass':str(t/2),
        'source_next':list(map(str,p)),'smoothed_next':list(map(str,q)),
        'KL_nats':str(kl),'rejected_upper_bound':str(rejected_bound),
        'original_expected_visits':str(visits)}

    a=[[F(0),F(10)],[F(0),F(0)]]
    assert mm(a,a)==[[0,0],[0,0]]
    aa=mm(a,transpose(a))
    assert aa==[[100,0],[0,0]]
    out['Stein_resolvent_counterexample']={
        'spectral_radius':0,'inverse_radius_margin':1,
        'resolvent_on_identity_operator_norm':101}

    R=F(2);w=F(3,2)
    assert all(abs(x/4+c)<w for x in [-R,R] for c in [F(-1,2),F(1,2)])
    assert F(1,4)+F(1,4)==F(1,2)
    out['transfer_counterexample']={
        'R':str(R),'w':str(w),'sigma':1,'L_of_one':'1/2',
        'V_n':'2**(-n)','positive_limit':False}

    ratio_squared=F(4,3*(6-3)+1)
    assert ratio_squared<1
    out['QR_constant_counterexample']={
        'k':3,'N':6,'claimed_ratio_squared':str(ratio_squared)}

    alpha=[[F(1),F(0)]];stop=[[F(1,2)],[F(1,2)]]
    ta=[[F(1,4),F(1,8)],[F(0),F(1,4)]]
    tb=[[F(0),F(1,8)],[F(1,4),F(0)]]
    p_ab=mm(mm(mm(alpha,ta),tb),stop)[0][0]
    p_ba=mm(mm(mm(alpha,tb),ta),stop)[0][0]
    mu_a=sum(mm(alpha,ta)[0])
    ell=transpose(mm(alpha,ta))
    ell=[[v[0]/mu_a] for v in ell]
    correct=mm(mm(transpose(ell),tb),stop)[0][0]
    assert p_ab!=p_ba and correct==p_ab/mu_a
    out['dual_orientation']={
        'p_ab':str(p_ab),'p_ba':str(p_ba),
        'correct_future_at_b_after_a':str(correct)}

    entropies=[]
    for m in [2,10,100]:
        probs=[1-F(1,m*m),F(1,2*m*m),F(1,2*m*m)]
        variance=sum(pv*y*y for pv,y in zip(probs,[0,m,-m]))
        assert variance==1
        entropies.append({'m':m,'variance':str(variance),
                          'entropy_bits':str(entropy_bits(probs))})
    assert entropy_bits([F(1,2),F(1,2)])==1
    assert Decimal(entropies[-1]['entropy_bits'])<Decimal('0.002')
    out['equal_variance_different_entropies']=entropies
    # Exact inequalities establish rounding of both rotated points to (1,1):
    assert F(1,2)**2<F(1,2)<F(3,2)**2
    assert F(1,2)**2<F(2)<F(3,2)**2
    out['orthogonal_quantization_collision']={
        'inputs':[[1,0],[2,0]],'rotation_degrees':45,
        'both_unit_grid_outputs':[1,1]}

    distributions=[tuple(F(v,5) for v in (i,j,5-i-j))
                   for i in range(6) for j in range(6-i)]
    checked=0
    for p0,q0 in product(distributions,repeat=2):
        eps=tv(p0,q0)
        for mask in range(1,8):
            idx=[i for i in range(3) if mask>>i&1]
            s=sum(p0[i] for i in idx);r=sum(q0[i] for i in idx)
            if s>eps:
                assert r>0
                pc=[p0[i]/s for i in idx];qc=[q0[i]/r for i in idx]
                assert tv(pc,qc)<=eps/s
                checked+=1
    out['conditional_TV_exact_grid']={
        'denominator':5,'alphabet_size':3,'verified_cases':checked}
    out['status']='all checks passed'
    return out

if __name__=='__main__':
    with localcontext() as ctx:
        ctx.prec=80
        result=run()
    result['script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    target=Path(__file__).with_name('patch_verification_results.json')
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
