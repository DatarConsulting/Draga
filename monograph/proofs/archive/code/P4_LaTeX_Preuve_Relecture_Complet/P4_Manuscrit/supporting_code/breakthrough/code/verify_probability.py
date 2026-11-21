"""Independent rational acceptance checker for the positive machines."""
import json,copy,itertools,math
from pathlib import Path
from fractions import Fraction as F
import numpy as np
import sympy as sp
import verify as v
ROOT=Path(__file__).resolve().parents[1]
def check(c):
    r=[[F(x) for x in row] for row in c['rows_exact']];K=c['alphabet_size'];D=2**c['Q'];G=2**c['belief_fractional_bits']
    z=[[F(x,D) for x in row] for row in c['coefficient_integers']]
    assert len(r)==len(z)==2 and all(len(x)==2*K+1 for x in r+z)
    assert all(sum(row)==1 and min(row)>0 for row in r+z)
    dr=max(sum(abs(a-b) for a,b in zip(x,y))/2 for x,y in zip(r,z));assert dr==F(c['coefficient_row_TV'])
    tau=min(x[-1] for x in z);tau0=min(x[-1] for x in r)
    assert tau==F(c['stopping_min']) and dr/tau0==F(c['coefficient_word_TV'])
    q=max(abs(z[0][2*a]*z[1][2*a+1]-z[0][2*a+1]*z[1][2*a])/min(sum(z[0][2*a:2*a+2]),sum(z[1][2*a:2*a+2]))**2 for a in range(K))
    assert q==F(c['contraction'])<1
    C=(sum(abs(sum(z[0][2*a:2*a+2])-sum(z[1][2*a:2*a+2])) for a in range(K))+abs(z[0][-1]-z[1][-1]))/2
    E=F(c['belief_error']);assert C==F(c['event_lipschitz']) and E>=q*E+F(1,2*G)
    delta=F(c['reduction_TV']);assert C*E/tau==F(c['machine_word_TV_bound'])
    assert delta+dr/tau0+C*E/tau==F(c['total_word_TV_bound'])<=F(1,100)
    assert c['state_bits']==c['belief_fractional_bits']+1 and c['state_bytes']==math.ceil(c['state_bits']/8)
    assert G<=2**22 and (G*D)*G<2**128
    rank=None
    if 'other_rows' in c:
        o=[[F(x) for x in row] for row in c['other_rows']]
        assert all(sum(row)==1 and min(row)>0 for row in o) and delta==F(1,1000)
        A=[]
        for a in range(K):A.append(v.block(np.array([row[2*a:2*a+2] for row in r],object),np.array([row[2*a:2*a+2] for row in o],object)))
        al=np.array([[(1-delta)/2,(1-delta)/2,delta/2,delta/2]],object);be=np.array([[row[-1]] for row in r+o],object)
        words=[()]+[w for d in range(1,4) for w in itertools.product(range(K),repeat=d)]
        rows=[];cols=[]
        for w in words:
            prod=v.eye(4)
            for a in w:prod=prod@A[a]
            rows.append((al@prod)[0]);cols.append((prod@be)[:,0])
        R=np.array(rows,object);Cmat=np.array(cols,object).T
        ri=list(sp.Matrix(R.T).rref()[1])[:4];ci=list(sp.Matrix(Cmat).rref()[1])[:4]
        assert len(ri)==len(ci)==4
        assert v.det(R[ri]) and v.det(Cmat[:,ci]);rank=4
    return {'id':c['id'],'exact_acceptance':True,'source_rank':rank,'TV_bound':float(F(c['total_word_TV_bound']))}
def run():
    d=json.loads((ROOT/'results/probability.json').read_text());records=[check(c) for c in d['synthetic']+d['trained']];bad=[]
    for key in ['row','contraction','belief_error']:
        c=copy.deepcopy(d['synthetic'][0])
        if key=='row':c['coefficient_integers'][0][0]+=1
        elif key=='contraction':c['contraction']='0'
        else:c['belief_error']='0'
        try:check(c)
        except AssertionError:bad.append(key)
    assert len(bad)==3
    (ROOT/'results/probability_verification.json').write_text(json.dumps({'all_passed':True,'records':records,'negative_controls_rejected':bad},indent=2))
    print('PASS',len(records),flush=True)
if __name__=='__main__':run()
