#!/usr/bin/env python3
"""Exact recognizer and arithmetic replay for DRAGA AAK v0.6.0.

The universal proofs are written in the companion note, not formalized.
No optimizer, floating-point SVD, or finite-context approximation is used.
"""
from dataclasses import dataclass
from collections import deque
from pathlib import Path
import argparse
import itertools
import json
import platform
import sympy as S


def require(ok, message):
    if ok != True:
        raise AssertionError(message)


@dataclass
class WFA:
    A: tuple
    alpha: S.Matrix
    beta: S.Matrix

    def eval(self, word):
        row = self.alpha.T
        for letter in word:
            row = row*self.A[letter]
        return S.expand((row*self.beta)[0])


def add(F, G, scale=1):
    require(len(F.A) == len(G.A), 'Matching alphabets')
    return WFA(tuple(S.diag(a,b) for a,b in zip(F.A,G.A)),
               F.alpha.col_join(scale*G.alpha), F.beta.col_join(G.beta))


def reduced(scale=S.Rational(33,736)):
    return WFA((S.diag(S.Rational(1,2),S.Rational(1,4)),
                S.diag(S.Rational(1,4),S.Rational(1,2))),
               S.Matrix([scale,scale]), S.ones(2,1))


def pure_words(m, weights=(1,1)):
    d=len(weights); n=d*(m-1)+2
    A=[S.zeros(n) for _ in weights]
    for i,weight in enumerate(weights):
        chain=[0]+[1+i*(m-1)+j for j in range(m-1)]+[n-1]
        for j,(p,q) in enumerate(zip(chain,chain[1:])):
            A[i][p,q]=weight if j == 0 else 1
    return WFA(tuple(A),S.eye(n)[:,0],S.eye(n)[:,-1])


def monomial(word, coefficient=1, alphabet_size=2):
    n=len(word)+1; A=[S.zeros(n) for _ in range(alphabet_size)]
    for i,a in enumerate(word):
        A[a][i,i+1]=1
    return WFA(tuple(A),coefficient*S.eye(n)[:,0],S.eye(n)[:,-1])


def word_string(w):
    return ''.join(chr(97+i) for i in w) or '<empty>'


def independent_rows(pairs, n):
    basis=S.zeros(0,n); selected=[]
    for row,word in pairs:
        trial=basis.col_join(row)
        if trial.rank() > basis.rows:
            basis=trial; selected.append((row,word))
    return selected


def homogeneous_certificate(F,m):
    """All-word exact verification, via degree spaces then invariant closure."""
    n=F.beta.rows
    current=independent_rows([(F.alpha.T,())],n)
    for degree in range(m+1):
        if degree != m:
            for row,w in current:
                value=(row*F.beta)[0]
                if value != 0:
                    return {'accepted':False,'witness':word_string(w),'value':str(value)}
        current=independent_rows([(row*A,w+(i,)) for row,w in current
                                 for i,A in enumerate(F.A)],n)
    # Current contains a basis for all forward states at length m+1.
    queue=deque(current); closure=S.zeros(0,n)
    while queue:
        row,w=queue.popleft(); value=(row*F.beta)[0]
        if value != 0:
            return {'accepted':False,'witness':word_string(w),'value':str(value)}
        trial=closure.col_join(row)
        if trial.rank() > closure.rows:
            closure=trial
            queue.extend((row*A,w+(i,)) for i,A in enumerate(F.A))
    return {'accepted':True,'tail_invariant_dimension':closure.rows}


def coefficient_norm_squared(F,m):
    Q=F.beta*F.beta.T
    for _ in range(m):
        Q=sum((A*Q*A.T for A in F.A),S.zeros(Q.rows))
    return S.factor((F.alpha.T*Q*F.alpha)[0])


def equivalent(F,G):
    D=add(F,G,-1); n=D.beta.rows
    queue=deque([(D.alpha.T,())]); basis=S.zeros(0,n)
    while queue:
        row,w=queue.popleft()
        if (row*D.beta)[0] != 0:
            return False
        trial=basis.col_join(row)
        if trial.rank() > basis.rows:
            basis=trial
            queue.extend((row*A,w+(i,)) for i,A in enumerate(D.A))
    return True


def recover(F,m):
    """Recognize C_{2,m} from a complete exact bounded WFA, and construct G.

    Boundedness is an input promise. Once membership passes, boundedness of
    the reconstructed G follows from H_G=H_F-H_h with h finite support.
    """
    require(m >= 4,'The two-state reconstruction requires m >= 4')
    short=[()]+[(i,) for i in range(len(F.A))]
    K=S.Matrix([[F.eval(p+q) for q in short] for p in short])
    r=K.rank()
    if r > 2:
        return {'accepted':False,'reason':'short Hankel rank exceeds two'}
    if r == 0:
        G=WFA(tuple(S.zeros(1) for _ in F.A),S.zeros(1,1),S.ones(1,1))
    else:
        selected=None
        for rows in itertools.combinations(range(len(short)),r):
            for cols in itertools.combinations(range(len(short)),r):
                M=K.extract(rows,cols)
                if M.det() != 0:
                    selected=(rows,cols,M); break
            if selected is not None:
                break
        rows,cols,M=selected
        P=[short[i] for i in rows]; Q=[short[j] for j in cols]
        A=tuple(M.inv()*S.Matrix([[F.eval(p+(i,)+q) for q in Q] for p in P])
                for i in range(len(F.A)))
        G=WFA(A,S.Matrix([F.eval(q) for q in Q]),
              M.inv()*S.Matrix([F.eval(p) for p in P]))
    residual=add(F,G,-1)
    membership=homogeneous_certificate(residual,m)
    if not membership['accepted']:
        return {'accepted':False,'reason':'residual is not homogeneous',
                'membership':membership}
    a2=coefficient_norm_squared(residual,m)
    return {'accepted':True,'G':G,'rank':r,'a_squared':a2,'membership':membership,
            'relative_gap_bound':S.Rational(8,m+5)}


def algebra_audit():
    m,r,a,t=S.symbols('m r a t')
    require(S.expand((m+1)*(a-t)-2*r*(a+t)
                     -((m+1-2*r)*a-(m+1+2*r)*t)) == 0,
            'Uniform lower-bound rearrangement')
    require(S.simplify(a-a*(m+1-2*r)/(m+1+2*r)
                       -4*r*a/(m+1+2*r)) == 0,'Gap formula')
    # Rank-two diagonal infinite series, exact all-word kernel Gram matrix.
    K=S.Matrix([[S.Rational(16,11),S.Rational(4,3)],
                [S.Rational(4,3),S.Rational(16,11)]])
    scale=S.Rational(33,736)
    require(sorted((scale*K).eigenvals()) == [S.Rational(1,184),S.Rational(1,8)],
            'Exact norm of the retained infinite series')
    low=S.Rational(1414213,1000000); high=S.Rational(707107,500000)
    require(low**2 < 2 < high**2,'Rational isolating interval for sqrt(2)')
    m0=397
    L=S.Rational(99,100)*S.sqrt(2)-S.Rational(1,1600)
    Llow=S.Rational(99,100)*low-S.Rational(1,1600)
    gapupper=high/100+S.Rational(1,1600)
    require(Llow > S.Rational(1399445,1000000),'Certified decimal lower bound')
    require(Llow > S.Rational(9,8),'Strict structured gap beyond sigma3')
    require(2*136**2 > 191**2,'Strict gap for every m >= 18')
    require(gapupper < S.Rational(14768,1000000),'Certified additive gap')
    require(S.simplify(((m0+1)*S.sqrt(2)-2*(S.sqrt(2)+S.Rational(1,8)))
                      / (m0+3)-L) == 0,'Refined lower certificate')
    return {'m':m0,'scale':str(scale),'target_rank':2*m0+2,'reduced_rank':2,
            'retained_norm':'1/8','residual_norm':'sqrt(2)',
            'lower_uniform':str(S.sqrt(2)*S.Rational(m0-3,m0+5)),
            'lower_refined':str(L),'sigma3_upper':'9/8',
            'additive_gap_upper':'14768/1000000','strict_gap':'PASS',
            'strict_family_all_lengths_m_at_least':18,
            'note':'Target rank uses the written independent-range proof; no 796x796 floating SVD.'}


def replay_recognition():
    cases=[]
    for m,weights,scale in [(4,(1,1),S.Rational(33,736)),
                             (5,(2,-1),S.Rational(-2,7))]:
        F=add(pure_words(m,weights),reduced(scale))
        rec=recover(F,m)
        require(rec['accepted'],'Recognition must pass')
        require(equivalent(rec['G'],reduced(scale)),'Exact recovered behavior')
        require(rec['a_squared'] == sum(S.Integer(x)**2 for x in weights),
                'Coefficient norm computed without word enumeration')
        require(F.A[0]*F.A[1] == F.A[1]*F.A[0],'Commutativity')
        cases.append({'m':m,'weights':list(weights),'scale':str(scale),
                      'a_squared':str(rec['a_squared']),'status':'PASS',
                      'membership':rec['membership']})
    # General class is not limited to pure words or commuting transitions.
    residual=add(monomial((0,1,0,1),2),monomial((1,0,0,1),-3))
    F=add(residual,reduced())
    rec=recover(F,4)
    require(rec['accepted'] and rec['a_squared']==13,'Noncommutative signed fixture')
    require(equivalent(rec['G'],reduced()),'Noncommutative recovery')
    # Rank-one and rank-zero branches.
    rank1=WFA((S.Matrix([[S.Rational(1,3)]]),S.Matrix([[S.Rational(1,4)]])),
              S.Matrix([2]),S.ones(1,1))
    for G,rank in [(rank1,1),(WFA((S.zeros(1),S.zeros(1)),S.zeros(1,1),S.ones(1,1)),0)]:
        rec=recover(add(monomial((0,0,0,0)),G),4)
        require(rec['accepted'] and rec['rank']==rank and equivalent(rec['G'],G),
                'Lower-rank recovery branch')
    base=add(pure_words(4),reduced())
    bad=add(base,monomial((1,)*5,S.Rational(1,1000)))
    require(all(base.eval(w)==bad.eval(w)
                for l in range(4) for w in itertools.product(range(2),repeat=l)),
            'Negative control is invisible to all reconstruction moments')
    rejected=recover(bad,4)
    require(not rejected['accepted'],'Global membership must reject a late perturbation')
    require(rejected['membership']['witness']=='bbbbb','Exact late witness')
    return cases,rejected


def replay_homogeneous_blocks():
    reports=[]
    for m in (2,3,4):
        values={w:S.Integer((7*i+3)%11-5)
                for i,w in enumerate(itertools.product(range(2),repeat=m))}
        a2=sum(x*x for x in values.values())
        for p in range(m+1):
            U=list(itertools.product(range(2),repeat=p))
            V=list(itertools.product(range(2),repeat=m-p))
            C=S.Matrix([[values[u+v] for v in V] for u in U])
            require(sum(x*x for x in C)==a2,'All split blocks have the same Frobenius norm')
        # Trace multiplicity against a bounded geometric rank-one Hankel.
        lam=[S.Rational(1,4),S.Rational(1,5)]
        total=0
        for p in range(m+1):
            for u in itertools.product(range(2),repeat=p):
                for v in itertools.product(range(2),repeat=m-p):
                    total += values[u+v]*S.prod(lam[i] for i in u+v)
        expected=(m+1)*sum(c*S.prod(lam[i] for i in w) for w,c in values.items())
        require(S.expand(total-expected)==0,'Hankel trace counts m+1 splittings')
        reports.append({'m':m,'coefficient_norm_squared':str(a2),'status':'PASS'})
    return reports


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=Path('certificates.json'))
    args=parser.parse_args()
    report={'version':'0.6.0','python':platform.python_version(),'sympy':S.__version__,
            'arithmetic':'EXACT','universal_proof':'WRITTEN_NOT_FORMALIZED',
            'independent_review':'NOT_DONE','novelty':'NOT_ESTABLISHED',
            'exact_optimum':'NOT_COMPUTED','global_binary_constant':'UNRESOLVED'}
    report['example']=algebra_audit()
    print('PASS: symbolic gap and exact strict-gap example',flush=True)
    report['recognition'],report['negative_control']=replay_recognition()
    print('PASS: full-WFA recognition, recovery, rank branches and late-word rejection',flush=True)
    report['homogeneous_blocks']=replay_homogeneous_blocks()
    print('PASS: homogeneous split identities and trace multiplicity',flush=True)
    args.output.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print('Written:',args.output)


if __name__=='__main__':
    main()
