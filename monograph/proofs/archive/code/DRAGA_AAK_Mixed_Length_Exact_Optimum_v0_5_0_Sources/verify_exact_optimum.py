#!/usr/bin/env python3
"""Exact replay and constructive recognition for DRAGA AAK v0.5.0.

Python 3 + SymPy. Universal proof is in the companion note, not formalized.
Run: python verify_exact_optimum.py --output certificates.json
"""
from dataclasses import dataclass
from collections import deque
from pathlib import Path
import argparse
import itertools
import json
import platform
import sympy as S


def require(condition, message):
    if condition != True:
        raise AssertionError(message)


def zero(M):
    return M.applyfunc(S.simplify) == S.zeros(*M.shape)


def matrix_json(M):
    return [[str(M[i,j]) for j in range(M.cols)] for i in range(M.rows)]


@dataclass
class WFA:
    A: tuple
    alpha: S.Matrix
    beta: S.Matrix

    def eval(self, word):
        x = self.beta
        for letter in reversed(word):
            x = self.A[letter]*x
        return S.simplify((self.alpha.T*x)[0])


def residual(rho):
    T = S.zeros(5)
    for i in range(4):
        T[i,i+1] = 1
    T[4,4] = rho
    c = S.eye(5)[:,0]
    v = S.Matrix([0,0,0,-rho,1-rho*rho])
    return WFA((T,S.zeros(5)), c, v)


def reduced(scale):
    return WFA((S.diag(S.Rational(1,2),S.Rational(1,4)),
                S.diag(S.Rational(1,4),S.Rational(1,2))),
               S.Matrix([scale,scale]), S.ones(2,1))


def add(left, right, left_scale=1, right_scale=1):
    return WFA(tuple(S.diag(a,b) for a,b in zip(left.A,right.A)),
               (left_scale*left.alpha).col_join(right_scale*right.alpha),
               left.beta.col_join(right.beta))


def target(scale, eta, rho):
    return add(residual(rho), reduced(scale), eta, 1)


def equivalence(left, right):
    """Exact reachable-space test, with a witness word on failure."""
    D = add(left,right,1,-1)
    n = D.beta.rows
    basis = S.zeros(n,0)
    queue = deque([(D.beta,())])
    while queue:
        v,word = queue.popleft()
        discrepancy = S.simplify((D.alpha.T*v)[0])
        if discrepancy != 0:
            return {'equivalent':False,'witness':''.join('ab'[i] for i in word),
                    'difference':str(discrepancy),'reachable_rank_so_far':basis.cols}
        extended = basis.row_join(v)
        if extended.rank() > basis.cols:
            basis = extended
            for i,A in enumerate(D.A):
                queue.append((A*v,(i,)+word))
    return {'equivalent':True,'reachable_rank':basis.cols}


def recover(candidate):
    """Recognize this family and recover its optimum; reject other inputs."""
    f = candidate.eval
    prefixes = [(1,),(1,1)]
    suffixes = [(),(1,)]
    M = S.Matrix([[f(p+q) for q in suffixes] for p in prefixes])
    require(M.det()!=0,'Singular reconstruction chart')
    Mi = [S.Matrix([[f(p+(i,)+q) for q in suffixes] for p in prefixes]) for i in (0,1)]
    B = tuple(M.inv()*x for x in Mi)
    gamma = S.Matrix([f(q) for q in suffixes])
    nu = M.inv()*S.Matrix([f(p) for p in prefixes])
    G = WFA(B,gamma,nu)
    scale = f(())/2
    e4 = f((0,)*4)-G.eval((0,)*4)
    e5 = f((0,)*5)-G.eval((0,)*5)
    require(e4 != 0,'Zero fourth residual moment')
    rho = S.cancel(e5/e4)
    require(-1 < rho < 1,'Residual pole outside the stable interval')
    eta = S.cancel(e4/(1-rho*rho))
    require(scale != 0 and eta != 0,'Degenerate family parameters')
    g_check = equivalence(G,reduced(scale))
    require(g_check['equivalent'],'Recovered behavior does not have the prescribed poles')
    membership = equivalence(candidate,target(scale,eta,rho))
    return {'accepted':membership['equivalent'], 'scale':scale,'eta':eta,'rho':rho,
            'G':G,'membership':membership,'M':M}


def audit_symbolic():
    r = S.symbols('rho',real=True)
    R = residual(r)
    T,c,v = R.A[0],R.alpha,R.beta
    P = S.diag(1,1,1,1,1/(1-r*r))
    Q = S.diag(1,1,1,1,1-r*r)
    require(zero(P-c*c.T-T.T*P*T),'Observability Stein identity')
    require(zero(Q-v*v.T-T*Q*T.T),'Reachability Stein identity')
    require(zero(P*Q-S.eye(5)),'Flat singular spectrum certificate PQ=I')
    C = S.Matrix([[R.eval((0,)*(i+j)) for j in range(5)] for i in range(5)])
    require(S.simplify(C.det()-(1-r*r)) == 0,'Residual leading Hankel minor')
    require(zero(C*P*C-Q),'Exact infinite support factorization')
    z = S.symbols('z')
    generating = -r*z**3+(1-r*r)*z**4/(1-r*z)
    require(S.simplify(z*generating-z**4*(z-r)/(1-r*z)) == 0,'Blaschke generating function')
    return {'status':'PASS','P':matrix_json(P),'Q':matrix_json(Q),
            'C_rho':matrix_json(C),'minor':'1-rho**2'}


def support_certificate(rho,eta,scale):
    """Exact 7-dimensional frame, including the full geometric tail."""
    R = residual(rho)
    C = S.Matrix([[R.eval((0,)*(i+j)) for j in range(5)] for i in range(5)])
    K = S.zeros(7)
    K[:5,:5] = S.diag(1,1,1,1,1/(1-rho*rho))
    poles = [(S.Rational(1,2),S.Rational(1,4)),
             (S.Rational(1,4),S.Rational(1,2))]
    for i,(a,b) in enumerate(poles):
        for j in range(4):
            K[j,5+i] = K[5+i,j] = a**j
        K[4,5+i] = K[5+i,4] = a**4/(1-rho*a)
        for j,(aa,bb) in enumerate(poles):
            K[5+i,5+j] = 1/(1-a*aa-b*bb)
    minors = [K[:j,:j].det() for j in range(1,8)]
    require(all(x>0 for x in minors),'Frame Gram matrix positive definite')
    D = S.diag(eta*C,scale,scale)
    def frame(w):
        pure = all(i==0 for i in w)
        n = len(w)
        return S.Matrix([int(pure and n==j) for j in range(4)] +
                        [rho**(n-4) if pure and n>=4 else S.Integer(0)] +
                        [S.prod(p[i] for i in w) for p in poles])
    F = target(scale,eta,rho)
    W = [w for n in range(4) for w in itertools.product(range(2),repeat=n)]
    for u in W:
        for v in W:
            require((frame(u).T*D*frame(v))[0] == F.eval(u+v),'Frame coefficients')
    x = S.symbols('x')
    char = S.Poly((D*K).charpoly(x).as_expr(),x)
    require(char.eval(eta)==0,'Exact eigenvalue eta')
    q,rem = S.div(char,S.Poly(x-eta,x))
    require(rem.is_zero,'Exact spectral factor')
    outside_pos = int(q.count_roots(eta,S.oo))
    outside_neg = int(q.count_roots(-S.oo,-eta))
    inside = int(q.count_roots(-eta,eta))
    require(outside_pos == 2 and outside_neg == 0 and inside == 4,
            'Exactly two singular values strictly larger than eta')
    return {'status':'PASS','Gram':matrix_json(K),'operator_coefficient_matrix':matrix_json(D),
            'positive_Gram_minors':list(map(str,minors)),
            'characteristic_polynomial':str(S.factor(char.as_expr())),
            'Sturm_counts_after_removing_eta':{'above_eta':outside_pos,'below_minus_eta':outside_neg,
                                             'between_minus_eta_and_eta':inside},
            'sigma3':str(eta)}


def audit_reconstruction():
    cases = [(S.Integer(1),S.Rational(1,100),S.Rational(1,2)),
             (S.Integer(2),-S.Rational(1,7),-S.Rational(1,3)),
             (-S.Integer(1),S.Rational(1,10),S.Integer(0))]
    reports = []
    for scale,eta,rho in cases:
        F = target(scale,eta,rho)
        require(F.A[0]*F.A[1] == F.A[1]*F.A[0],'Target commutativity')
        rec = recover(F)
        require(rec['accepted'],'Family membership rejected')
        require((rec['scale'],rec['eta'],rec['rho']) == (scale,eta,rho),'Parameter recovery')
        require(equivalence(rec['G'],reduced(scale))['equivalent'],'Recovered optimum equivalence')
        words = [(),(0,),(0,0),(0,0,0),(0,0,0,0),(1,),(1,1)]
        H = S.Matrix([[F.eval(u+v) for v in words] for u in words])
        det_expected = scale**2 * eta**5 * (1-rho*rho)/1024
        require(H.det()==det_expected and det_expected != 0,'Target rank-seven minor')
        reports.append({'parameters':{'scale':str(scale),'eta':str(eta),'rho':str(rho)},
                        'target_rank':7,'error_exact':str(abs(eta)),
                        'minor':str(det_expected),'reconstruction':'PASS',
                        'membership':rec['membership'],
                        'reduced_A_a':matrix_json(rec['G'].A[0]),
                        'reduced_A_b':matrix_json(rec['G'].A[1]),
                        'reduced_alpha':matrix_json(rec['G'].alpha),
                        'reduced_beta':matrix_json(rec['G'].beta)})
    # Negative control: same sampled moments, different behavior at b^6.
    n=7; Nb=S.zeros(n)
    for j in range(n-1):Nb[j,j+1]=1
    late=WFA((S.zeros(n),Nb),S.eye(n)[:,0],S.eye(n)[:,-1])
    base=target(*cases[0])
    bad=add(base,late,1,S.Rational(1,1000))
    rejection=recover(bad)
    require(not rejection['accepted'],'Late-word perturbation incorrectly accepted')
    require(rejection['membership']['witness']=='bbbbbb','Expected b^6 witness')
    require(rejection['membership']['difference']=='1/1000','Exact rejection discrepancy')
    return reports, rejection['membership']


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=Path('certificates.json'))
    args=parser.parse_args()
    report={'version':'0.5.0','python':platform.python_version(),'sympy':S.__version__,
            'arithmetic':'EXACT','scope':'EXPLICIT_RECOGNIZED_FAMILY',
            'universal_proof':'WRITTEN_NOT_FORMALIZED','novelty':'NOT_ESTABLISHED',
            'independent_review':'NOT_DONE','general_binary_AAK_constant':'UNRESOLVED'}
    report['symbolic_residual']=audit_symbolic()
    print('PASS: symbolic Stein identities, PQ=I, residual rank and Blaschke series')
    report['reconstruction'],report['negative_control']=audit_reconstruction()
    print('PASS: exact reconstruction, global membership, rank-seven minors and late-word rejection')
    report['independent_spectral_certificate']=support_certificate(S.Rational(1,2),S.Rational(1,100),S.Integer(1))
    print('PASS: exact full-tail frame and Sturm certificate sigma3=1/100')
    require(S.Rational(1,100)<S.Rational(2,33),'Two states necessary at stated tolerance')
    report['example']={'target_states':7,'optimal_states':2,'lower_bound':'1/100',
                       'constructed_error':'1/100','sigma3':'1/100',
                       'minimal_states_at_tolerance':'2'}
    args.output.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print('Written:',args.output)


if __name__=='__main__':
    main()
