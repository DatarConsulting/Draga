#!/usr/bin/env python3
"""Exact finite checks for DRAGA AAK v0.4.0 (Python 3 + SymPy).

Universal proofs are written in the companion note, not formalized here.
Run: python verify_homogeneous.py --output certificates.json
"""
import argparse
import itertools
import json
import platform
from pathlib import Path
import sympy as s


def require(condition, message):
    if condition != True:
        raise AssertionError(message)


def family(d, m):
    n = d * (m - 1) + 2
    idx = {(i, j): 1 + i * (m - 1) + j - 1
           for i in range(d) for j in range(1, m)}
    A = []
    for i in range(d):
        M = s.zeros(n)
        chain = [0] + [idx[i, j] for j in range(1, m)] + [n - 1]
        for u, v in zip(chain, chain[1:]):
            M[u, v] = 1
        A.append(M)
    small = s.zeros(n)
    small[0, -1] = small[-1, 0] = s.sqrt(d)
    for i in range(d):
        for j in range(1, m):
            small[idx[i, j], idx[i, m-j]] = 1
    return A, small, idx


def audit_family():
    replay = []
    z = s.symbols('z')
    for d, m in [(1, 2), (1, 4), (2, 2), (2, 3), (2, 4), (2, 7),
                 (2, 12), (3, 3), (4, 3)]:
        A, H, idx = family(d, m)
        n = H.rows
        for i, M in enumerate(A):
            require((M**m)[0, -1] == 1, 'Pure word coefficient')
            require(M**(m+1) == s.zeros(n), 'Nilpotence')
            require(all((M**j)[0, -1] == 0 for j in range(m)), 'No shorter word')
            for j, N in enumerate(A):
                if i != j:
                    require(M*N == N*M == s.zeros(n), 'Cross products and commutativity')
        expected = s.expand((z-d)**2 * (z-1)**(d*(m-1)))
        require((H.T*H).charpoly(z).as_expr().expand() == expected, 'Squared spectrum')
        require(H.det() != 0, 'Minimal rank')
        # Full word section, independent of the reduced support formula.
        if (d, m) in [(1, 4), (2, 2), (2, 3), (2, 4), (3, 3)]:
            W = [w for j in range(m+1) for w in itertools.product(range(d), repeat=j)]
            def coeff(w):
                return int(len(w) == m and len(set(w)) == 1)
            full = s.Matrix([[coeff(u+v) for v in W] for u in W])
            J = s.zeros(len(W), n)
            J[W.index(()), 0] = 1
            for i in range(d):
                for j in range(1, m):
                    J[W.index((i,)*j), idx[i, j]] = 1
                J[W.index((i,)*m), -1] = 1/s.sqrt(d)
            require(s.simplify(J.T*J) == s.eye(n), 'Support isometry')
            require(s.simplify(J*H*J.T) == full, 'Full word Hankel support')
            # A stable rank-one series, to check the exact trace counting.
            b = [s.Rational(i+1, 10*d) for i in range(d)]
            def g(w):
                return s.prod(b[i] for i in w)
            inner = sum(full[u, v]*g(W[u]+W[v])
                        for u in range(len(W)) for v in range(len(W)))
            coeff_inner = sum(coeff(w)*g(w) for w in W if len(w) == m)
            require(inner == (m+1)*coeff_inner, 'Trace multiplicity m+1')
            Gnorm = 1/(1-sum(x*x for x in b))
            require(inner**2 <= d*Gnorm**2, 'Operator/nuclear duality: rank-one fixture')
            require(coeff_inner**2 <= d*Gnorm**2/s.Integer(m+1)**2,
                    'Normalized coefficient trace bound')
        replay.append({'d': d, 'm': m, 'rank': n, 'status': 'PASS'})
    return replay


def audit_arbitrary_homogeneous():
    # Signed, noncommutative coefficient fixtures ensure the block argument
    # does not rely on the pure-word family or commutativity.
    replay = []
    for d, m in [(2, 2), (2, 3), (2, 4), (3, 3)]:
        Wm = list(itertools.product(range(d), repeat=m))
        values = {w: s.Integer((i*7+3) % 11 - 5) for i, w in enumerate(Wm)}
        a2 = sum(x*x for x in values.values())
        require(a2 > 0, 'Nonzero homogeneous fixture')
        for p in range(m+1):
            U = list(itertools.product(range(d), repeat=p))
            V = list(itertools.product(range(d), repeat=m-p))
            C = s.Matrix([[values[u+v] for v in V] for u in U])
            require(sum(x*x for x in C) == a2, 'Frobenius norm independent of split')
            if p == 1:
                require(s.trace(C*C.T) == a2 and C.rows == d, 'Largest eigenvalue at least trace/d')
        replay.append({'d': d, 'm': m, 'coefficient_norm_squared': str(a2), 'status': 'PASS'})
    return replay


def audit_algebra():
    m, k, a, A, t = s.symbols('m k a A t')
    require(s.expand((m+1)*(a-t)-k*(A+t)-((m+1)*a-k*A-(m+1+k)*t)) == 0,
            'Trace lower-bound rearrangement')
    x = s.symbols('x', positive=True)
    require(s.limit((x-1)/(x+3), x, s.oo) == 1, 'Asymptotic sharpness')
    N = s.symbols('N', positive=True, integer=True)
    require(s.simplify(((m-1)/(m+3)).subs(m, 4*N-3)-(1-1/N)) == 0,
            'Exact table sequence')
    examples = []
    for length, ratio in [(37, s.Rational(9,10)), (397, s.Rational(99,100)),
                          (3997, s.Rational(999,1000))]:
        require(s.Rational(length-1, length+3) == ratio, 'Table coefficient')
        require(2*ratio**2 > 1, 'Table lower bound stronger than sigma3')
        examples.append({'m': length, 'target_rank': 2*length,
                         'sigma3': '1', 'ratio_lower': str(s.sqrt(2)*ratio),
                         'ratio_upper': 'sqrt(2)'})
    return examples


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('certificates.json'))
    args = parser.parse_args()
    report = {'version': '0.4.0', 'python': platform.python_version(), 'sympy': s.__version__,
              'arithmetic': 'EXACT', 'universal_proof': 'WRITTEN_NOT_FORMALIZED',
              'global_binary_constant': 'UNRESOLVED', 'independent_review': 'NOT_DONE',
              'bibliographic_novelty': 'NOT_ESTABLISHED'}
    report['family'] = audit_family()
    print('PASS: family realizations, commutativity, spectra, full word support and trace counting')
    report['homogeneous_fixtures'] = audit_arbitrary_homogeneous()
    print('PASS: arbitrary signed homogeneous fixtures, all split Frobenius identities')
    report['table'] = audit_algebra()
    print('PASS: general algebra, limit and exact numerical examples')
    args.output.write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print('Written:', args.output)


if __name__ == '__main__':
    main()
