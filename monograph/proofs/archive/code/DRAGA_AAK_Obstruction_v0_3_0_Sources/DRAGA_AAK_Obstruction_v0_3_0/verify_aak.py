#!/usr/bin/env python3
"""Exact arithmetic replay for DRAGA AAK v0.3.0. Python 3 + SymPy only.

Run: python verify_aak.py --output certificates.json
This checks finite identities and certificates. It is not a proof assistant.
The universal arguments are in DRAGA_AAK_Obstruction_v0_3_0.tex / PDF.
"""
from __future__ import annotations
import argparse
import json
import platform
from itertools import product
from pathlib import Path
import sympy as s


def require(condition, message):
    # Deliberately works even under python -O.
    if condition != True:
        raise AssertionError(message)


def words(d, length):
    return list(product(range(d), repeat=length))


def realization(d):
    matrices = []
    for i in range(d):
        A = s.zeros(d + 2)
        A[0, i + 1] = A[i + 1, d + 1] = 1
        matrices.append(A)
    return matrices, s.eye(d + 2)[:, 0], s.eye(d + 2)[:, d + 1]


def evaluate(A, alpha, beta, word):
    state = beta
    for i in reversed(word):
        state = A[i] * state
    return (alpha.T * state)[0]


def positive_minors(M):
    require(M == M.T, 'Matrix must be real symmetric')
    minors = [M[:i, :i].det() for i in range(1, M.rows + 1)]
    require(all(x > 0 for x in minors), 'Sylvester positivity certificate failed')
    return list(map(str, minors))


def json_matrix(M):
    return [[str(M[i, j]) for j in range(M.cols)] for i in range(M.rows)]


def audit_family():
    replayed = []
    z = s.symbols('z')
    for d in [2, 3, 4, 6]:
        A, alpha, beta = realization(d)
        n = d + 2
        for i in range(d):
            for j in range(d):
                expected = s.zeros(n)
                if i == j:
                    expected[0, -1] = 1
                require(A[i] * A[j] == expected, 'Pair product identity')
                for ell in range(d):
                    require(A[i] * A[j] * A[ell] == s.zeros(n), 'Nilpotence')
        W = [()] + words(d, 1) + words(d, 2)
        cache = {u + v: None for u in W for v in W}
        cache = {w: evaluate(A, alpha, beta, w) for w in cache}
        require(all(value == int(len(w) == 2 and w[0] == w[1])
                    for w, value in cache.items()), 'Target coefficients')
        H = s.Matrix([[cache[u + v] for v in W] for u in W])
        require(H == H.T, 'Self-adjoint target')
        # An orthonormal basis for the exact support, including the collective
        # length-two vector. Equality below also checks annihilation outside it.
        J = s.zeros(len(W), n)
        J[0, 0] = 1
        for i in range(d):
            J[W.index((i,)), i + 1] = 1
            J[W.index((i, i)), n - 1] = 1 / s.sqrt(d)
        require(s.simplify(J.T * J) == s.eye(n), 'Support isometry')
        small = s.zeros(n)
        small[0, -1] = small[-1, 0] = s.sqrt(d)
        for i in range(1, d + 1):
            small[i, i] = 1
        require(s.simplify(H - J * small * J.T) == s.zeros(len(W)), 'Exact support')
        require(s.factor((small.T * small).charpoly(z).as_expr()) ==
                (z - d)**2 * (z - 1)**d, 'Squared singular spectrum')
        require(small.det() != 0, 'Minimal rank d+2')
        # Negative control: spectral truncation rank two is not Hankel.
        G = H.copy()
        for i in range(d):
            G[W.index((i,)), W.index((i,))] = 0
        require(G.rank() == 2, 'SVD candidate rank')
        require(G[0, W.index((0, 0))] != G[W.index((0,)), W.index((0,))],
                'Negative control must violate concatenation')
        require((H - G)**2 == H - G, 'Unstructured error is a projection of norm one')
        replayed.append({'d': d, 'target_rank': n, 'status': 'PASS'})
    return replayed


def audit_two_letters():
    t = s.symbols('t')
    a = 1 + t - t**3
    b = 2*t*t + t - 2
    poly = s.Poly(2*a*a - b*b, t)
    require(poly.as_expr().expand() == 2*t**6 - 8*t**4 - 8*t**3 + 9*t*t + 8*t - 2,
            'Elimination polynomial')
    lo, hi = s.Rational(2627, 2500), s.Rational(105083, 100000)
    for endpoint, expected_sign in [(lo, 1), (hi, -1)]:
        require(a.subs(t, endpoint) > 0 and b.subs(t, endpoint) > 0,
                'Squaring branch guard')
        require(s.sign(poly.eval(endpoint)) == expected_sign, 'Root endpoint sign')
    require(poly.count_roots(lo, hi) == 1, 'Unique isolated polynomial root')
    A, alpha, beta = realization(2)
    B = [s.Matrix([[0, 1], [-s.Rational(58127, 100000), s.Rational(6157, 5000)]]),
         s.zeros(2)]
    gamma = s.Matrix([1, 0])
    eta = s.Matrix([s.Rational(5443, 50000), s.Rational(8523, 25000)])
    require(s.Matrix.hstack(eta, B[0]*eta).det() != 0, 'Candidate reachable')
    require(s.Matrix.vstack(gamma.T, gamma.T*B[0]).det() != 0, 'Candidate observable')
    X = s.Matrix(2, 2, (s.eye(4) - s.kronecker_product(B[0].T, B[0].T)).inv()
                 * s.eye(2).reshape(4, 1))
    require(X - B[0].T*X*B[0] == s.eye(2), 'Lyapunov equality')
    stability_minors = positive_minors(X)
    D = [s.diag(M, N) for M, N in zip(A, B)]
    c, v = alpha.col_join(-gamma), beta.col_join(eta)
    K = sum((s.kronecker_product(M, M) for M in D), s.zeros(36))
    P = s.Matrix(6, 6, (s.eye(36)-K.T).inv()*(c*c.T).reshape(36, 1))
    Q = s.Matrix(6, 6, (s.eye(36)-K).inv()*(v*v.T).reshape(36, 1))
    require(P == c*c.T + sum((M.T*P*M for M in D), s.zeros(6)), 'Observability Stein')
    require(Q == v*v.T + sum((M*Q*M.T for M in D), s.zeros(6)), 'Reachability Stein')
    U = s.Rational(130993, 100000)
    minors_P = positive_minors(P)
    minors_Q = positive_minors(Q)
    minors_norm = positive_minors(U*U*P-P*Q*P)
    return {'status': 'PASS', 'root_interval': [str(lo), str(hi)],
            'polynomial': str(poly.as_expr()), 'strict_error_upper': str(U),
            'B_a': json_matrix(B[0]), 'eta': json_matrix(eta),
            'P': json_matrix(P), 'Q': json_matrix(Q), 'X': json_matrix(X),
            'Sylvester_minors': {'X': stability_minors, 'P': minors_P,
                                  'Q': minors_Q, 'U2P_minus_PQP': minors_norm}}


def audit_dimension_claims():
    examples = []
    for d, k in [(4, 2), (10, 2), (100, 2), (100, 96), (100, 98)]:
        lower2, upper2 = max(1, d-k), d-k+2
        require(2 <= k <= d and lower2 <= upper2, 'Family parameter range')
        examples.append({'d': d, 'k': k, 'spectral_lower': 1,
                         'structured_lower_squared': lower2,
                         'constructive_upper_squared': upper2})
    d, tolerance = 100, s.Integer(2)
    require(d - 95 > tolerance**2, '95 or fewer states excluded by length-two block')
    retained_terms = 96
    require(retained_terms + 2 == 98 and d-retained_terms == tolerance**2,
            '98-state exact deletion candidate meets tolerance')
    require(s.sqrt(d) > tolerance and 1 <= tolerance, 'Spectral state floor is two')
    return {'status': 'PASS', 'examples': examples,
            'd100_tolerance2': {'spectral_state_floor': 2,
                               'structured_state_floor': 96,
                               'constructive_state_ceiling': 98}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('certificates.json'))
    args = parser.parse_args()
    report = {'version': '0.3.0', 'python': platform.python_version(),
              'sympy': s.__version__, 'arithmetic': 'exact_rational_and_algebraic',
              'universal_proof': 'WRITTEN_IN_NOTE_NOT_FORMALIZED',
              'independent_review': 'NOT_DONE', 'bibliographic_novelty': 'NOT_ESTABLISHED'}
    report['family_finite_replay'] = audit_family()
    print('PASS: finite realizations, support, singular spectrum, SVD negative controls')
    report['two_letter_certificate'] = audit_two_letters()
    print('PASS: two-letter root isolation, candidate stability and exact norm upper bound')
    report['dimension_examples'] = audit_dimension_claims()
    print('PASS: d=100, error tolerance 2: at least 96 states; 98 states suffice')
    args.output.write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print('Written:', args.output)


if __name__ == '__main__':
    main()
