#!/usr/bin/env python3
"""Finite rational checks for the AAK fixed-(d,k) proof note.

No floating-point arithmetic, external packages, or quantifier elimination.
These checks exercise identities, not the infinite-dimensional theorem.
"""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import json


def matrix(rows):
    return [[F(x) for x in row] for row in rows]


def zeros(n, m=None):
    return [[F(0) for _ in range(n if m is None else m)] for _ in range(n)]


def eye(n):
    a = zeros(n)
    for i in range(n):
        a[i][i] = F(1)
    return a


def transpose(a):
    return [list(row) for row in zip(*a)]


def add(a, b):
    return [[x + y for x, y in zip(ar, br)] for ar, br in zip(a, b)]


def scale(a, c):
    return [[F(c) * x for x in row] for row in a]


def mul(a, b):
    bt = transpose(b)
    return [[sum((x * y for x, y in zip(row, col)), F(0)) for col in bt] for row in a]


def outer(v):
    return [[x * y for y in v] for x in v]


def trace(a):
    return sum((a[i][i] for i in range(len(a))), F(0))


def determinant(a):
    a = [row[:] for row in a]
    value = F(1)
    for i in range(len(a)):
        pivot = next((j for j in range(i, len(a)) if a[j][i]), None)
        if pivot is None:
            return F(0)
        if pivot != i:
            a[i], a[pivot] = a[pivot], a[i]
            value *= -1
        value *= a[i][i]
        for j in range(i + 1, len(a)):
            factor = a[j][i] / a[i][i]
            for l in range(i, len(a)):
                a[j][l] -= factor * a[i][l]
    return value


def solve(a, b):
    augmented = [row[:] + [rhs] for row, rhs in zip(a, b)]
    n = len(a)
    for i in range(n):
        pivot = next((j for j in range(i, n) if augmented[j][i]), None)
        if pivot is None:
            raise ValueError("Singular exact system")
        augmented[i], augmented[pivot] = augmented[pivot], augmented[i]
        divisor = augmented[i][i]
        augmented[i] = [x / divisor for x in augmented[i]]
        for j in range(n):
            if j != i:
                factor = augmented[j][i]
                augmented[j] = [x - factor * y for x, y in zip(augmented[j], augmented[i])]
    return [row[-1] for row in augmented]


def positive_map(transitions, x, adjoint):
    result = zeros(len(x))
    for a in transitions:
        at = transpose(a)
        term = mul(mul(at, x), a) if adjoint else mul(mul(a, x), at)
        result = add(result, term)
    return result


def stein_matrix(transitions, adjoint=True):
    n = len(transitions[0])
    operator = eye(n * n)
    for i in range(n):
        for j in range(n):
            unit = zeros(n)
            unit[i][j] = F(1)
            image = positive_map(transitions, unit, adjoint)
            for r in range(n):
                for c in range(n):
                    operator[r * n + c][i * n + j] -= image[r][c]
    return operator


def stein(transitions, rhs, adjoint=True):
    n = len(rhs)
    operator = stein_matrix(transitions, adjoint)
    solution = solve(operator, [x for row in rhs for x in row])
    return [solution[i * n:(i + 1) * n] for i in range(n)]


def characteristic_plus_coeffs(a):
    """Coefficients e_1,...,e_n in det(s I + a), exactly."""
    n = len(a)
    b = eye(n)
    coefficients = []
    for j in range(1, n + 1):
        ab = mul(a, b)
        cj = -trace(ab) / j
        coefficients.append(((-1) ** j) * cj)
        b = add(ab, scale(eye(n), cj))
    return coefficients


def psd_via_coeffs(a):
    assert a == transpose(a), "Symmetry is required"
    return all(c >= 0 for c in characteristic_plus_coeffs(a))


def psd_via_all_principal_minors(a):
    assert a == transpose(a)
    n = len(a)
    for size in range(1, n + 1):
        for indices in combinations(range(n), size):
            minor = [[a[i][j] for j in indices] for i in indices]
            if determinant(minor) < 0:
                return False
    return True


def block_diag(a, b):
    n, k = len(a), len(b)
    result = zeros(n + k)
    for i in range(n):
        for j in range(n):
            result[i][j] = a[i][j]
    for i in range(k):
        for j in range(k):
            result[n + i][n + j] = b[i][j]
    return result


REPORT = []


def check(name, condition, details=None):
    REPORT.append({"check": name, "passed": bool(condition), "details": details})
    if not condition:
        raise AssertionError(name)


def certify(transitions, p, q, t, label):
    n = len(p)
    t = F(t)
    s = add(scale(p, t * t), scale(mul(mul(p, q), p), -1))
    by_coeffs = psd_via_coeffs(s)
    by_minors = psd_via_all_principal_minors(s)
    check(label + ": PSD encodings agree", by_coeffs == by_minors)
    kp = stein_matrix(transitions, True)
    kq = stein_matrix(transitions, False)
    check(label + ": Stein coefficient matrices transpose", kp == transpose(kq))
    delta = determinant(kp)
    check(label + ": Stein denominator nonzero", delta != 0)
    for denominator in (delta, F(-3)):
        u, v = scale(p, denominator), scale(q, denominator)
        cleared = scale(add(scale(u, t * t * denominator * denominator),
                            scale(mul(mul(u, v), u), -1)), denominator)
        check(label + ": denominator identity " + str(denominator),
              cleared == scale(s, denominator ** 4))
        check(label + ": PSD preserved " + str(denominator),
              psd_via_coeffs(cleared) == by_coeffs)
    return by_coeffs


def error_gramians(target, alpha, beta, candidate, gamma, nu, label):
    transitions = [block_diag(a, b) for a, b in zip(target, candidate)]
    c = [F(x) for x in alpha] + [-F(x) for x in gamma]
    b = [F(x) for x in beta] + [F(x) for x in nu]
    p = stein(transitions, outer(c), True)
    q = stein(transitions, outer(b), False)
    check(label + ": observability Stein identity", p == add(outer(c), positive_map(transitions, p, True)))
    check(label + ": reachability Stein identity", q == add(outer(b), positive_map(transitions, q, False)))
    check(label + ": P and Q symmetric PSD", psd_via_coeffs(p) and psd_via_coeffs(q))
    return transitions, p, q


def main():
    j = matrix([[0, 1, 0], [0, 0, 1], [0, 0, 0]])
    target = [j, mul(j, j)]
    alpha, beta = [1, 0, 0], [1, 0, 1]

    # Complete finite support is (epsilon, a, b, aa); nilpotency kills longer words.
    h = matrix([[1, 0, 1, 1], [0, 1, 0, 0], [1, 0, 0, 0], [1, 0, 0, 0]])
    check("target squared singular values 4,1,1,0",
          characteristic_plus_coeffs(mul(h, h)) == [F(6), F(9), F(4), F(0)])
    h_epsilon = [row[:] for row in h]
    h_epsilon[0][0] = F(0)
    check("epsilon candidate squared singular values 2,2,1,0",
          characteristic_plus_coeffs(mul(h_epsilon, h_epsilon)) == [F(5), F(8), F(4), F(0)])

    transitions, p, q = error_gramians(target, alpha, beta, [matrix([[0]]), matrix([[0]])],
                                      [0], [0], "zero candidate")
    check("target norm <=2", certify(transitions, p, q, 2, "threshold 2"))
    check("target norm >3/2", not certify(transitions, p, q, F(3, 2), "threshold 3/2"))

    transitions, p, q = error_gramians(target, alpha, beta, [matrix([[0]]), matrix([[0]])],
                                      [1], [1], "epsilon candidate")
    check("epsilon candidate error >7/5", not certify(transitions, p, q, F(7, 5), "epsilon 7/5"))
    check("epsilon candidate error <=3/2", certify(transitions, p, q, F(3, 2), "epsilon 3/2"))

    transitions, p, q = error_gramians(target, alpha, beta, target, alpha, beta, "identical candidate")
    check("identical candidate P singular", determinant(p) == 0)
    check("identical candidate error <=0", certify(transitions, p, q, 0, "identical t=0"))

    stable = [matrix([[F(1, 2)]]), matrix([[F(1, 4)]])]
    x = stein(stable, eye(1), True)
    check("scalar stable Lyapunov solution", x == matrix([[F(16, 11)]]))
    check("scalar stable solution positive", determinant(x) > 0)
    unstable = [matrix([[F(6, 5)]]), matrix([[0]])]
    x = stein(unstable, eye(1), True)
    check("scalar unstable Lyapunov solution", x == matrix([[F(-25, 11)]]))
    check("scalar unstable rejected", determinant(x) < 0)

    for d in (2, 3, 4):
        word = (0,)
        output = {}
        for letter in range(d):
            h = word + (letter,)
            shifted = h[:-1] if h[-1] == letter else None
            if shifted is not None:
                output[shifted] = output.get(shifted, F(0)) + 1
        lhs = sum((c * c for c in output.values()), F(0))
        rhs = F(d)
        check("backward row contraction fails d=" + str(d), lhs == d*d and lhs > rhs,
              {"lhs_squared": str(lhs), "rhs_squared": str(rhs)})

    # Exact regression of the PSD characterization, including zero eigenvalues.
    for idx, a in enumerate((matrix([[0, 0], [0, 1]]), matrix([[1, 2], [2, 1]]),
                             matrix([[1, 1], [1, 1]]), matrix([[0, 1], [1, 0]]))):
        check("PSD characterization sample " + str(idx),
              psd_via_coeffs(a) == psd_via_all_principal_minors(a))

    result = {"status": "PASS", "arithmetic": "exact rational Fractions",
              "checks": len(REPORT), "global_QE_implemented": False,
              "general_proof_formalized": False, "results": REPORT}
    path = Path(__file__).with_name("exact_identity_checks.json")
    path.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in result.items() if k != "results"}, ensure_ascii=False))


if __name__ == "__main__":
    main()
