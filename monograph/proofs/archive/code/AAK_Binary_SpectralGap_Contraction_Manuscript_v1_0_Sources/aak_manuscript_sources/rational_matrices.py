"""Exact rational matrix helpers (Python standard library)."""
from fractions import Fraction as F
from itertools import combinations

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

