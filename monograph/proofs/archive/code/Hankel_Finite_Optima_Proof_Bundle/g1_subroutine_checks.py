#!/usr/bin/env python3
"""Reproduce small exact-arithmetic checks for the G1 construction.

Run with Python 3 and SymPy (mpmath is a SymPy dependency):
    python g1_subroutine_checks.py
    python g1_subroutine_checks.py --output /path/to/checks.json

All pass/fail decisions use rational arithmetic.  The infinite Hankel
Hilbert--Schmidt sums and Gramians are evaluated by exact Stein equations;
no finite word block is substituted for the infinite operator.  The
high-precision eigenvalue calculation is a separately labelled diagnostic.

This checks subroutines and examples, not global optimality or the complete
net.  The parameter-net count is reported as an UNEXECUTED upper bound.
"""

from __future__ import annotations

import argparse
import itertools
import json
import platform
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path

import mpmath as mp
import sympy as sp


Q = sp.Rational
ASSERTIONS: list[str] = []


def require(name: str, condition: object) -> None:
    if not bool(condition):
        raise AssertionError(name)
    ASSERTIONS.append(name)


def frac(value: object) -> Fraction:
    value = sp.cancel(value)
    return Fraction(int(sp.numer(value)), int(sp.denom(value)))


def qstr(value: object) -> str:
    return str(frac(value))


def matrix_json(value: sp.Matrix) -> list[list[str]]:
    return [[qstr(value[i, j]) for j in range(value.cols)]
            for i in range(value.rows)]


def matrix_equal(left: sp.Matrix, right: sp.Matrix) -> bool:
    return all(sp.cancel(x) == 0 for x in left - right)


def vec(matrix: sp.Matrix) -> sp.Matrix:
    """Column-major vectorization, matching conjugate(A) tensor A."""
    return sp.Matrix([matrix[i, j]
                      for j in range(matrix.cols)
                      for i in range(matrix.rows)])


def unvec(vector: sp.Matrix, n: int) -> sp.Matrix:
    return sp.Matrix(n, n, lambda i, j: vector[i + n*j])


def cp(letters: tuple[sp.Matrix, ...], matrix: sp.Matrix,
       dual: bool = False) -> sp.Matrix:
    terms = ((a.H*matrix*a if dual else a*matrix*a.H) for a in letters)
    return sum(terms, sp.zeros(matrix.rows))


def cp_matrix(letters: tuple[sp.Matrix, ...],
              dual: bool = False) -> sp.Matrix:
    n = letters[0].rows
    chosen = tuple(a.H for a in letters) if dual else letters
    return sum((sp.kronecker_product(a.conjugate(), a) for a in chosen),
               sp.zeros(n*n))


def stein(letters: tuple[sp.Matrix, ...], source: sp.Matrix,
          dual: bool = False) -> sp.Matrix:
    n = source.rows
    result = unvec((sp.eye(n*n) - cp_matrix(letters, dual)).inv()*vec(source), n)
    result = result.applyfunc(sp.cancel)
    require("exact Stein equation", matrix_equal(
        result, source + cp(letters, result, dual)))
    return result


def psd_exact(matrix: sp.Matrix) -> bool:
    """All principal minors; used only for matrices of order at most four."""
    if not matrix_equal(matrix, matrix.H):
        return False
    n = matrix.rows
    for size in range(1, n + 1):
        for indices in itertools.combinations(range(n), size):
            minor = sp.cancel(matrix.extract(indices, indices).det())
            if minor.is_real is not True or minor < 0:
                return False
    return True


@dataclass(frozen=True)
class WFA:
    letters: tuple[sp.Matrix, ...]
    c: sp.Matrix
    b: sp.Matrix

    @property
    def n(self) -> int:
        return self.b.rows

    def json(self) -> dict:
        return {"state_dimension": self.n,
                "alphabet_size": len(self.letters),
                "letters": [matrix_json(a) for a in self.letters],
                "readout_c": matrix_json(self.c),
                "boundary_b": matrix_json(self.b)}


def difference(first: WFA, second: WFA) -> WFA:
    require("equal alphabet size", len(first.letters) == len(second.letters))
    letters = tuple(sp.diag(a, b) for a, b in zip(first.letters, second.letters))
    return WFA(letters, first.c.row_join(-second.c), first.b.col_join(second.b))


def damp(model: WFA, gamma: sp.Rational) -> WFA:
    return WFA(tuple(gamma*a for a in model.letters), model.c, model.b)


def gramian_data(model: WFA, tag: str) -> tuple[sp.Matrix, sp.Matrix, sp.Rational]:
    p = stein(model.letters, model.b*model.b.H)
    q = stein(model.letters, model.c.H*model.c, dual=True)
    require(tag + ": reachability Gramian PSD", psd_exact(p))
    require(tag + ": observability Gramian PSD", psd_exact(q))
    hs_sq = sp.cancel(sp.trace(p*q))
    # Independent weighted sum: (I-Phi)^(-2)(bb*) counts all word splits.
    weighted = stein(model.letters, p)
    weighted_hs_sq = sp.cancel((model.c*weighted*model.c.H)[0])
    require(tag + ": two-Stein weighted HS equals trace(PQ)", hs_sq == weighted_hs_sq)
    return p, q, hs_sq


def moment_check(model: WFA, tag: str) -> dict:
    moments = []
    current = model.b*model.b.H
    for _ in range(4):
        current = stein(model.letters, current)
        moments.append(current)
    hs_sq = sp.cancel((model.c*moments[1]*model.c.H)[0])
    j_sq = sp.cancel((model.c*(2*cp(model.letters, moments[2])
                        + 6*cp(model.letters, cp(model.letters, moments[3])))
                        *model.c.H)[0])

    # Separate scalar generating function F(z)=sum_t energy_t z^t.
    # F(1)+F'(1)=sum(t+1)energy_t.
    # F'''(1)+4F''(1)+2F'(1)=sum t^2(t+1)energy_t.
    z = sp.Symbol("z")
    resolvent = (sp.eye(model.n**2) - z*cp_matrix(model.letters)).inv()
    fz = sp.cancel((model.c*unvec(resolvent*vec(model.b*model.b.H), model.n)
                    *model.c.H)[0])
    derivatives = [sp.cancel(sp.diff(fz, z, j).subs(z, 1)) for j in range(4)]
    hs_independent = sp.cancel(derivatives[0] + derivatives[1])
    j_independent = sp.cancel(derivatives[3] + 4*derivatives[2] + 2*derivatives[1])
    require(tag + ": four-Stein HS formula", hs_sq == hs_independent)
    require(tag + ": four-Stein J2 formula", j_sq == j_independent)
    return {"model": model.json(), "generating_function_energy": str(fz),
            "four_stein_matrices": [matrix_json(x) for x in moments],
            "HS_squared": qstr(hs_sq), "J2_squared": qstr(j_sq),
            "generating_function_HS_squared": qstr(hs_independent),
            "generating_function_J2_squared": qstr(j_independent),
            "all_equalities_exact": True}


def power_le_root(value: Fraction, target: Fraction, degree: int) -> bool:
    return value.numerator**degree * target.denominator <= (
        target.numerator * value.denominator**degree)


def root_bracket(target: Fraction, degree: int, upper: Fraction,
                 tolerance: Fraction) -> tuple[Fraction, Fraction, int]:
    """Rational bisection for a positive scalar root; no floating point."""
    require("root target nonnegative", target >= 0)
    if target == 0:
        return Fraction(0), Fraction(0), 0
    require("root initial upper bound", not power_le_root(upper, target, degree)
            or upper**degree == target)
    lower = Fraction(0)
    iterations = 0
    while upper - lower > tolerance:
        middle = (lower + upper)/2
        if power_le_root(middle, target, degree):
            lower = middle
        else:
            upper = middle
        iterations += 1
    require("root lower rational power check", lower**degree <= target)
    require("root upper rational power check", upper**degree >= target)
    return lower, upper, iterations


def mp_fraction(value: object) -> mp.mpf:
    value = frac(value)
    return mp.mpf(value.numerator)/value.denominator


def diagnostic_norm(p: sp.Matrix, q: sp.Matrix) -> str:
    """Floating calculation, explicitly excluded from every assertion."""
    with mp.workdps(80):
        product = p*q
        matrix = mp.matrix([[mp_fraction(product[i, j])
                             for j in range(product.cols)]
                            for i in range(product.rows)])
        eigenvalues = mp.eig(matrix, left=False, right=False)
        maximum = max(mp.re(value) for value in eigenvalues)
        return mp.nstr(mp.sqrt(maximum), 65)


def trace_certificate(p: sp.Matrix, q: sp.Matrix, exponent: int,
                      uniform_bound: Fraction, tag: str) -> dict:
    dimension = p.rows
    trace = frac(sp.trace((p*q)**exponent))
    require(tag + ": trace nonnegative", trace >= 0)
    degree = 2*exponent
    tolerance = Fraction(1, 2**42)
    lo0, lo1, iterations_lower = root_bracket(
        trace/dimension, degree, 2*uniform_bound, tolerance)
    up0, up1, iterations_upper = root_bracket(
        trace, degree, 2*uniform_bound, tolerance)
    require(tag + ": ordered norm interval", lo0 <= up1)
    # PSD P,Q make all eigenvalues of PQ nonnegative.  Thus
    # (tr((PQ)^p)/m)^(1/(2p)) <= ||H|| <= tr((PQ)^p)^(1/(2p)).
    return {"state_dimension_for_error": dimension, "p_trace": exponent,
            "root_degree": degree, "trace_PQ_power": str(trace),
            "lower_exact_endpoint_root_bracket": [str(lo0), str(lo1)],
            "upper_exact_endpoint_root_bracket": [str(up0), str(up1)],
            "certified_infinite_norm_interval": [str(lo0), str(up1)],
            "interval_width": str(up1-lo0),
            "bisections": {"lower": iterations_lower, "upper": iterations_upper},
            "each_root_bracket_width_at_most": str(tolerance),
            "rational_power_comparisons_passed": True}


def ceil_sqrt_rational(value: object) -> int:
    value = frac(value)
    integer = sp.integer_nthroot(value.numerator//value.denominator, 2)[0]
    while integer*integer*value.denominator < value.numerator:
        integer += 1
    return int(integer)


def supplementary_word_check(first: WFA, second: WFA, gamma: sp.Rational,
                             eta: sp.Rational, boundary_bound: sp.Rational) -> list[dict]:
    rows = []
    for length in range(6):
        energy = 0
        for word in itertools.product(range(len(first.letters)), repeat=length):
            first_product, second_product = sp.eye(first.n), sp.eye(second.n)
            for letter in word:
                first_product *= first.letters[letter]
                second_product *= second.letters[letter]
            coefficient = (first.c*first_product*first.b
                           - second.c*second_product*second.b)[0]
            energy += coefficient*sp.conjugate(coefficient)
        bound = boundary_bound*eta*(2*gamma**length + (
            length*gamma**(length-1) if length else 0))
        require(f"length {length}: coefficient-vector Lipschitz estimate",
                sp.cancel(energy) <= bound**2)
        rows.append({"length": length, "words_enumerated": len(first.letters)**length,
                     "exact_coefficient_energy": qstr(energy),
                     "squared_length_bound": qstr(bound**2)})
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=Path(__file__).with_suffix(".json"))
    args = parser.parse_args()

    scalar = WFA((sp.Matrix([[Q(1, 3)]]),), sp.Matrix([[Q(3, 4)]]),
                 sp.Matrix([Q(2, 3)]))
    model = WFA((sp.Matrix([[Q(1, 3), Q(1, 5)], [0, Q(1, 4)]]),
                 sp.Matrix([[0, Q(1, 6)], [Q(1, 7), 0]])),
                sp.Matrix([[1, -Q(1, 3)]]), sp.Matrix([Q(2, 3), Q(1, 2)]))
    require("noncommuting input letters", model.letters[0]*model.letters[1]
            != model.letters[1]*model.letters[0])
    for name, current in (("scalar", scalar), ("noncommuting", model)):
        require(name + ": strict joint column stability", psd_exact(
            sp.eye(current.n) - sum((a.H*a for a in current.letters),
                                   sp.zeros(current.n))))
        # The exact principal minors are strictly positive here, so the
        # transition stack has norm strictly below one.
        require(name + ": positive stability determinant", (
            sp.eye(current.n) - sum((a.H*a for a in current.letters),
                                   sp.zeros(current.n))).det() > 0)
    moments = {"scalar_one_letter": moment_check(scalar, "scalar"),
               "noncommuting_two_letter": moment_check(model, "noncommuting")}
    scalar_x, scalar_amplitude = Q(1, 9), Q(1, 2)
    require("scalar HS closed form", frac(moments["scalar_one_letter"]["HS_squared"])
            == frac(scalar_amplitude**2/(1-scalar_x)**2))
    require("scalar J2 closed form", frac(moments["scalar_one_letter"]["J2_squared"])
            == frac(scalar_amplitude**2*2*scalar_x*(1+2*scalar_x)/(1-scalar_x)**4))

    gamma = Q(3, 4)
    error = difference(model, damp(model, gamma))
    radial_p, radial_q, radial_hs_sq = gramian_data(error, "radial error")
    j_sq = Q(moments["noncommuting_two_letter"]["J2_squared"])
    moment_bias_sq = sp.cancel((1-gamma)**2*j_sq)
    j = max(1, ceil_sqrt_rational(j_sq))
    rational_bias_bound = frac((1-gamma)*j)
    require("radial infinite HS bias inequality", radial_hs_sq <= moment_bias_sq)
    radial_intervals = [trace_certificate(radial_p, radial_q, p,
                                         rational_bias_bound, f"radial p={p}")
                        for p in (16, 32)]
    for certificate in radial_intervals:
        certified_upper = Fraction(certificate["certified_infinite_norm_interval"][1])
        require("radial trace interval implies input bias bound",
                certified_upper <= rational_bias_bound)
        # This additional comparison is sharper than the rounded J bound.
        require("radial trace upper squared below raw moment bound",
                certified_upper**2 <= frac(moment_bias_sq))
    radial = {"gamma": qstr(gamma), "J_rational_upper_bound": j,
              "infinite_error_HS_squared": qstr(radial_hs_sq),
              "raw_input_moment_bias_squared": qstr(moment_bias_sq),
              "rational_operator_bias_bound_(1-gamma)J": str(rational_bias_bound),
              "infinite_HS_inequality_exact": True,
              "P": matrix_json(radial_p), "Q": matrix_json(radial_q),
              "trace_power_certificates": radial_intervals,
              "high_precision_norm_DIAGNOSTIC_ONLY": diagnostic_norm(radial_p, radial_q)}

    first = WFA((sp.Matrix([[Q(1, 3), Q(1, 10)], [0, Q(1, 4)]]),
                 sp.Matrix([[0, Q(1, 8)], [Q(1, 7), 0]])),
                sp.Matrix([[Q(2, 3), Q(1, 5)]]), sp.Matrix([1, Q(1, 3)]))
    second = WFA((first.letters[0] + sp.diag(Q(1, 100), -Q(1, 120)),
                  first.letters[1] + sp.Matrix([[0, Q(1, 150)], [-Q(1, 200), 0]])),
                 first.c + sp.Matrix([[Q(1, 100), -Q(1, 150)]]),
                 first.b + sp.Matrix([Q(1, 75), -Q(1, 100)]))
    gamma, eta, boundary_bound = Q(1, 2), Q(1, 50), Q(2)
    mu = 1-gamma**2
    for name, current in (("first", first), ("second", second)):
        require(name + ": joint column gamma bound", psd_exact(
            gamma**2*sp.eye(current.n) - sum((a.H*a for a in current.letters),
                                           sp.zeros(current.n))))
        require(name + ": readout unit bound", (current.c*current.c.H)[0] <= 1)
        require(name + ": boundary B bound", (current.b.H*current.b)[0] <= boundary_bound**2)
    transition_delta = sum(((a-b).H*(a-b) for a, b in zip(first.letters, second.letters)),
                           sp.zeros(first.n))
    require("joint transition difference eta bound", psd_exact(
        eta**2*sp.eye(first.n) - transition_delta))
    require("readout difference eta bound", (
        (first.c-second.c)*(first.c-second.c).H)[0] <= eta**2)
    require("boundary difference B eta bound", (
        (first.b-second.b).H*(first.b-second.b))[0] <= (boundary_bound*eta)**2)
    lip_error = difference(first, second)
    lip_p, lip_q, lip_hs_sq = gramian_data(lip_error, "Lipschitz error")
    lip_bound = sp.cancel(5*boundary_bound*eta/mu**2)
    require("infinite joint-column Lipschitz HS inequality", lip_hs_sq <= lip_bound**2)
    lip_intervals = [trace_certificate(lip_p, lip_q, p, frac(lip_bound),
                                     f"Lipschitz p={p}") for p in (16, 32)]
    for certificate in lip_intervals:
        require("Lipschitz trace interval implies operator bound", Fraction(
            certificate["certified_infinite_norm_interval"][1]) <= frac(lip_bound))
    lipschitz = {"first_model": first.json(), "second_model": second.json(),
                 "gamma": qstr(gamma), "mu": qstr(mu), "eta": qstr(eta),
                 "B": qstr(boundary_bound),
                 "infinite_error_HS_squared": qstr(lip_hs_sq),
                 "rational_infinite_operator_bound_5Beta_over_mu_squared": qstr(lip_bound),
                 "infinite_HS_inequality_exact": True,
                 "P": matrix_json(lip_p), "Q": matrix_json(lip_q),
                 "trace_power_certificates": lip_intervals,
                 "supplementary_finite_length_checks_not_used_for_infinite_certificate":
                     supplementary_word_check(first, second, gamma, eta, boundary_bound),
                 "high_precision_norm_DIAGNOSTIC_ONLY": diagnostic_norm(lip_p, lip_q)}

    # Report the theorem's net size; do NOT enumerate it or its norm ledger.
    epsilon, k, d = Q(1, 10), 1, len(model.letters)
    m = max(1, ceil_sqrt_rational(Q(moments["noncommuting_two_letter"]["HS_squared"])))
    j = max(1, ceil_sqrt_rational(j_sq))
    theorem_gamma = 1-epsilon/(3*j)
    theorem_mu = 1-theorem_gamma**2
    theorem_eta = epsilon*theorem_mu**2/(15*(2*m))
    coordinate_dimension = 2*d*k*k + 4*k
    v_a = int(sp.ceiling(sp.sqrt(2*d*k*k)))
    base_ceiling = int(sp.ceiling(3+2*v_a/theorem_eta))
    count_bound = base_ceiling**coordinate_dimension
    simpler_base = int(sp.ceiling(3+240*v_a*m*j*j/epsilon**3))
    simpler_count_bound = simpler_base**coordinate_dimension
    z = Q(m) + Q(2*m)/theorem_mu
    state_dimension_error = model.n + k
    ell = max(1, (state_dimension_error-1).bit_length())
    theoretical_p = int(sp.ceiling(18*z*ell/epsilon))
    grid = {"status": "NOT EXECUTED: upper-bound calculation only",
            "full_grid_candidates_evaluated": 0,
            "input_state_dimension_n": model.n, "reduced_budget_k": k, "alphabet_size_d": d,
            "epsilon": qstr(epsilon), "M": m, "J": j,
            "gamma": qstr(theorem_gamma), "mu": qstr(theorem_mu),
            "eta": qstr(theorem_eta), "real_coordinate_dimension": coordinate_dimension,
            "v_A": v_a, "rounded_up_grid_base": base_ceiling,
            "integer_candidate_count_upper_bound": str(count_bound),
            "decimal_digits_in_count_bound": len(str(count_bound)),
            "simplified_integer_candidate_count_upper_bound": str(simpler_count_bound),
            "Z_uniform_error_norm_bound": qstr(z),
            "trace_exponent_required_by_theorem": theoretical_p,
            "theorem_trace_exponent_executed": False,
            "demonstration_trace_exponents_actually_executed": [16, 32]}

    report = {"schema_version": 1,
              "purpose": "Exact small-example checks of G1 subroutines; no global grid optimization",
              "versions": {"python": platform.python_version(), "sympy": sp.__version__,
                           "mpmath": mp.__version__},
              "arithmetic": {"certification": "SymPy exact rational matrices and fractions.Fraction scalar powers",
                             "floating_point_pass_fail_decisions": 0,
                             "diagnostic_decimal_precision": 80},
              "moment_formula_checks": moments,
              "radial_damping_check": radial, "joint_column_lipschitz_check": lipschitz,
              "unexecuted_parameter_grid_bound": grid,
              "counts": {"input_moment_examples": 2,
                         "infinite_error_realizations": 2,
                         "trace_power_certificates": 4,
                         "full_parameter_grid_candidates_evaluated": 0,
                         "supplementary_words_enumerated": sum(d**t for t in range(6)),
                         "exact_assertions_passed": len(ASSERTIONS)},
              "assertions_passed": ASSERTIONS,
              "status": "PASS"}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"PASS: {len(ASSERTIONS)} exact checks; wrote {args.output}")
    print(f"Full parameter net not executed; candidate-count upper bound has {len(str(count_bound))} decimal digits.")


if __name__ == "__main__":
    main()
