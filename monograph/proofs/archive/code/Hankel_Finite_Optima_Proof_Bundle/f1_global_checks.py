#!/usr/bin/env python3
"""High-precision diagnostics for the exact rank-one binary Hankel optimum.

Run:
    python f1_global_checks.py --output f1_global_checks_results.json

Dependency: mpmath. No word-horizon truncation is used. The infinite Hankel
error norm is evaluated through two finite Stein equations on the exact
target-plus-competitor realization. These computations are diagnostics,
not interval certificates and not the proof of global optimality.

Convention: f(w) = c A_w b, ordinary counting-measure word-Hankel entries.
"""

from __future__ import annotations

import argparse
import json
import platform
from pathlib import Path

import mpmath as mp


Q_CASES = ("0.0001", "0.001", "0.01", "0.1", "0.25", "0.5", "0.75", "0.9", "0.99", "0.9999")


def adj(a):
    return a.transpose_conj()


def maxabs(a):
    return max(abs(a[i, j]) for i in range(a.rows) for j in range(a.cols))


def herm(a):
    return (a + adj(a)) / 2


def number(value):
    return mp.nstr(value, mp.mp.dps)


def complex_number(value):
    return {"real": number(mp.re(value)), "imag": number(mp.im(value))}


def solve_selected_root(q):
    """Bisection of a strictly increasing difference on (q,sqrt(q))."""
    r2 = 1 - q

    def difference(u):
        left = q * (u - q) ** 2 / (r2 ** 2 * u ** 2)
        right = r2 * (q - u ** 2) / (q * (1 - u ** 2))
        return left - right

    lo, hi = q, mp.sqrt(q)
    if not (difference(lo) < 0 and difference(hi) > 0):
        # Rounding can make the analytic zero at the upper endpoint slightly
        # negative; the left contribution there is strictly positive.
        if not (difference(lo) < 0 and abs(difference(hi)) < mp.eps):
            raise AssertionError("Root bracket failed.")
    iterations = 0
    for iterations in range(1, 5 * mp.mp.dps + 65):
        middle = (lo + hi) / 2
        if middle == lo or middle == hi:
            break
        if difference(middle) < 0:
            lo = middle
        else:
            hi = middle
    u = (lo + hi) / 2
    return u, iterations, hi - lo, difference(u)


def exact_branch_parameters(q):
    r2, r = 1 - q, mp.sqrt(1 - q)
    u, iterations, width, monotone_residual = solve_selected_root(q)
    t = (q ** 2 + (1 - 2 * q) * u ** 2) / (q * (1 - u ** 2))
    x = (1 - u) / r
    y = mp.sqrt(u * (1 - u) / r2)
    m = (u - q) / r2
    tau = u / (y * (1 + u))
    c = m * tau
    left = q ** 2 * (1 - u ** 2) * (u - q) ** 2
    right = r2 ** 3 * u ** 2 * (q - u ** 2)
    direct_a = r2 ** 2 * t ** 2 + (q ** 2 + 2 * q - 1) * (t - q)
    direct_rhs = 4 * q ** 3 * (t - q) * (1 - 2 * q + q * t)
    return {"q": q, "r": r, "u": u, "t": t, "x": x, "y": y,
            "m": m, "tau": tau, "c": c,
            "root_iterations": iterations, "root_bracket_width": width,
            "monotone_root_residual": monotone_residual,
            "u_quartic_absolute_residual": abs(left - right),
            "u_quartic_relative_residual": abs(left - right) / max(abs(left), abs(right)),
            "t_quartic_absolute_residual": abs(direct_a ** 2 - direct_rhs),
            "t_quartic_relative_residual": abs(direct_a ** 2 - direct_rhs) /
                                            max(abs(direct_a ** 2), abs(direct_rhs)),
            "direct_unsquared_selector": direct_a,
            "direct_t_interval_lower": 2 * q / (1 + q)}


def stein(letters, source, side):
    """Solve all n^2 scalar equations explicitly; row-major vectorization."""
    n = source.rows
    system = mp.eye(n * n)
    rhs = mp.matrix(n * n, 1)
    for i in range(n):
        for j in range(n):
            ij = i * n + j
            rhs[ij] = source[i, j]
            for k in range(n):
                for l in range(n):
                    if side == "row":
                        coefficient = sum(a[i, k] * mp.conj(a[j, l]) for a in letters)
                    elif side == "column":
                        coefficient = sum(mp.conj(a[k, i]) * a[l, j] for a in letters)
                    else:
                        raise ValueError(side)
                    system[ij, k * n + l] -= coefficient
    vector = mp.lu_solve(system, rhs)
    raw = mp.matrix(n, n)
    for i in range(n):
        for j in range(n):
            raw[i, j] = vector[i * n + j]
    symmetry_residual = maxabs(raw - adj(raw))
    result = herm(raw)
    if side == "row":
        residual = result - source - sum((a * result * adj(a) for a in letters), mp.zeros(n))
    else:
        residual = result - source - sum((adj(a) * result * a for a in letters), mp.zeros(n))
    relative = maxabs(residual) / max(mp.mpf(1), maxabs(result), maxabs(source))
    return result, relative, symmetry_residual


def infinite_error_squared(parameters, theta):
    """Stable target2 plus character1, using its exact infinite Gramians."""
    q, r = parameters["q"], parameters["r"]
    x, y, c = parameters["x"], parameters["y"], parameters["c"]
    phase = mp.exp(1j * theta)
    lambda_b = y * phase
    character_coefficient = c / phase
    aa = mp.diag([r, 0, x])
    ab = mp.matrix([[0, 1, 0], [0, 0, 0], [0, 0, lambda_b]])
    boundary = mp.matrix([0, 1, 1])
    readout = mp.matrix([[q, 0, -character_coefficient]])
    letters = (aa, ab)
    reachability, mr, msym = stein(letters, boundary * adj(boundary), "row")
    observability, qr, qsym = stein(letters, adj(readout) * readout, "column")
    # M=L L*, so the nonzero eigenvalues of M Q equal those of L* Q L.
    lower = mp.cholesky(reachability)
    squared_singulars, _ = mp.eighe(herm(adj(lower) * observability * lower))
    reachability_eigenvalues, _ = mp.eighe(reachability)
    observability_eigenvalues, _ = mp.eighe(observability)
    norm_squared = mp.re(squared_singulars[squared_singulars.rows - 1])
    return {"actual_squared_error": norm_squared,
            "actual_error_norm": mp.sqrt(norm_squared),
            "relative_stein_residual": max(mr, qr),
            "gramian_hermitian_residual": max(msym, qsym),
            "reachability_gramian_operator_norm": mp.re(reachability_eigenvalues[reachability.rows - 1]),
            "observability_gramian_operator_norm": mp.re(observability_eigenvalues[observability.rows - 1]),
            "reachability_gramian_minimum_eigenvalue": mp.re(reachability_eigenvalues[0]),
            "observability_gramian_minimum_eigenvalue": mp.re(observability_eigenvalues[0]),
            "lambda_b": lambda_b, "coefficient_c": character_coefficient}


def eigenvalue_2x2(a, d, offdiag_magnitude):
    return (a + d + mp.sqrt((a - d) ** 2 + 4 * offdiag_magnitude ** 2)) / 2


def run_case(q_string, tolerance):
    p = exact_branch_parameters(mp.mpf(q_string))
    q, u, t, m, y = (p[name] for name in ("q", "u", "t", "m", "y"))
    alpha_squared = q * m / u ** 2
    beta_squared = alpha_squared * y ** 2
    right_squared = eigenvalue_2x2(1 - beta_squared, q * (1 - m),
                                  mp.sqrt(q * beta_squared * m))
    left_squared = eigenvalue_2x2(1 - alpha_squared, q * (1 - beta_squared),
                                 mp.sqrt(q * alpha_squared * beta_squared))
    left_slack_formula = ((1 - u) / ((1 - q) ** 2 * u) *
                          (q * u ** 2 + (1 - q) * u - q ** 2))
    normalization_error = abs(p["x"] ** 2 + y ** 2 + m - 1)
    phases = (("0", mp.mpf(0)), ("pi/3", mp.pi / 3), ("pi", mp.pi))
    phase_records = []
    actuals = []
    largest_norm_difference = mp.mpf(0)
    largest_stein_residual = mp.mpf(0)
    for label, theta in phases:
        measured = infinite_error_squared(p, theta)
        difference = abs(measured["actual_squared_error"] - t)
        largest_norm_difference = max(largest_norm_difference, difference)
        largest_stein_residual = max(largest_stein_residual, measured["relative_stein_residual"])
        actuals.append(measured["actual_squared_error"])
        passed = (difference < tolerance and measured["relative_stein_residual"] < tolerance and
                  measured["reachability_gramian_minimum_eigenvalue"] > 0 and
                  measured["observability_gramian_minimum_eigenvalue"] > 0)
        phase_records.append({"theta": label, "theta_radians": number(theta),
                              **{key: number(value) for key, value in measured.items()
                                 if key not in ("lambda_b", "coefficient_c")},
                              "lambda_b": complex_number(measured["lambda_b"]),
                              "coefficient_c": complex_number(measured["coefficient_c"]),
                              "absolute_squared_error_difference_from_formula": number(difference),
                              "passed": bool(passed)})
    phase_invariance = max(abs(value - actuals[0]) for value in actuals)
    selected = (q < u < mp.sqrt(q) and p["direct_t_interval_lower"] < t < 1 and
                p["direct_unsquared_selector"] > 0 and m > 0 and left_squared < t)
    passed = (selected and normalization_error < tolerance and abs(right_squared - t) < tolerance and
              p["u_quartic_relative_residual"] < tolerance and
              p["t_quartic_relative_residual"] < tolerance and phase_invariance < tolerance and
              all(record["passed"] for record in phase_records))
    output = {"q_input": q_string,
              **{key: number(value) for key, value in p.items() if key != "root_iterations"},
              "root_iterations": p["root_iterations"],
              "geometric_normalization_error": number(normalization_error),
              "character_sum_squared_moduli": number(p["x"] ** 2 + y ** 2),
              "competitor_hankel_norm_abs_c_over_m": number(p["c"] / m),
              "R_squared": number(right_squared), "L_squared": number(left_squared),
              "right_active_slack_t_minus_R_squared": number(t - right_squared),
              "left_strict_slack_t_minus_L_squared": number(t - left_squared),
              "left_schur_slack_one_minus_m_minus_q_v_y_squared": number(left_slack_formula),
              "maximum_phase_norm_invariance_difference": number(phase_invariance),
              "phases": phase_records, "passed": bool(passed)}
    return output, largest_norm_difference, largest_stein_residual, phase_invariance


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default="f1_global_checks_results.json")
    parser.add_argument("--dps", type=int, default=90)
    args = parser.parse_args()
    if args.dps < 70:
        parser.error("--dps must be at least 70")
    mp.mp.dps = args.dps
    tolerance = mp.mpf("1e-60")
    results = [run_case(q, tolerance) for q in Q_CASES]
    cases = [result[0] for result in results]
    summary = {"all_passed": all(case["passed"] for case in cases),
               "q_cases": len(cases), "phase_cases_per_q": 3,
               "infinite_hankel_norm_evaluations": 3 * len(cases),
               "maximum_formula_vs_infinite_squared_norm_difference": number(max(r[1] for r in results)),
               "maximum_relative_stein_residual": number(max(r[2] for r in results)),
               "maximum_phase_invariance_difference": number(max(r[3] for r in results))}
    payload = {"status": "High-precision diagnostics; not interval certificates or a proof of globality.",
               "environment": {"python": platform.python_version(), "mpmath": mp.__version__},
               "decimal_precision": args.dps, "comparison_tolerance": number(tolerance),
               "method": "Monotone selected-root bisection plus two finite Stein equations for the infinite Hankel error norm; no word-horizon truncation.",
               "phase_convention": "lambda_b=y exp(i theta), c_theta=c exp(-i theta)",
               "summary": summary, "cases": cases}
    Path(args.output).write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    if not summary["all_passed"]:
        raise SystemExit("A diagnostic failed; inspect the JSON.")


if __name__ == "__main__":
    main()
