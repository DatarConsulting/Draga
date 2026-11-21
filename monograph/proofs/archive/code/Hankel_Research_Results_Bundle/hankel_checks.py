#!/usr/bin/env python3
"""Floating-point diagnostics for the word-Hankel structural-oracle manuscript.

Dependencies: numpy, scipy. Run: python hankel_checks.py --output hankel_checks_results.json
No floating-point computation in this file is offered as a proof or interval certificate.

Sequential convention: x_{wa}=A_a x_w and f(w)=c A_reverse(w) b.
The transpose word-Hankel factorization has reachability/observability Gramians
M=bb*+sum A_a M A_a*, N=c*c+sum A_a* N A_a. Reversal/permutation
of word indices does not change the operator norm; the standard error-Gramian
formula consequently computes the same word-Hankel norm in this convention.
"""
from __future__ import annotations
import argparse
import json
import platform
from pathlib import Path
import numpy as np
import scipy
from scipy.linalg import block_diag, solve


ATOL = 5e-9


def adj(a):
    return a.conj().T


def herm(a):
    return (a + adj(a)) / 2


def psd_sqrt(a):
    vals, vecs = np.linalg.eigh(herm(a))
    scale = max(1.0, float(np.max(np.abs(vals))))
    if float(vals.min()) < -1e-8 * scale:
        raise AssertionError(f"PSD failure beyond roundoff: {vals.min()}")
    return (vecs * np.sqrt(np.maximum(vals, 0.0))) @ adj(vecs)


def stein(letters, source, side="row"):
    """Solve a finite Stein equation by a Kronecker linear system, no word enumeration."""
    n = source.shape[0]
    dtype = np.result_type(source.dtype, *[a.dtype for a in letters])
    kernel = np.zeros((n*n, n*n), dtype=dtype)
    if side == "row":
        for a in letters:
            kernel += np.kron(a.conj(), a)
    elif side == "column":
        for a in letters:
            kernel += np.kron(a.T, adj(a))
    else:
        raise ValueError(side)
    vec = solve(np.eye(n*n, dtype=dtype)-kernel,
                source.reshape(-1, order="F"), assume_a="gen")
    result = herm(vec.reshape((n,n), order="F"))
    if side == "row":
        residual = result-source-sum(a @ result @ adj(a) for a in letters)
    else:
        residual = result-source-sum(adj(a) @ result @ a for a in letters)
    relative = float(np.linalg.norm(residual, "fro") /
                     max(1., np.linalg.norm(result, "fro"), np.linalg.norm(source, "fro")))
    return result, relative


def hankel_norm(letters, b, c):
    """Infinite operator norm from finite exact-in-principle error Gramians."""
    m, mr = stein(letters, b @ adj(b), "row")
    n, nr = stein(letters, adj(c) @ c, "column")
    ns = psd_sqrt(n)
    val = float(np.linalg.eigvalsh(herm(ns @ m @ ns)).max())
    return np.sqrt(max(0., val)), max(mr, nr)


def product_sum_bound(letters, side, target=0.2, max_horizon=2048):
    """Return sum_t r_t or c_t upper bound using one finite CP block test.

    If q_L<1, submultiplicativity gives sum_t q_t <=
    (sum_{t=0}^{L-1} q_t)/(1-q_L). q_t is evaluated by CP iteration.
    """
    k = letters[0].shape[0]
    power = np.eye(k, dtype=letters[0].dtype)
    partial = 0.0
    for horizon in range(1, max_horizon+1):
        qt = np.sqrt(max(0., float(np.linalg.eigvalsh(herm(power)).max())))
        partial += qt
        if side == "row":
            power = herm(sum(a @ power @ adj(a) for a in letters))
        else:
            power = herm(sum(adj(a) @ power @ a for a in letters))
        qL = np.sqrt(max(0., float(np.linalg.eigvalsh(power).max())))
        if qL <= target:
            return {"upper_bound": float(partial/(1-qL)),
                    "horizon": horizon, "block_norm": float(qL),
                    "partial_sum": float(partial)}
    raise AssertionError(f"No strict finite product certificate by horizon {max_horizon}")


def certificate(letters, b, c, k):
    n = letters[0].shape[0]
    d = len(letters)
    m, mr = stein(letters, b @ adj(b), "row")
    obs, nr = stein(letters, adj(c) @ c, "column")
    oval, ovec = np.linalg.eigh(herm(obs))
    if oval.min() <= 0:
        raise AssertionError("The test input must be observable/minimal.")
    root = (ovec*np.sqrt(oval)) @ adj(ovec)
    invroot = (ovec*(1/np.sqrt(oval))) @ adj(ovec)
    canonical = [root @ a @ invroot for a in letters]
    bc = root @ b
    cc = c @ invroot
    reach = herm(root @ m @ root)
    eig, u = np.linalg.eigh(reach)
    order = np.argsort(eig)[::-1]
    eig = np.maximum(eig[order], 0.0)
    u = u[:, order]
    v, w = u[:, :k], u[:, k:]
    singulars = np.sqrt(eig)
    epsilon = float(singulars[k])
    retained = singulars[:k]
    C = np.diag(retained)
    Ci = np.diag(1/retained)
    reduced = [adj(v) @ a @ v for a in canonical]
    weighted = [Ci @ a @ C for a in reduced]
    b_reduced = adj(v) @ bc
    c_reduced = cc @ v
    eta_blocks = [adj(v) @ a @ w for a in canonical]
    tau_blocks = [adj(w) @ a @ v @ C for a in canonical]
    eta = np.sqrt(max(0., float(np.linalg.eigvalsh(herm(
        sum(e @ adj(e) for e in eta_blocks))).max())))
    tau = np.sqrt(max(0., float(np.linalg.eigvalsh(herm(
        sum(adj(t) @ t for t in tau_blocks))).max())))
    kr = product_sum_bound(reduced, "row")
    kc = product_sum_bound(weighted, "column")
    D = epsilon*(1+eta*kr["upper_bound"])+tau*kc["upper_bound"]
    C_oracle = 1+np.sqrt(d)*kr["upper_bound"]+(np.sqrt(d)+1)*kc["upper_bound"]
    error_letters = [block_diag(a, ar) for a, ar in zip(letters, reduced)]
    error_b = np.concatenate([b, b_reduced], axis=0)
    error_c = np.concatenate([c, -c_reduced], axis=1)
    actual, er = hankel_norm(error_letters, error_b, error_c)
    passivity_residual = float(np.linalg.norm(
        np.eye(n)-adj(cc) @ cc-sum(adj(a) @ a for a in canonical), "fro"))
    tau_upper = (np.sqrt(d)+1)*epsilon
    tolerance = ATOL*max(1., D, actual, tau_upper)
    passed = bool(actual <= D+tolerance and tau <= tau_upper+tolerance and
                  D <= C_oracle*epsilon+tolerance)
    return {"n": n, "k": k, "d": d,
            "singular_values": singulars.tolist(), "epsilon": epsilon,
            "eta": float(eta), "tau": float(tau),
            "tau_upper": float(tau_upper),
            "tau_over_epsilon": float(tau/epsilon),
            "Kr": kr, "Kc": kc, "D": float(D),
            "actual_infinite_hankel_error": float(actual),
            "actual_over_D": float(actual/D),
            "oracle_constant_upper": float(C_oracle),
            "gramian_relative_residual": max(mr, nr, er),
            "canonical_passivity_residual_fro": passivity_residual,
            "observability_condition_number": float(oval.max()/oval.min()),
            "passed": passed}


def random_case(n, d, field, seed):
    rng = np.random.default_rng(seed)
    def draw(shape):
        out = rng.normal(size=shape)
        if field == "complex":
            out = (out+1j*rng.normal(size=shape))/np.sqrt(2)
        return out
    letters = [draw((n,n))/np.sqrt(n) for _ in range(d)]
    rownorm = float(np.linalg.eigvalsh(herm(sum(a @ adj(a) for a in letters))).max())
    colnorm = float(np.linalg.eigvalsh(herm(sum(adj(a) @ a for a in letters))).max())
    # Both one-step CP norms are below 0.65^2; canonical coordinates can be nonnormal.
    letters = [a*(0.65/np.sqrt(max(rownorm, colnorm))) for a in letters]
    b, c = draw((n,1)), draw((1,n))
    b /= np.linalg.norm(b)
    c /= np.linalg.norm(c)
    return letters, b, c


def binary_section(r, k, horizon):
    """Enumerate only potentially nonzero rows/columns, not the full free monoid.

    This is the actual P_{<=L} H P_{<=L} finite section with all zero rows/columns
    omitted. Compare its singular values with independently derived finite formulas.
    """
    rows = ["a"*i for i in range(horizon+1)]
    rows += ["a"*i+"b"*j for j in range(1,k+1)
             for i in range(horizon-j+1)]
    cols = ["a"*i+"b"*k for i in range(horizon-k+1)]
    cols += ["b"*j for j in range(k)]
    q = 1-r*r
    def coefficient(word):
        if word.endswith("b"*k):
            initial = word[:-k]
            if all(letter == "a" for letter in initial):
                return q*r**len(initial)
        return 0.0
    matrix = np.array([[coefficient(u+v) for v in cols] for u in rows])
    measured = np.linalg.svd(matrix, compute_uv=False)[:k+1]
    finite_formula = [np.sqrt((1-r**(2*(horizon+1)))*
                             (1-r**(2*(horizon-k+1))))]
    finite_formula += [np.sqrt(q*(1-r**(2*(horizon-k+s+1)))) for s in range(k)]
    finite_formula = np.sort(finite_formula)[::-1]
    full = np.array([1.]+[np.sqrt(q)]*k)
    discrepancy = float(np.max(np.abs(measured-finite_formula)))
    return {"r": r, "k": k, "horizon": horizon,
            "nonzero_section_shape": list(matrix.shape),
            "measured_nonzero_singulars": measured.tolist(),
            "finite_section_formula": finite_formula.tolist(),
            "full_operator_singulars": full.tolist(),
            "finite_section_formula_max_abs_error": discrepancy,
            "full_minus_finite_singulars": (full-measured).tolist(),
            "passed": bool(discrepancy < 2e-12 and np.all(measured <= full+2e-12))}


def sharp_one_state(h):
    """Rank-three infinite-norm check of the explicit nearoptimal competitor.

    Basis vectors are v, e_epsilon, p. Their Gram matrix and the coefficient
    matrix of (H-G)^*(H-G) are exact analytic formulas evaluated in floats.
    """
    s = h/(1-h)
    A = 1-h*h
    gram = np.array([[1.,0.,h],[0.,1.,np.sqrt(h)],[h,np.sqrt(h),1.]])
    matrix = np.array([[1.,0.,-h],[0.,s*s,-s*h**1.5],
                       [-h,-s*h**1.5,h]])
    root = psd_sqrt(gram)
    actual_squared = float(np.linalg.eigvalsh(herm(root @ matrix @ root)).max())
    lower = A+h**5-7*h**6
    upper = A+h**5+4*h**6
    return {"h":h,"q":s*s,"explicit_error_squared":actual_squared,
            "universal_optimum_lower_bound":lower,
            "explicit_competitor_upper_bound":upper,
            "passed":bool(lower-5e-15 <= actual_squared <= upper+5e-15)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default="hankel_checks_results.json")
    args = parser.parse_args()
    random_results = []
    # All retained dimensions in n=3,...,8, over both real and complex scalars.
    # Two independent seed families give 108 separate reduction certificates.
    for family in (0, 1):
        for field_index, field in enumerate(("real", "complex")):
            for n in range(3,9):
                d = 2 if family == 0 else 3
                seed = 20261024+10000*family+1000*field_index+17*n
                letters, b, c = random_case(n, d, field, seed)
                for k in range(1,n):
                    record = certificate(letters,b,c,k)
                    record.update({"field":field,"seed":seed,"family":family})
                    random_results.append(record)
    binary = [binary_section(r,k,L) for r in (0.2,0.75,0.97,0.995)
              for k in (1,2,3) for L in (32,96)]
    # Three states; reading a then b produces coefficient 0.1, and the empty
    # coefficient is 1. All other coefficients are zero.
    aa = np.zeros((3,3)); aa[1,0] = 1
    ab = np.zeros((3,3)); ab[2,1] = 1
    b = np.array([[1.],[0.],[0.]])
    c = np.array([[1.,0.,0.1]])
    impulse_ab = certificate([aa,ab],b,c,1)
    impulse_ab["series"] = "f(epsilon)=1, f(ab)=0.1, all other coefficients zero"
    sharp = [sharp_one_state(1/j) for j in (64,80,100,128,200)]
    all_records = random_results+binary+[impulse_ab]+sharp
    summary = {
        "all_passed": all(record["passed"] for record in all_records),
        "random_reductions":len(random_results),"binary_sections":len(binary),
        "sharp_one_state_checks":len(sharp),
        "random_fields":["real","complex"],"random_n":[3,4,5,6,7,8],
        "random_d":[2,3], "random_k":"all integers 1 through n-1",
        "max_actual_over_D":max(r["actual_over_D"] for r in random_results),
        "max_tau_over_epsilon":max(r["tau_over_epsilon"] for r in random_results),
        "max_gramian_relative_residual":max(r["gramian_relative_residual"] for r in random_results),
        "max_canonical_passivity_residual_fro":max(r["canonical_passivity_residual_fro"] for r in random_results),
        "max_binary_formula_abs_error":max(r["finite_section_formula_max_abs_error"] for r in binary),
        "impulse_ab_epsilon":impulse_ab["epsilon"],
        "impulse_ab_D":impulse_ab["D"],
        "impulse_ab_actual_error":impulse_ab["actual_infinite_hankel_error"]}
    results = {
        "status":"floating-point diagnostics only; not proofs, interval bounds, or a novelty certificate",
        "environment":{"python":platform.python_version(),"numpy":np.__version__,"scipy":scipy.__version__},
        "tolerance":ATOL,"seed_construction":"20261024+10000*family+1000*field_index+17*n",
        "methods":{"infinite_norm":"finite Stein Gramians, never a horizon truncation",
                   "stability_constants":"finite CP-map product block test; target block norm <= 0.2",
                   "binary_sections":"all words of row/column length <= L, zero rows/columns omitted"},
        "summary":summary,"random_cases":random_results,"binary_cases":binary,
        "impulse_ab":impulse_ab,"sharp_one_state":sharp}
    Path(args.output).write_text(json.dumps(results,indent=2)+"\n")
    print(json.dumps(summary,indent=2))
    if not summary["all_passed"]:
        raise SystemExit("A floating-point diagnostic failed; inspect the JSON.")


if __name__ == "__main__":
    main()
