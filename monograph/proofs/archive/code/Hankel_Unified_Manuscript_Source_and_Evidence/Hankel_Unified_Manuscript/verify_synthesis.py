#!/usr/bin/env python3
"""Independent checks of displayed formulas; standard library only.

No historical C, model training, or GDN experiment is rerun by this script.
Rational checks are exact. The quartic branch checks are explicitly numerical.
"""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import json
import math
import platform
import sys
import time

ROOT = Path(__file__).resolve().parent


def transpose(a):
    return list(map(list, zip(*a)))


def mm(a, b):
    return [[sum((x*y for x, y in zip(row, col)), F(0))
             for col in zip(*b)] for row in a]


def plus(*matrices):
    return [[sum(mat[i][j] for mat in matrices)
             for j in range(len(matrices[0][0]))]
            for i in range(len(matrices[0]))]


def gramian_check():
    rho = F(1, 2)
    transitions = [[[rho, rho*c/2], [F(0), rho/2]] for c in (0, 1)]
    p = [[F(4, 21), F(2, 21)], [F(2, 21), F(8, 7)]]
    q = [[F(2), F(1, 3)], [F(1, 3), F(4, 21)]]
    bb = [[F(0), F(0)], [F(0), F(1)]]
    cc = [[F(1), F(0)], [F(0), F(0)]]
    assert p == plus(bb, *(mm(mm(a, p), transpose(a)) for a in transitions))
    assert q == plus(cc, *(mm(mm(transpose(a), q), a) for a in transitions))
    pq = mm(p, q)
    hs2 = sum(pq[i][i] for i in range(2))
    assert hs2 == F(292, 441) < F(25, 36)
    return {"stein_identities": 2, "hilbert_schmidt_squared": str(hs2)}


def exact_graph(a, m):
    assert a >= 1 and m >= a+2
    rho = F(1, 2**a)
    L = (a+1)*m+2
    eps = F(1, 2**L)
    T = (L-2)//a+1
    widths = [min(t, L-1-a*t) for t in range(T)]
    offsets = []
    nstates = 0
    for n in widths:
        offsets.append(nstates)
        nstates += 2**n
    absorbing = nstates
    nstates += 1
    assert 2**(m+1)+1 <= nstates < 2**(m+2)
    if a == 1:
        assert nstates == 2**(m+2)-2

    def cell(t, b):
        n = widths[t]
        lo_k = b*2**(t-n)
        hi_k = (b+1)*2**(t-n)-1
        return F(lo_k, 2**((a+1)*t)), F(hi_k, 2**((a+1)*t))

    def decode(state_id):
        if state_id == absorbing:
            return None
        for t, offset in enumerate(offsets):
            if state_id < offset+2**widths[t]:
                return t, state_id-offset
        raise AssertionError("invalid ID")

    transition_count = 0
    output_checks = 0
    for t, n in enumerate(widths):
        for b in range(2**n):
            sid = offsets[t]+b
            assert decode(sid) == (t, b)
            lo, hi = cell(t, b)
            y = F(b, 2**((a+1)*t)) if t <= m else (2*b+1)*eps
            assert abs(lo-y) <= eps and abs(hi-y) <= eps
            output_checks += 2
            for c in (0, 1):
                if t+1 == T:
                    new_id = absorbing
                    dest_lo, dest_hi = F(0), 2*eps
                else:
                    nn = widths[t+1]
                    bnew = c*2**(nn-1) + (b >> (n-nn+1))
                    new_id = offsets[t+1]+bnew
                    assert decode(new_id) == (t+1, bnew)
                    dest_lo, dest_hi = cell(t+1, bnew)
                new_lo = rho*(lo+c*rho**t)/2
                new_hi = rho*(hi+c*rho**t)/2
                assert dest_lo <= new_lo <= new_hi <= dest_hi
                assert 0 <= new_id < nstates
                transition_count += 1
    assert rho**T <= 2*eps
    output_checks += 2
    transition_count += 2  # two absorbing transitions
    return {"a": a, "m": m, "L": L, "states": nstates,
            "transitions": transition_count, "endpoint_checks": output_checks,
            "persistent_bits": m+2}


def packing(a, m):
    rho = F(1, 2**a)
    eps = F(1, 2**((a+1)*m+2))
    witness = [(rho**t, F(k, 2**((a+1)*t)))
               for t in range(m+1) for k in range(2**t)]
    witness += [(rho**(m+1), F(0)),
                (rho**(m+1), F(2**(m+1)-1, 2**((a+1)*(m+1))))]
    assert len(witness) == 2**(m+1)+1
    smallest = None
    pair_count = 0
    for (z, x), (zz, xx) in combinations(witness, 2):
        dist = max(abs(x-xx), rho*abs(x-xx+z-zz)/2)
        assert dist >= F(16, 5)*eps
        assert dist-2*(eps/4) > 2*(eps+eps/4)
        smallest = dist if smallest is None else min(smallest, dist)
        pair_count += 1
    return {"a": a, "m": m, "prefixes": len(witness),
            "pairs": pair_count, "minimum_distance_over_epsilon": str(smallest/eps)}


def shell_checks():
    rows = []
    for d in (2, 3, 4):
        for T in (2, 3, 5, 8):
            delta = 0.125
            b = delta*d**(-T/2)
            root = math.hypot(1-delta, 2*b)
            leading = ((1-delta)+root)/2
            trailing = 2*b*b/(root+(1-delta))
            assert trailing < b < leading
            values = sorted([leading, trailing]+[b]*(T-1), reverse=True)
            assert values[1] == b
            M = delta*T
            g = (1+M)*math.log1p(M)-M*math.log(M)
            h = -delta*math.log(delta)-(1-delta)*math.log1p(-delta)
            aopt = 1/(1+M)
            bopt = M/(d*(1+M))
            direct = ((1-delta)*math.log((1-delta)/aopt)
                      + delta*(math.log(delta)-T*math.log(d)
                               -math.log(aopt)-T*math.log(bopt)))
            assert abs(direct-(g-h)) < 1e-12
            rows.append({"alphabet": d, "length": T, "sigma2": b,
                         "kl_formula_absolute_error": abs(direct-(g-h))})
    return rows


def quartic_checks():
    rows = []
    for q in (0.001, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99):
        def poly(u):
            return q*q*(1-u*u)*(u-q)**2-(1-q)**3*u*u*(q-u*u)
        lo, hi = q, math.sqrt(q)
        for _ in range(100):
            mid = (lo+hi)/2
            if poly(mid) > 0:
                hi = mid
            else:
                lo = mid
        u = (lo+hi)/2
        t = (q*q+(1-2*q)*u*u)/(q*(1-u*u))
        aq = (1-q)**2*t*t+(q*q+2*q-1)*(t-q)
        residual = aq*aq-4*q**3*(t-q)*(1-2*q+q*t)
        alpha = (1-u)/math.sqrt(1-q)
        beta = math.sqrt(u*(1-u)/(1-q))
        omega = (u-q)/(1-q)
        assert q < u < math.sqrt(q)
        assert 2*q/(1+q) < t < 1
        assert aq > 0 and abs(residual) < 1e-10
        assert abs(alpha*alpha+beta*beta-(1-omega)) < 1e-12
        assert alpha*alpha+beta*beta < 1
        rows.append({"q": q, "u": u, "squared_error_branch": t,
                     "quartic_absolute_residual": abs(residual)})
    return rows


def main():
    started = time.perf_counter()
    result = {"scope": "Independent synthesis checks; no historical C or LLM rerun",
              "python": sys.version, "platform": platform.platform(),
              "gramian": gramian_check()}
    result["exact_graphs"] = [exact_graph(a, m) for a in range(1, 5)
                              for m in range(a+2, 10)]
    result["exact_packings"] = [packing(a, m) for a in range(1, 4)
                                for m in range(a+2, 7)]
    result["numerical_shell_checks"] = shell_checks()
    result["numerical_quartic_checks"] = quartic_checks()
    result["summary"] = {
        "graph_configurations": len(result["exact_graphs"]),
        "states": sum(x["states"] for x in result["exact_graphs"]),
        "transitions": sum(x["transitions"] for x in result["exact_graphs"]),
        "packing_configurations": len(result["exact_packings"]),
        "witness_pairs": sum(x["pairs"] for x in result["exact_packings"]),
        "shell_cases": len(result["numerical_shell_checks"]),
        "quartic_cases": len(result["numerical_quartic_checks"]),
        "all_assertions_passed": True,
    }
    result["elapsed_seconds"] = time.perf_counter()-started
    (ROOT/"verification_results.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result["summary"], indent=2))


if __name__ == "__main__":
    main()
