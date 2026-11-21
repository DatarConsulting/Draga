#!/usr/bin/env python3
"""Exact finite certificates for the binary fixed-rank divergent-ratio family.

All computations use rational Fractions. The universal lower bound is proved
in the accompanying note; this script verifies its explicit finite identities.
"""
from fractions import Fraction as F
from itertools import product
from math import comb
from pathlib import Path
from datetime import datetime, timezone
import json

from rational_matrices import (
    zeros, eye, transpose, add, scale, mul, outer, positive_map,
    characteristic_plus_coeffs,
)


CHECKS = []


def check(name, condition):
    CHECKS.append({"check": name, "passed": bool(condition)})
    if not condition:
        raise AssertionError(name)


def diag(values):
    out = zeros(len(values))
    for i, value in enumerate(values):
        out[i][i] = F(value)
    return out


def row_times(row, a):
    return [sum((row[i] * a[i][j] for i in range(len(row))), F(0))
            for j in range(len(row))]


def evaluate(a, b, alpha, beta, word):
    row = alpha[:]
    for letter in word:
        row = row_times(row, a if letter == "a" else b)
    return sum((x * y for x, y in zip(row, beta)), F(0))


def expected_value(word, k, r, q):
    if len(word) >= k and word.endswith("b" * k) and "b" not in word[:-k]:
        return q * r ** (len(word) - k)
    return F(0)


def binomial(n, j):
    return comb(n, j) if 0 <= j <= n else 0


def kronecker_sum(a, b):
    n = len(a)
    out = zeros(n * n)
    for i, j, p, q in product(range(n), repeat=4):
        out[i*n+j][p*n+q] = a[i][p] * a[j][q] + b[i][p] * b[j][q]
    return out


def characteristic_coeffs(a):
    return [F(1)] + [(-1)**j * c for j, c in enumerate(characteristic_plus_coeffs(a), start=1)]


def matrix_polynomial(a, coeffs_descending):
    out = zeros(len(a))
    for c in coeffs_descending:
        out = add(mul(out, a), scale(eye(len(a)), c))
    return out


def main():
    cases = []
    for k in (1, 2, 3, 4):
        central = comb(2*k, k)
        check(f"Vandermonde k={k}", sum(comb(k, j)**2 for j in range(k+1)) == central)
        for N in (2, 10, 100, 10000):
            r = F(N-1, N)
            q = 1 - r*r
            n = k+1
            a, b = zeros(n), zeros(n)
            a[0][0] = r
            for i in range(k):
                b[i][i+1] = F(1)
            alpha = [q] + [F(0)]*k
            beta = [F(0)]*k + [F(1)]
            p = scale(eye(n), q)
            gram_q = diag([1/q] + [F(1)]*k)
            label = f"k={k},N={N}"

            check(label + ": observability Stein", p == add(outer(alpha), positive_map([a,b], p, True)))
            check(label + ": reachability Stein", gram_q == add(outer(beta), positive_map([a,b], gram_q, False)))
            pq = mul(p, gram_q)
            check(label + ": exact squared singular values", pq == diag([F(1)] + [q]*k))
            expected_coeffs = [binomial(k,j)*q**j + binomial(k,j-1)*q**(j-1)
                               for j in range(1,n+1)]
            check(label + ": spectral polynomial", characteristic_plus_coeffs(pq) == expected_coeffs)
            cp = kronecker_sum(a,b)
            check(label + ": upper triangular CP matrix",
                  all(cp[i][j] == 0 for i in range(n*n) for j in range(i)))
            check(label + ": CP eigenvalues r^2 and zero",
                  [cp[i][i] for i in range(n*n)] == [r*r] + [F(0)]*(n*n-1) and r*r < 1)

            words = ["".join(w) for length in range(k+3) for w in product("ab", repeat=length)]
            check(label + ": exact word behavior",
                  all(evaluate(a,b,alpha,beta,w) == expected_value(w,k,r,q) for w in words))
            check(label + ": unit compressed target coefficient", q*q/(1-r*r)**2 == 1)
            check(label + ": rational ratio certificate", F(N*N, central*(2*N-1)) == 1/(central*q))

            # Proper probability automaton: all outgoing transitions + stopping mass sum to one.
            prob_b = [row[:] for row in b]
            prob_b[0][1] = 1-r
            prob_alpha = [F(1)] + [F(0)]*k
            check(label + ": proper probability rows",
                  all(sum(a[i])+sum(prob_b[i])+beta[i] == 1 for i in range(n)))
            check(label + ": probability target is scaled family",
                  all(evaluate(a,prob_b,prob_alpha,beta,w) == expected_value(w,k,r,q)/(1+r)
                      for w in words))

            cases.append({"k":k,"N":N,"r":str(r),"q":str(q),
                          "hankel_norm_squared":"1","small_singular_value_squared":str(q),
                          "structured_error_lower_bound_squared":str(F(1,central)),
                          "ratio_lower_bound_squared":str(1/(central*q))})

        # B=(1/2)I+(1/4)J is a contraction: ||B||<=3/4.
        contraction = scale(eye(k), F(1,2))
        for i in range(k-1):
            contraction[i][i+1] = F(1,4)
        descending = characteristic_coeffs(contraction)
        check(f"Cayley-Hamilton contraction k={k}", matrix_polynomial(contraction,descending) == zeros(k))
        coeffs = list(reversed(descending))
        check(f"coefficient norm bound k={k}", sum(c*c for c in coeffs) <= central)
        gamma = [F(1)] + [F(0)]*(k-1)
        nu = [F(0)]*(k-1) + [F(1)]
        r = F(99,100)
        q = 1-r*r
        power = eye(k)
        etas = []
        for j in range(k+1):
            moment = sum((x*y for x,y in zip(row_times(gamma,power),nu)),F(0))
            etas.append(q * moment / (1-r/F(4))**2)
            power = mul(power,contraction)
        check(f"weighted compressed recurrence k={k}", sum(c*e for c,e in zip(coeffs,etas)) == 0)

        for M in (1, 10, 100):
            N = 2*central*M*M+1
            check(f"arbitrarily large prescribed ratio k={k},M={M}",
                  F(N*N,central*(2*N-1)) > M*M)

        # Finite-support corollary: exact tail mass and strictly positive lower endpoint.
        for N in (10,20):
            r=F(N-1,N); q=1-r*r; L=N*N
            eta=(1+r)*r**(L+1)
            check(f"finite-support geometric tail k={k},N={N}",
                  q*r**(L+1)/(1-r) == eta)
            check(f"finite-support positive lower endpoint k={k},N={N}", eta*eta < F(1,central))

    result={"status":"PASS","generated_at_utc":datetime.now(timezone.utc).isoformat(),
            "arithmetic":"exact rational Fractions","checks":len(CHECKS),
            "parameter_cases":len(cases),"scope":"Finite identities; universal bound is proved analytically",
            "proof_assistant_formalization":False,"global_optimizer_executed":False,
            "cases":cases,"results":CHECKS}
    Path(__file__).with_name('binary_exact_certificates.json').write_text(
        json.dumps(result,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(json.dumps({key:result[key] for key in ('status','checks','parameter_cases','arithmetic')}))


if __name__=='__main__':
    main()
