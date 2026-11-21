#!/usr/bin/env python3
"""Exact interior-word certificates for the failed row-contractive lift.

This verifies finite identities. The infinite-space conclusions are proved
analytically in the manuscript. Python standard library only.
"""
from fractions import Fraction
from itertools import product
from pathlib import Path
import json


def words(alphabet, minimum, maximum):
    return [w for length in range(minimum, maximum + 1)
            for w in product(alphabet, repeat=length)]


def backward(vector, letter):
    output = {}
    for word, value in vector.items():
        if len(word) >= 2 and word[-1] == letter:
            prefix = word[:-1]
            output[prefix] = output.get(prefix, Fraction(0)) + value
    return output


def squared_norm(vector):
    return sum((value * value for value in vector.values()), Fraction(0))


def main():
    checks = []

    def check(name, condition):
        checks.append({"check": name, "passed": bool(condition)})
        if not condition:
            raise AssertionError(name)

    for dimension in (2, 3, 4):
        alphabet = tuple(range(dimension))
        for prefix_length in (1, 2, 3):
            prefix = tuple(0 for _ in range(prefix_length))
            inputs = [{prefix + (letter,): Fraction(1)} for letter in alphabet]
            output = {}
            for letter, vector in zip(alphabet, inputs):
                shifted = backward(vector, letter)
                check(f"d={dimension},len={prefix_length},letter={letter}: output",
                      shifted == {prefix: Fraction(1)})
                for word, value in shifted.items():
                    output[word] = output.get(word, Fraction(0)) + value
            check(f"d={dimension},len={prefix_length}: input norm squared",
                  sum(squared_norm(v) for v in inputs) == dimension)
            check(f"d={dimension},len={prefix_length}: output norm squared",
                  squared_norm(output) == dimension * dimension)
            check(f"d={dimension},len={prefix_length}: strict failure",
                  squared_norm(output) > sum(squared_norm(v) for v in inputs))

        # Interior block: codomain words have lengths 1,...,L-1; domain has
        # lengths 1,...,L. Every right-created codomain word lies in the
        # domain, so the row Gramian is exactly d*I on this finite block.
        for maximum in (2, 3):
            domain = words(alphabet, 1, maximum)
            codomain = words(alphabet, 1, maximum - 1)
            index = {word: i for i, word in enumerate(codomain)}
            gram = [[0] * len(codomain) for _ in codomain]
            for letter in alphabet:
                for word in domain:
                    if len(word) >= 2 and word[-1] == letter:
                        i = index[word[:-1]]
                        gram[i][i] += 1
            expected = [[dimension if i == j else 0
                         for j in range(len(codomain))]
                        for i in range(len(codomain))]
            check(f"d={dimension},L={maximum}: exact interior row Gramian",
                  gram == expected)

    result = {
        "status": "PASS",
        "arithmetic": "exact integers and rational Fractions",
        "checks": len(checks),
        "scope": "Interior-word identities and finite row-Gramian blocks",
        "proof_assistant_formalization": False,
        "results": checks,
    }
    Path(__file__).with_name("contraction_exact_certificates.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in
                      ("status", "checks", "arithmetic")}))


if __name__ == "__main__":
    main()
