"""Compare the actual Lean ordinal implementation with independent exact interval arithmetic.

For thresholds 1 through 12, subdivide real intervals using Fraction and number
them left to right. Check every one-letter transition and every midpoint. For
larger thresholds, represent an address by its newest-first binary word and
use a small suffix-count table to identify the ordinal. Test the output error
at both endpoints of the initial cell and at every subsequent step.
"""
from dataclasses import dataclass
from fractions import Fraction
import json
from pathlib import Path
import random
import subprocess

PROJECT = Path(__file__).resolve().parents[1] / "lean-real"


@dataclass(frozen=True)
class Cell:
    word: str
    left: Fraction
    width: Fraction


def small_cells(n, prefix="", left=Fraction(0), width=Fraction(1)):
    if n <= 0:
        yield Cell(prefix, left, width)
    else:
        yield from small_cells(n - 2, prefix + "0", left, width / 4)
        yield from small_cells(n - 1, prefix + "1", left + width / 2, width / 2)


def trim(word, threshold):
    total = 0
    for i, bit in enumerate(word, start=1):
        total += 2 if bit == "0" else 1
        if total >= threshold:
            return word[:i]
    raise ValueError("Reference address did not reach the stopping threshold")


def interval(word):
    left, width = Fraction(0), Fraction(1)
    for bit in word:
        if bit == "1":
            left += width/2
            width /= 2
        else:
            width /= 4
    return Cell(word, left, width)


def counts(n):
    table = [1, 2]
    while len(table) <= n:
        table.append(table[-2] + table[-1])
    return table


def ordinal(word, n, sizes):
    index = 0
    for bit in word:
        if bit == "1":
            index += sizes[max(0, n-2)]
            n -= 1
        else:
            n -= 2
    return index


def step(x, bit):
    return (1+x)/2 if bit == "1" else x/4


def compare():
    inputs, expected = [], []
    small_edges = large_steps = 0

    def add(n, word, letters, ids=None):
        nonlocal large_steps
        sizes = counts(n)
        def number(w):
            return ids[w] if ids is not None else ordinal(w, n, sizes)
        inputs.append(f"{n} {number(word)} {letters or '-'}")
        start = interval(word)
        low, high = start.left, start.left + start.width
        error = Fraction(1, 2**(n+1))
        for bit in [None, *letters]:
            if bit is not None:
                word = trim(bit+word,n)
                low, high = step(low,bit), step(high,bit)
                large_steps += int(ids is None)
            cell = interval(word)
            midpoint = cell.left + cell.width/2
            if max(abs(midpoint-low),abs(midpoint-high)) > error:
                raise ValueError("Reference interval violated the contract")
            expected.append((number(word),midpoint.numerator,midpoint.denominator))

    for n in range(1,13):
        cells = list(small_cells(n))
        ids = {cell.word:i for i,cell in enumerate(cells)}
        if len(ids) != counts(n)[n]:
            raise ValueError("Duplicate address or incorrect count")
        for cell in cells:
            for letter in "01":
                add(n,cell.word,letter,ids)
                small_edges += 1
    rng = random.Random(31302)
    for n in (32, 64, 128):
        seeds = ["0"*n,"1"*n,"01"*n,"10"*n]
        seeds += ["".join(rng.choice("01") for _ in range(n)) for _ in range(4)]
        for seed in seeds:
            address = trim(seed,n)
            for letters in ("", "0"*64, "1"*64, "01"*48,
                            "".join(rng.choice("01") for _ in range(256))):
                add(n,address,letters)
    run = subprocess.run(["lake","env","lean","--run","Trace.lean"],cwd=PROJECT,
                         text=True,input="\n".join(inputs)+"\n",capture_output=True,check=True)
    actual = [tuple(map(int,line.split())) for line in run.stdout.splitlines()]
    if len(actual) != len(expected):
        raise ValueError(f"Output count mismatch: {len(actual)} vs {len(expected)}")
    for row,(a,b) in enumerate(zip(actual,expected),start=1):
        if a != b:
            raise ValueError(f"Ordinal or exact rational mismatch at output {row}: {a} vs {b}")
    return dict(exhaustive_thresholds=[1,12],exhaustive_edges=small_edges,
                large_thresholds=[32,64,128],large_trace_transitions=large_steps,
                traces=len(inputs),exact_readouts=len(expected),
                states_at_128=str(counts(128)[128]),bits_at_128=(counts(128)[128]-1).bit_length())


if __name__ == "__main__":
    print(json.dumps(compare(),indent=2))
