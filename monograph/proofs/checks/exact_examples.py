"""Independent finite witnesses for three arguments in the monograph.

Every comparison uses integers or fractions. Each test specifies a finite family.
For the quantified claims, consult the cited proof and the Lean formalization.
"""
from fractions import Fraction as F
from itertools import product
import json


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def block_budget():
    """§2.7.6. Enumerate all two-state, two-column tables with row mass at most D.
    The two one-hot rows attain equality, detecting an erroneous state-free bound.
    """
    count = 0
    for d in range(1, 11):
        rows = [(a, b) for a in range(d+1) for b in range(d+1-a)]
        for (a, b), (c, e) in product(rows, repeat=2):
            require(max(a,c) + max(b,e) <= 2*d, "Block budget")
            count += 1
        require(sum(max(x, y) for x, y in zip((d,0),(0,d))) == 2*d, "Sharpness")
    return count


def stopping_intervals(threshold):
    """§3.13.1. Descend the actual affine images, including the overshooting leaves.
    Zero has cost two and maps x to x/4. One has cost one and maps x to (1+x)/2.
    """
    stack = [(0, F(0), F(1))]
    while stack:
        cost, a, width = stack.pop()
        if cost >= threshold:
            yield a, a + width
        else:
            stack.append((cost+1, a+width/2, width/2))
            stack.append((cost+2, a, width/4))


def fibonacci_witness():
    """Check endpoint spacing with a closed-interval convention and exact dyadics."""
    f = [0, 1]
    for _ in range(129):
        f.append(f[-1]+f[-2])
    total_leaves = 0
    for n in range(1, 17):
        theta = F(1,2**n)
        intervals = list(stopping_intervals(n))
        require(len(intervals) == f[n+2], "Terminal count")
        require(all(b-a in (theta, theta/2) for a,b in intervals), "Terminal width")
        endpoints = [b for _,b in intervals]
        require(all(y-x >= 3*theta/2 for x,y in zip(endpoints,endpoints[1:])), "Packing")
        total_leaves += len(intervals)
    bits = (f[130]-1).bit_length()
    require(bits == 90, "Threshold 128 bit budget")
    return {"enumerated_leaves": total_leaves, "F130": str(f[130]), "bits_at_128": bits}


def active_maps():
    """§3.15.9. For a finite (k,v), test code bijectivity, exact bit length and
    every update on every map. The digit zero represents an absent key. Separating
    key queries give the matching lower bound when all maps are reachable.
    """
    tested, transitions = 0, 0
    for k in range(1, 6):
        for v in range(1, 5):
            radix = v+1
            place = [radix**i for i in range(k)]
            size = radix**k
            codes = set()
            for m in product(range(radix),repeat=k):
                code = sum(t*p for t,p in zip(m,place))
                require(tuple((code//p)%radix for p in place)==m, "Decoding")
                require(code not in codes, "Distinct maps collided")
                codes.add(code)
                for i in range(k):
                    for t in range(radix):
                        new = m[:i] + (t,) + m[i+1:]
                        require(code+(t-m[i])*place[i]==sum(x*p for x,p in zip(new,place)), "Update")
                        transitions += 1
            bits=(size-1).bit_length()
            require(len(codes)==size and max(codes)==size-1 and 2**(bits-1)<size<=2**bits, "Bit width")
            tested += size
    return {"maps":tested,"individual_updates":transitions}


if __name__ == "__main__":
    print(json.dumps(dict(substochastic_tables=block_budget(),stopping=fibonacci_witness(),
                          active_maps=active_maps()),indent=2))
