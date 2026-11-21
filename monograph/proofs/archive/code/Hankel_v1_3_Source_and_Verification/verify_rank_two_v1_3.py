#!/usr/bin/env python3
"""Exact algebra and finite checks for the version-1.3 additions.

Requires Python 3 and SymPy. These checks do not replace the proofs.
"""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import hashlib
import json
import sympy as s

R = s.Rational
I = s.eye(2)
one = s.ones(2, 1)
swap = s.Matrix([[0, 1], [1, 0]])
A0 = s.Matrix([[2, 1], [1, 1]]) / 10
A1 = swap * A0 * swap
alpha = s.Matrix([[R(1, 2), R(1, 2)]])
stop = one / 2
r = (s.sqrt(5) - 1) / 2
gamma = (7 - 3 * s.sqrt(5)) / 2
h = 9 - 4 * s.sqrt(5)
S = s.Matrix([[r, 1-r], [1-r, r]])
G0 = s.Matrix([[(3+s.sqrt(5))/20, 0], [R(1,5), (3-s.sqrt(5))/20]])
G1 = swap * G0 * swap

def zero(matrix):
    return all(s.simplify(x) == 0 for x in matrix)

assert zero(S*A0-G0*S) and zero(S*A1-G1*S)
assert zero((G0+G1)*one-one/2)
assert (alpha*(I-A0-A1).inv()*stop)[0] == 1
assert s.Matrix.vstack(alpha, alpha*A0).det() == -R(1,40)
assert s.Matrix.hstack(stop, A0*stop).det() == -R(1,40)
assert (A0+A1).eigenvals() == {R(1,2): 1, R(1,10): 1}
assert zero(s.Matrix([[1,-1]])*A0-s.Matrix([[R(1,10),0]]))
assert zero(s.Matrix([[1,-1]])*A1-s.Matrix([[0,-R(1,10)]]))

t = s.symbols("t")
F0 = (1+t)/(2+t)
second = -(1+s.sqrt(5))/2
assert s.simplify((F0-r)/(F0-second)-gamma*(t-r)/(t-second)) == 0
assert s.simplify((R(1,2)-r)/(R(1,2)-second)+h) == 0
assert zero(
    s.Matrix([[R(8,13),R(5,13)]])
    +R(2,13)*alpha
    -R(15,13)*s.Matrix([[R(3,5),R(2,5)]]))
z = s.symbols("z", positive=True)
a_expr = r-s.sqrt(5)*z/(1+z)
assert s.simplify((1-a_expr-a_expr**2)/25-z/(5*(1+z)**2)) == 0

def clipped_matrices(a):
    B0 = s.Matrix([
        [(2+a)/10, 0],
        [(-a*a+3*a-1)/(10*(2*a-1)),
         (-a*a+4*a-2)/(10*(2*a-1))]])
    return B0, swap*B0*swap

golden_rows = []
a = R(1,2)
for L in range(1,13):
    a = s.cancel((1+a)/(2+a))
    B0, B1 = clipped_matrices(a)
    Q = B0+B1
    assert all(x >= 0 for x in B0) and all(x >= 0 for x in B1)
    assert zero(Q*one-one/2) and zero(one.T*Q-one.T/2)
    x = alpha*(I-Q).inv()
    assert x == s.Matrix([[1,1]])
    dictionary = s.Matrix([[a,1-a],[1-a,a]])
    k = (1-a-a*a)/10
    assert zero(dictionary*A0-B0*dictionary-s.Matrix([[k,-k],[0,0]]))
    assert zero(dictionary*A1-B1*dictionary-s.Matrix([[0,0],[-k,k]]))
    J = 2*k/R(5)
    assert J == (1-a-a*a)/25 > 0
    beta = s.Matrix([[R(1,4),R(1,4)]])
    v = beta*(I-Q).inv()
    assert v == s.Matrix([[R(1,2),R(1,2)]])
    assert zero(v-beta-v*Q)
    assert not zero(x-beta-x*Q)
    golden_rows.append({"depth":L,"upper_coordinate":str(a),
                        "defect":str(J),"tv_certificate":str(J/2),
                        "anchored_tv_certificate":str(J/4)})

# Compare endpoint iteration with exhaustive actual-prefix coordinates.
Ap = A0
Am = s.Matrix([[1,2],[2,1]])/20
tau_m = (I-Ap-Am)*one
assert tau_m == s.Matrix([R(11,20),R(13,20)])
assert (alpha*(I-Ap-Am).inv()*tau_m)[0] == 1
assert s.Matrix.vstack(alpha,alpha*Ap).det() != 0
assert s.Matrix.hstack(tau_m,Ap*tau_m).det() != 0
fp = lambda x: (1+x)/(2+x)
fm = lambda x: (2-x)/3
assert fm(F(1,2)) == F(1,2)
assert fm(fp(F(1,2))) == F(7,15)
lo = hi = F(1,2)
all_states = {F(1,2)}
frontier = {F(1,2)}
envelope_rows = []
for depth in range(1,10):
    candidates = [F(1,2),fp(lo),fp(hi),fm(lo),fm(hi)]
    lo,hi = min(candidates),max(candidates)
    frontier = {f(x) for x in frontier for f in (fp,fm)}
    all_states |= frontier
    assert lo == min(all_states) and hi == max(all_states)
    envelope_rows.append({"depth":depth,"lower":str(lo),"upper":str(hi),
                          "distinct_prefix_states":len(all_states)})

# Uniform clipped-ratio bound, including successors outside the anchors.
clip = lambda x: min(F(1), max(F(0), x))
ratio_cases = 0
outside_cases = 0
for left,gap in product((F(0),F(1,4),F(1,2)), (F(1,4),F(1,2))):
    right = left+gap
    if right > 1:
        continue
    eps = gap/8
    for value in [F(i,8) for i in range(9)]:
        original = clip((value-left)/gap)
        for dl,dr,dv in product((-eps,F(0),eps),repeat=3):
            le,re,ve = left+dl,right+dr,value+dv
            if not (0<=le<=1 and 0<=re<=1 and 0<=ve<=1):
                continue
            perturbed = clip((ve-le)/(re-le))
            assert abs(perturbed-original) <= 2*eps/(gap-2*eps)
            ratio_cases += 1
            outside_cases += int(value<left or value>right)

# Verify the local kernel TV bounds and downward-rounding properness.
def fraction(x):
    return F(int(s.numer(x)),int(s.denom(x)))

def tv(a,b):
    return sum(abs(x-y) for x,y in zip(a,b))/2

rounding_cases = 0
for depth in (1,2,4):
    aa=R(1,2)
    for _ in range(depth):
        aa=s.cancel((1+aa)/(2+aa))
    Bs=clipped_matrices(aa)
    for state in range(2):
        true_trans=[fraction(Bs[a][state,j]) for a in range(2) for j in range(2)]
        true_kernel=[F(1,2)]+true_trans
        masses=[F(1,2),sum(true_trans[:2]),sum(true_trans[2:])]
        weights=[true_trans[1]/masses[1],true_trans[3]/masses[2]]
        for sign,sign0,sign1,bits in product((-1,1),(-1,1),(-1,1),range(1,11)):
            eps_mass=F(1,50); eps_weight=F(1,50)
            estimated_masses=[masses[0]+sign*F(1,100),
                              masses[1]-sign*F(1,100),masses[2]]
            estimated_weights=[clip(weights[0]+sign0*eps_weight),
                               clip(weights[1]+sign1*eps_weight)]
            prelim=[]
            for letter in range(2):
                mass=estimated_masses[letter+1]
                w=estimated_weights[letter]
                prelim += [mass*(1-w),mass*w]
            grid=F(1,2**bits)
            rounded=[(v//grid)*grid for v in prelim]
            kernel=[1-sum(rounded)]+rounded
            assert sum(kernel)==1 and all(v>=0 for v in kernel)
            assert kernel[0]>=F(1,2)-eps_mass/2
            assert tv(true_kernel,kernel)<=eps_mass/2+eps_weight+4*grid
            rounding_cases+=1

results={
    "status":"all_checks_passed",
    "scope":"Exact algebra and finite instances; not proof-assistant verification or model training.",
    "sympy_version":s.__version__,
    "boundary_similarity":"exact",
    "properness_and_rank_witnesses":"exact",
    "cross_ratio_identity":"exact",
    "golden_clipping":golden_rows,
    "mixed_orientation_envelope":envelope_rows,
    "clipped_ratio_cases":ratio_cases,
    "outside_interval_ratio_cases":outside_cases,
    "kernel_rounding_cases":rounding_cases,
    "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
}
out=Path(__file__).with_name("rank_two_verification_results.json")
out.write_text(json.dumps(results,indent=2)+"\n")
print(json.dumps({k:v for k,v in results.items()
                  if k not in ("golden_clipping","mixed_orientation_envelope")},indent=2))
