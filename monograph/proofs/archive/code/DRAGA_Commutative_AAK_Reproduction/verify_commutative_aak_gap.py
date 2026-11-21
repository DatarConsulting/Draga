#!/usr/bin/env python3
"""Exact certificates for DRAGA_Commutative_AAK_Gap.pdf.

Default: Python standard library only. Optional --numerical uses NumPy/SciPy.
Finite tests verify algebraic identities and certificates; the universal
lower-bound argument is written in the accompanying mathematical note.
"""
from fractions import Fraction as F
from itertools import product
import argparse
import json


def matmul(a, b):
    return [[sum((x*y for x, y in zip(row, col)), F(0))
             for col in zip(*b)] for row in a]


def determinant(a):
    a = [row[:] for row in a]
    ans = F(1)
    for j in range(len(a)):
        p = next((i for i in range(j, len(a)) if a[i][j]), None)
        if p is None:
            return F(0)
        if p != j:
            a[j], a[p] = a[p], a[j]
            ans = -ans
        pivot = a[j][j]
        ans *= pivot
        for i in range(j+1, len(a)):
            fac = a[i][j] / pivot
            for k in range(j+1, len(a)):
                a[i][k] -= fac*a[j][k]
    return ans


def padd(a, b):
    return [(a[i] if i < len(a) else F(0)) +
            (b[i] if i < len(b) else F(0))
            for i in range(max(len(a), len(b)))]


def pscale(a, s):
    return [s*x for x in a]


def pmul(a, b):
    out = [F(0)]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return out


def cube_reduce(p):
    p = p[:] + [F(0)]*max(0, 3-len(p))
    for i in range(len(p)-1, 2, -1):
        p[i-3] += 2*p[i]
        p[i] = F(0)
    return p[:3]


def exact_checks():
    J = [[F(x) for x in row] for row in [[0,1,0],[0,0,1],[0,0,0]]]
    J2 = matmul(J,J)
    assert matmul(J,J2) == [[0]*3 for _ in range(3)]
    assert matmul(J,J2) == matmul(J2,J)
    alpha = [[F(1),F(0),F(0)]]
    beta = [[F(1)],[F(0)],[F(1)]]
    reach = alpha + matmul(alpha,J) + matmul(alpha,J2)
    obs = [
        [beta[i][0],matmul(J,beta)[i][0],matmul(J2,beta)[i][0]] for i in range(3)]
    assert determinant(reach) == 1 and determinant(obs) == -1
    coefficients = {"":F(1),"b":F(1),"aa":F(1)}
    # J^3=0 proves the all-length support statement. This loop checks indexing.
    for n in range(6):
        for letters in product("ab",repeat=n):
            word="".join(letters); row=alpha
            for a in word:
                row=matmul(row,J if a=="a" else J2)
            assert matmul(row,beta)[0][0] == coefficients.get(word,F(0))
    words=["","a","b","aa"]
    H=[[coefficients.get(u+v,F(0)) for v in words] for u in words]
    for v, lam in [([2,0,1,1],2),([0,1,0,0],1),([1,0,-1,-1],-1),([0,0,1,-1],0)]:
        assert matmul(H,[[F(x)] for x in v]) == [[F(lam*x)] for x in v]
    # Exact coefficient identity in y=t-1/2, not evaluation at a grid.
    target=[F(1,16),F(-3,2),F(33,2),F(6),F(-3)]
    square=pscale(pmul([F(-1,17),F(1)],[F(-1,17),F(1)]),F(51,4))
    positive=pscale(pmul([F(0),F(0),F(1)],pmul([F(1,2),F(1)],[F(5,2),F(-1)])),F(3))
    assert padd(padd(square,[F(5,272)]),positive)==target
    assert F(4,3)+(F(5,1224)/F(4,3)) == F(727,544)
    # Upper certificate on the complete invariant five-dimensional space.
    b,c,s=F(13,24),F(197,200),F(29,25)
    d=b**4/(1-b*b); D=[F(1)]*4+[d]; k=[F(1),F(0),b,F(0),F(1)]
    H5=[row+[F(0)] for row in H]+[[F(0)]*5]
    L=[[D[i]*H5[i][j]-c*D[i]*k[i]*k[j]*D[j] for j in range(5)] for i in range(5)]
    certificates={}
    for sign in [1,-1]:
        C=[[(s*D[i] if i==j else F(0))+sign*L[i][j] for j in range(5)] for i in range(5)]
        assert C==[list(r) for r in zip(*C)]
        minors=[determinant([row[:n] for row in C[:n]]) for n in range(1,6)]
        assert all(x>0 for x in minors)
        certificates[str(sign)]=[str(x) for x in minors]
    # Isospectral unary model in Q[r]/(r^3-2).
    c0,c1,c2=[F(2),F(-1)],[F(-1),F(0),F(1)],[F(0),F(1)]
    trace=padd(c0,c2)
    trace_squared=padd(padd(pmul(c0,c0),pscale(pmul(c1,c1),F(2))),pscale(pmul(c2,c2),F(3)))
    det=pscale(pmul(pmul(c2,c2),c2),F(-1))
    assert cube_reduce(trace)==[2,0,0]
    assert cube_reduce(trace_squared)==[6,0,0]
    assert cube_reduce(det)==[-2,0,0]
    # Explicit optimal unary approximant: b0=1/r=r^2/2 and c=2-r.
    b0=[F(0),F(0),F(1,2)]
    b02=cube_reduce(pmul(b0,b0))
    b06=cube_reduce(pmul(pmul(b02,b02),b02))
    assert b06==[F(1,4),0,0]
    Kv=[padd(padd(c0,pmul(c1,b0)),pmul(c2,b02)),
        padd(c1,pmul(c2,b0)),c2]
    for left,right in zip(Kv,[[F(2)],pscale(b0,F(2)),pscale(b02,F(2))]):
        assert cube_reduce(left)==cube_reduce(right)
    assert cube_reduce(c0)==cube_reduce(pscale(padd([F(1)],pscale(b02,F(-1))),F(2)))
    # Proper probability realization for f/3.
    Aa=[[F(0),F(1,3),F(0)],[F(0),F(0),F(1)],[F(0)]*3]
    Ab=[[F(0),F(0),F(1,3)],[F(0)]*3,[F(0)]*3]
    term=[F(1,3),F(0),F(1)]
    assert all(term[i]+sum(Aa[i])+sum(Ab[i])==1 for i in range(3))
    return {"status":"EXACT_CERTIFICATES_PASS", "lower_squared":"727/544",
            "upper":"29/25", "upper_is_strict":True,
            "upper_candidate":{"a":"0","b":str(b),"c":str(c)},
            "upper_leading_principal_minors":certificates,
            "H_singular_values":[2,1,1], "unary_singular_values":[2,1,1],
            "scope":"Algebraic identities and finite certificates; universal proof in accompanying note."}


def numerical_checks():
    import numpy as np
    from scipy.optimize import differential_evolution
    H=np.array([[1,0,1,1],[0,1,0,0],[1,0,0,0],[1,0,0,0]],float)
    out={"status":"HEURISTIC_NOT_GLOBAL_CERTIFICATE"}
    for is_complex in [False,True]:
        def objective(v):
            p,q=(v[0]+1j*v[1],v[2]+1j*v[3]) if is_complex else v
            m=1-abs(p)**2-abs(q)**2
            if m<=0:return 20-m
            w=H@np.array([1,p,q,p*p])
            M=H@H-m*np.outer(w,w.conj())
            return float(np.sqrt(max(0,np.linalg.eigvalsh(M)[-1])))
        opt=differential_evolution(objective,[(-1,1)]*(4 if is_complex else 2),seed=19,tol=1e-10)
        out["complex" if is_complex else "real"]={"value":opt.fun,"parameters":opt.x.tolist()}
    # Check the characteristic reduction at fixed complex test points.
    for p,q in [(0.2+0.1j,0.3-0.2j),(0.1j,0.8),(0,0.999)]:
        x=abs(p)**2;t=abs(q);m=1-x-t*t;z=q+p*p
        w=H@np.array([1,p,q,p*p]);M=H@H-m*np.outer(w,w.conj())
        ev=np.linalg.eigvalsh(M);lam=4/3
        Q=(2/3)*(m*(4*z.real-abs(z)**2-4*x)-4/3)
        assert abs(Q-(lam-ev[-1])*(lam-ev[1]))<1e-12
        assert ev[-1]>=727/544-1e-12
    return out


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--numerical",action="store_true")
    args=parser.parse_args()
    result={"exact":exact_checks()}
    if args.numerical:result["numerical"]=numerical_checks()
    print(json.dumps(result,indent=2))
