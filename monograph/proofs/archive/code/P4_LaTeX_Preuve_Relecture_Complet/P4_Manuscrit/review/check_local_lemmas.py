"""Narrow, independently written rational checks for the P4 referee review.

These checks do not prove the full theorem or replay historical campaigns.
They use only Python's standard library and do not import the codec.
"""
from fractions import Fraction as F
from random import Random
from pathlib import Path
import json


def verify_splits():
    checks = 0
    for den in range(2, 33):
        for k in range(1, den):
            q = F(k, den)
            b = min(q, 1-q)
            loss = (F(3, 2)-b)/b
            for capacity in range(1, 100):
                x = (capacity-1)*q + F(1, 2)
                left = x.numerator // x.denominator
                right = capacity-1-left
                assert left >= q*capacity-q-F(1, 2)
                assert right >= (1-q)*capacity-(1-q)-F(1, 2)
                assert q*loss+q+F(1, 2) <= loss
                assert (1-q)*loss+(1-q)+F(1, 2) <= loss
                if capacity > 1+1/(2*b):
                    assert left >= 1 and right >= 1
                checks += 1
    return checks


def multiply(X, Y):
    a,b,c,d = X
    A,B,C,D = Y
    return a*A+b*C, a*B+b*D, c*A+d*C, c*B+d*D


def determinant(X):
    a,b,c,d = X
    return a*d-b*c


def point(X, x):
    a,b,c,d = X
    return (a*x+b*(1-x))/((a+c)*x+(b+d)*(1-x))


def mass(X):
    return F(4*determinant(X), sum(X)**2)


def verify_projective_identities():
    rng = Random(20261103)
    recursions = 0
    for _ in range(400):
        X = (1,0,0,1)
        for _ in range(rng.randrange(1,12)):
            while True:
                Y = tuple(rng.randrange(1,9) for _ in range(4))
                if determinant(Y) > 0:
                    break
            t = F(X[0]+X[2], sum(X))
            a,b,c,d = Y
            denominator = (a+b)*t+(c+d)*(1-t)
            G = (a*t+c*(1-t))/denominator
            g = F(determinant(Y))/denominator**2
            product = multiply(X,Y)
            assert F(product[0]+product[2],sum(product)) == G
            assert mass(product) == mass(X)*g
            X = product
            recursions += 1
        left,right = F(1,4),F(3,4)
        t = F(X[0]+X[2],sum(X))
        z = lambda x: t*x+(1-t)*(1-x)
        assert point(X,right)-point(X,left) == (right-left)*mass(X)/(4*z(left)*z(right))
    return recursions,400


if __name__ == '__main__':
    split_count = verify_splits()
    recursions,widths = verify_projective_identities()
    result = {
        'passed': True,
        'seed': 20261103,
        'split_checks': split_count,
        'dual_recursion_checks': recursions,
        'width_identity_checks': widths,
        'production_code_imported': False,
        'full_P4_theorem_verified': False,
        'historical_campaigns_replayed': False,
    }
    target = Path(__file__).with_name('local_lemma_checks.json')
    target.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result))
