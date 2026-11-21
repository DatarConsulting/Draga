"""Exact polynomial data for a necessary four-dimensional compression."""
import itertools
import sympy as sy

def monomials(n,k):
    out=[]
    for d in range(k+1):
        for ind in itertools.combinations_with_replacement(range(n),d):
            a=[0]*n
            for i in ind:a[i]+=1
            out.append(tuple(a))
    return out

def build():
    x=sy.symbols('x0:13');u=x[12]
    K=sy.Matrix([[x[0],x[1],x[2],x[3]+x[6]],
                 [x[1],x[3],x[4],x[7]],
                 [x[2],x[5],x[6],x[8]],
                 [x[3]+x[6],x[9],x[10],x[11]]])
    F=sy.Matrix([[0,0,0,2],[0,1,0,sy.Rational(1,100)],
                 [0,0,1,sy.Rational(1,100)],[2,sy.Rational(1,100),sy.Rational(1,100),0]])
    M=sy.diag(1,1,1,2);E=F-K;S=u*M-E.T*M.inv()*E
    subsets=list(itertools.combinations(range(4),3))
    minors=[sy.expand(K.extract(a,b).det()) for a in subsets for b in subsets]
    ball=sy.Integer(80)-sum(z*z for z in x[:12])
    upper=sy.Rational(6533,5000)**2-u
    return x,K,F,M,S,minors,ball,upper
