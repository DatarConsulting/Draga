"""Exact response of the degree-three Brevig--Saksman extremizer.

No global classification or publication-priority claim is encoded here.
Arithmetic takes place in Q(alpha), 64 alpha^4-30 alpha^2+3=0.
"""
import json
from pathlib import Path
import sympy as sp
from sympy.polys.matrices import DomainMatrix as DM

alpha = sp.sqrt(15-sp.sqrt(33))/8
K = sp.QQ.algebraic_field(alpha)
aa = K.from_sympy(alpha)
bb = 6*aa-K.convert(sp.Rational(64,3))*aa**3
cc = 2-8*aa**2
zero, one = K.zero, K.one

def mat(rows):
    rr=[[K.convert(v) for v in row] for row in rows]
    return DM(rr,(len(rr),len(rr[0])),K)

def scalar(X):
    assert X.shape == (1,1)
    return X[0,0].element

def scale(X, x):
    return X.scalarmul(x)

m=3
T=mat([[0,1,0],[0,0,1],[-aa,-cc,-bb]])
l=mat([[1,0,0]])
I=DM.eye((m,m),K).to_dense()
tp=[I]
for _ in range(m):
    tp.append(tp[-1]*T)

# Solve Q-T^T Q T=l^T l exactly, row-major vectorization.
TR=T.to_list()
rows=[]
rhs=[]
for i in range(m):
    for j in range(m):
        rows.append([
            (one if (a==i and b==j) else zero)-TR[a][i]*TR[b][j]
            for a in range(m) for b in range(m)
        ])
        rhs.append([one if (i==0 and j==0) else zero])
sol=mat(rows).inv()*mat(rhs)
flat=[sol[i,0].element for i in range(m*m)]
Q=mat([flat[m*i:m*(i+1)] for i in range(m)])
assert (Q-T.transpose()*Q*T-l.transpose()*l).is_zero_matrix

row=l*tp[m]
x0=Q.inv()*row.transpose()
G=scalar(row*x0)
cs=[scalar(l*tp[j]*x0) for j in range(m+1)]

S0=zero
v0=mat([[0],[0],[0]])
M0=DM.zeros((m,m),K).to_dense()
for j in range(m):
    S0+=cs[j]*scalar(x0.transpose()*Q*tp[m-j]*x0)
for j in range(m-1):
    v0+=scale(tp[m-2-j]*x0,cs[j])
for j in range(m+1):
    M0+=tp[m-j]*x0*l*tp[j]
K0=scale(M0+Q.inv()*M0.transpose()*Q,K.convert(sp.Rational(1,2)))
assert (Q*K0-K0.transpose()*Q).is_zero_matrix

b=(cs[m-1]**2-S0+scalar(v0.transpose()*Q*K0.inv()*v0))/G
H0=Q*K0
cert_Q=mat([
    [(342+768*aa**2)/333,(1116*aa-4224*aa**3)/333,(528-2304*aa**2)/333],
    [(1116*aa-4224*aa**3)/333,(798-576*aa**2)/333,(1980*aa-6272*aa**3)/333],
    [(528-2304*aa**2)/333,(1980*aa-6272*aa**3)/333,(858-192*aa**2)/333],
])
assert (Q-cert_Q).is_zero_matrix
cert_H=mat([
    [K.convert(sp.Rational(21,32))+7*aa**2/16,aa/2-3*aa**3,K.convert(sp.Rational(5,8))-9*aa**2/4],
    [aa/2-3*aa**3,K.convert(sp.Rational(7,4))-5*aa**2/2,3*aa-20*aa**3/3],
    [K.convert(sp.Rational(5,8))-9*aa**2/4,3*aa-20*aa**3/3,3-19*aa**2/3],
])
assert (H0-cert_H).is_zero_matrix
assert S0==K.convert(sp.Rational(6669,65536))-K.convert(sp.Rational(1809,32768))*aa**2
assert cs[2]**2==K.convert(sp.Rational(6075,32768))-K.convert(sp.Rational(7767,16384))*aa**2
inverse_term=scalar(v0.transpose()*Q*K0.inv()*v0)
assert inverse_term==K.convert(sp.Rational(171,65536))+K.convert(sp.Rational(5625,32768))*aa**2
assert b==K.convert(sp.Rational(75,496))-K.convert(sp.Rational(103,248))*aa**2
Q_minors=[Q.extract(range(k),range(k)).det() for k in range(1,4)]
H_minors=[H0.extract(range(k),range(k)).det() for k in range(1,4)]
expected_Q=[
    K.convert(sp.Rational(38,37))+K.convert(sp.Rational(256,111))*aa**2,
    K.convert(sp.Rational(244,111))+K.convert(sp.Rational(896,333))*aa**2,
    K.convert(sp.Rational(968,333))+K.convert(sp.Rational(11008,999))*aa**2,
]
expected_H=[
    K.convert(sp.Rational(21,32))+K.convert(sp.Rational(7,16))*aa**2,
    K.convert(sp.Rational(1287,1024))-K.convert(sp.Rational(915,512))*aa**2,
    K.convert(sp.Rational(153,64))-K.convert(sp.Rational(189,32))*aa**2,
]
assert Q_minors==expected_Q
assert H_minors==expected_H
Gamma2=sp.simplify(K.to_sympy(G))
response=sp.simplify(K.to_sympy(b))
expected_G=1-sp.Rational(1153,1)/(2*(1103+33*sp.sqrt(33)))
assert sp.simplify(Gamma2-expected_G)==0

result={
    'Gamma2':str(Gamma2),
    'response_b':str(response),
    'response_b_decimal':str(sp.N(response,70)),
    'boundary_error_squared':str(sp.simplify(1-Gamma2)),
    'checks':[
        'exact Stein identity',
        'exact Q-selfadjoint response matrix',
        'exact agreement with Brevig--Saksman leading constant',
        'all displayed Q and H0 matrix entries',
        'all displayed scalar contractions',
        'all displayed principal-minor formulas',
        'exact quartic-field response reduction',
    ],
    'scope':'Response of the known extremizer; globality uses the extremizer-classification bridge.'
}
print(json.dumps(result,indent=2))
Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
