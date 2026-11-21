"""Rebuild the three-state complex Hessian from the manuscript's local chart.
Exact symbolic computation; this checks an algebraic step, not the global proof.
Run: python3 verify_complex_hessian.py
"""
from pathlib import Path
import sympy as s
import json,time
start=time.time()
I=s.I
x1,y1,y2,x3,y3,x4,y4=s.symbols('x1 y1 y2 x3 y3 x4 y4',real=True)
h=s.symbols('h0:16',real=True)
variables=(x1,y1,y2,x3,y3,x4,y4)+h
zero={v:0 for v in variables}
z0=s.Matrix([0,1,0,0]); v=s.Matrix([x1+I*y1,I*y2,x3+I*y3,x4+I*y4])
norm2=sum(vv*s.conjugate(vv) for vv in v).expand()
N=s.Matrix([[0,0],[-s.Rational(4,7),s.sqrt(5)/7],[-5*s.sqrt(5)/14,-s.Rational(3,7)],[-s.sqrt(7)/14,s.sqrt(35)/7]])
HH=s.diag(*h[:4]); k=4
for i in range(4):
 for j in range(i+1,4):
  HH[i,j]=h[k]+I*h[k+1]; HH[j,i]=h[k]-I*h[k+1]; k+=2
V=[s.zeros(4,1).row_join(N), z0.row_join(I*HH*N), v.row_join(-(HH*HH*N+z0*z0.T*N)/2)]
B=[a[:3,:] for a in V]; ell=[a[3:4,:] for a in V]
z=[z0[:3,:],v[:3,:],-norm2*z0[:3,:]/2]
J=[]
for degree in range(3):
 value=0
 for a in range(degree+1):
  for b in range(degree-a+1):
   for c in range(degree-a-b+1):
    d=degree-a-b-c
    value+=(ell[a]*B[b]*B[c]*z[d])[0]
 J.append(s.expand(value))
print('Built second-order chart',flush=True)
quad=s.expand(J[1]*s.conjugate(J[1])+s.conjugate(J[0])*J[2]+J[0]*s.conjugate(J[2]))
assert s.simplify(J[0]*s.conjugate(J[0]))==s.Rational(7,16)
assert s.expand(s.conjugate(J[0])*J[1]+J[0]*s.conjugate(J[1]))==0
a=quad.subs(zero)
g=s.Matrix([s.diff(quad,vv).subs(zero) for vv in variables])
H=s.Matrix([[s.diff(quad,u,vv) for vv in variables] for u in variables])
assert a==-s.Rational(23,112)
expected_g=-3*x1/4-h[5]/7-5*s.sqrt(5)*h[7]/56-s.sqrt(7)*h[9]/56
assert s.expand((g.T*s.Matrix(variables))[0]-expected_g)==0
print('Checked constant and gradient',flush=True)
# Split the symmetric matrix into connected coordinate blocks before charpoly.
unused=set(range(23)); blocks=[]
while unused:
 block={unused.pop()}; todo=list(block)
 while todo:
  i=todo.pop()
  found={j for j in unused if H[i,j]!=0}
  unused-=found; block|=found; todo.extend(found)
 blocks.append(sorted(block))
lam=s.symbols('lambda'); actual=1
for block in blocks:
 actual*=H.extract(block,block).charpoly(lam).as_expr()
actual=s.factor(actual,extension=s.sqrt(5))
expected=lam**10*(8*lam+7)**4*(16*lam**2+42*lam+25)**2*(10976*lam**2+56620*lam+53525)*(32*lam**3+308*lam**2+655*lam+350)/s.Integer(368293445632)
assert s.expand(actual-expected)==0
xi=s.zeros(23,1); xi[0]=-s.Rational(6,7); xi[7+5]=-s.Rational(2,7); xi[7+7]=-s.sqrt(5)/35
assert (H*xi+g).applyfunc(s.simplify)==s.zeros(23,1)
assert s.simplify(a+(g.T*xi)[0]/2)==s.Rational(1,7)
assert H==H.T
report={'scope':'Exact check of the local complex Hessian and optimized coefficient; not a global optimality proof.',
 'passed':True,'coordinates':23,'block_sizes':[len(b) for b in blocks],
 'hessian_rank':13,'constant':str(a),'optimized_coefficient':'1/7',
 'characteristic_polynomial':str(s.factor(expected)),
 'sympy':s.__version__,'elapsed_seconds':time.time()-start}
Path(__file__).with_name('complex_hessian_results.json').write_text(json.dumps(report,indent=2))
print(json.dumps(report,indent=2))
