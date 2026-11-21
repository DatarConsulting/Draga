import sympy as s
import json
sq=s.sqrt;I=s.I
D=s.Matrix([[-s.Rational(4,7),sq(5)/7],[-5*sq(5)/14,-s.Rational(3,7)]])
l=s.Matrix([[-sq(7)/14,sq(35)/7]])
N=s.zeros(4,2);N[1:3,:]=D;N[3,:]=l
z=s.Matrix([0,1,0,0]);zb=z[:3,0]
B0=s.zeros(3,3);B0[:,1:3]=N[:3,:]
e0=s.zeros(1,3);e0[:,1:3]=N[3,:]
xs=s.symbols('x1 y1 y2 x3 y3 x4 y4',real=True)
a,b,c,d,e,f,g=xs
v=s.Matrix([a+I*b,I*c,d+I*e,f+I*g])
vnorm=a*a+b*b+c*c+d*d+e*e+f*f+g*g
hs=s.symbols('h0:16',real=True)
H=s.zeros(4,4);p=0
for i in range(4):H[i,i]=hs[p];p+=1
for i in range(4):
 for j in range(i+1,4):
  H[i,j]=hs[p]+I*hs[p+1];H[j,i]=hs[p]-I*hs[p+1];p+=2
K=I*H
N1=K*N;N2=(K*K*N-z*z.T*N)/2
B1=s.zeros(3,3);B1[:,0]=zb;B1[:,1:3]=N1[:3,:]
B2=s.zeros(3,3);B2[:,0]=v[:3,0];B2[:,1:3]=N2[:3,:]
e1=s.zeros(1,3);e1[0,0]=z[3];e1[:,1:3]=N1[3,:]
e2=s.zeros(1,3);e2[0,0]=v[3];e2[:,1:3]=N2[3,:]
assert s.simplify(N.T*N-s.eye(2)) == s.zeros(2)
J0=s.simplify((e0*B0*B0*zb)[0]);assert J0 == sq(7)/4
print('J0',J0,flush=True)
J1=s.expand((e1*B0*B0*zb+e0*(B1*B0+B0*B1)*zb+e0*B0*B0*v[:3,0])[0]);print('J1 real',s.simplify(s.re(J1)),flush=True)
J2=s.expand((e2*B0*B0*zb+e1*(B1*B0+B0*B1)*zb+e1*B0*B0*v[:3,0]+e0*(B2*B0+B0*B2+B1*B1)*zb+e0*(B1*B0+B0*B1)*v[:3,0])[0]-vnorm*J0/2)
F2=s.expand(s.re(J1*s.conjugate(J1))+2*J0*s.re(J2))
vs=xs+hs
Hess=s.hessian(F2,vs);grad=s.Matrix([s.diff(F2,x).subs(dict.fromkeys(vs,0)) for x in vs]);base=s.simplify(F2.subs(dict.fromkeys(vs,0)))
print('base',base,'nonzero gradient',[(str(x),s.simplify(y)) for x,y in zip(vs,grad) if y!=0],flush=True)
print('Hrank...',flush=True);rank=Hess.rank();print('rank',rank,flush=True);assert rank == 13
solution=s.linsolve((Hess,-grad));sol=next(iter(solution));free=set().union(*[x.free_symbols for x in sol]);part=s.Matrix([s.simplify(x.subs({f:0 for f in free})) for x in sol]);value=s.simplify(base+(grad.T*part)[0]/2)
assert value == s.Rational(1,7)
assert s.simplify(Hess*part+grad) == s.zeros(23,1)
print('max coeff',value,'solution',part.T,flush=True)
# record exact matrix and verify semidefiniteness by characteristic polynomial
cp=s.factor(Hess.charpoly().as_expr());print('charpoly',cp,flush=True)
y=s.Symbol('lambda');expected=y**10*(8*y+7)**4*(16*y*y+42*y+25)**2*(10976*y*y+56620*y+53525)*(32*y**3+308*y*y+655*y+350)/s.Integer(368293445632)
assert s.expand(cp-expected) == 0
print('Exact checks passed: isometry, boundary value, Hessian rank, critical response, characteristic polynomial.',flush=True)
open('research_v4_numeric/k3_hessian_exact.txt','w').write('variables='+str(vs)+'\nbase='+str(base)+'\ngrad='+str(grad.T)+'\nHessian='+str(Hess)+'\ncharpoly='+str(cp)+'\nparticular='+str(part.T)+'\nmaxcoef='+str(value)+'\n')
