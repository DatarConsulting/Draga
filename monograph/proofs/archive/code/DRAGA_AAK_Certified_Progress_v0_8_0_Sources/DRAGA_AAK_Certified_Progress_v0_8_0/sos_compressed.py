"""Order-two moment SDP on a necessary compressed Hankel rank condition.

Numerical discovery only. Exact checking is handled separately.
"""
from pathlib import Path
import itertools,json,time
import numpy as np
import scipy.sparse as sp
import sympy as sy
import cvxpy as cp
ROOT=Path(__file__).resolve().parent

from polynomial_model import build,monomials

def main():
    x,K,F,M,S,minors,ball,upper=build();n=len(x)
    mons=monomials(n,4);index={a:i for i,a in enumerate(mons)}
    b1=monomials(n,1);b2=monomials(n,2)
    def terms(p):return sy.Poly(p,*x).terms()
    def local(mat,basis):
        rows=[];cols=[];values=[];d=mat.rows;nb=len(basis);D=d*nb
        for i in range(d):
          for j in range(d):
            tt=terms(mat[i,j])
            for a,ma in enumerate(basis):
              for b,mb in enumerate(basis):
                row=(i*nb+a)*D+j*nb+b
                for mon,c in tt:
                  if c:
                    ex=tuple(aa+bb+cc for aa,bb,cc in zip(ma,mb,mon))
                    rows.append(row);cols.append(index[ex]);values.append(float(c))
        return sp.coo_matrix((values,(rows,cols)),shape=(D*D,len(mons))).tocsr(),D
    maps=[]
    for name,mat,basis in [('moment',sy.eye(1),b2),('norm',S,b1),
                           ('ball',sy.Matrix([ball]),b1),('u',sy.Matrix([x[-1]]),b1),
                           ('upper',sy.Matrix([upper]),b1)]:
        A,D=local(mat,basis);maps.append((name,A,D))
    rows=[];cols=[];vals=[];eq_polys=[]
    for p in minors:
      for mon in b1:
        multiplier=sy.prod(z**i for z,i in zip(x,mon));pol=sy.Poly(p*multiplier,*x)
        row=len(eq_polys);eq_polys.append(pol.as_expr())
        for m,c in pol.terms():
          if c:rows.append(row);cols.append(index[m]);vals.append(float(c))
    Aeq=sp.coo_matrix((vals,(rows,cols)),shape=(len(eq_polys),len(mons))).tocsr()
    y=cp.Variable(len(mons));constraints=[y[0]==1,Aeq@y==0]
    for name,A,D in maps:constraints.append(cp.reshape(A@y,(D,D),order='C')>>0)
    ui=index[tuple([0]*12+[1])]
    problem=cp.Problem(cp.Minimize(y[ui]),constraints)
    print('Built',len(mons),'moments; cones',[z[2] for z in maps],flush=True)
    value=problem.solve(solver='SCS',eps=2e-8,max_iters=100000,verbose=False)
    print('RESULT',problem.status,value,'sqrt',np.sqrt(max(value,0)),flush=True)
    data={'status':problem.status,'objective':float(value),'sqrt_objective':float(np.sqrt(max(value,0))),
          'normalization_dual':float(constraints[0].dual_value),
          'equality_dual':constraints[1].dual_value.tolist(),
          'cones':{name:constraints[i+2].dual_value.tolist() for i,(name,_,_) in enumerate(maps)},
          'cvxpy_version':cp.__version__,'solver':'SCS',
          'solve_seconds':problem.solver_stats.solve_time,'iterations':problem.solver_stats.num_iters}
    (ROOT/'sos_compressed_discovery.json').write_text(json.dumps(data,indent=2))
    np.savez_compressed(ROOT/'sos_compressed_moments.npz',y=y.value)

if __name__=='__main__':main()
