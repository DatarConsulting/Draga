#!/usr/bin/env python3
"""Exact algebra checks and high-precision diagnostics for the v4 generalization.

Numeric results are diagnostics, not interval certificates or a proof of global
optimality. Infinite Hankel norms use the 7-state target-minus-competitor error
realization and its full Stein Gramians; no finite word horizon is substituted.

Dependencies: mpmath, sympy.
"""
from __future__ import annotations
import argparse
import json
import platform
from pathlib import Path
import mpmath as mp
import sympy as sp

K3_CASES = ('0.01','0.001','0.0001','0.00001','0.000001')
F1_CASES = ('0.0001','0.01','0.1','0.5','0.9','0.99')

def adj(x): return x.transpose_conj()
def maxabs(x): return max(abs(x[i,j]) for i in range(x.rows) for j in range(x.cols))
def herm(x): return (x+adj(x))/2
def num(x): return mp.nstr(mp.re(x),mp.mp.dps)

def exact_algebra_checks():
    """Independent exact Q-metric response and chart passivity identities."""
    T=sp.Matrix([[0,1],[-sp.Rational(1,2),-1]])
    l=sp.Matrix([[1,0]])
    Q=sp.Matrix([[sp.Rational(8,5),sp.Rational(4,5)],[sp.Rational(4,5),sp.Rational(12,5)]])
    gamma=sp.sqrt(7)/4
    x=sp.Matrix([-1,-3])/(2*sp.sqrt(7))
    a=[sp.simplify((l*T**j*x)[0]) for j in range(3)]
    S=sp.simplify(sum(a[j]*(x.T*Q*T**(2-j)*x)[0] for j in range(2)))
    v=sp.simplify(a[0]*x)
    R=sum((T**(2-j)*x*l*T**j for j in range(3)),sp.zeros(2))
    K=sp.simplify((R+Q.inv()*R.T*Q)/2)
    response=sp.simplify((v.T*Q*K.inv()*v)[0])
    b=sp.simplify(a[1]**2-gamma*S+gamma*response)
    assert Q-T.T*Q*T == l.T*l
    assert sp.simplify((x.T*Q*x)[0]) == 1
    assert a[2] == gamma
    assert sp.simplify(K.T*Q-Q*K) == sp.zeros(2)
    assert S == 23*sp.sqrt(7)/196
    assert response == 3*sp.sqrt(7)/196
    assert b == sp.Rational(1,7)
    # Abstract exact chart identity, with p=N^*z and s^2+r^2=1.
    s,r,p1,p2=sp.symbols('s r p1 p2',real=True)
    p=sp.Matrix([p1,p2])
    gram=sp.zeros(3)
    gram[0,0]=s*s;gram[0,1:3]=s*r*p.T;gram[1:3,0]=s*r*p
    gram[1:3,1:3]=sp.eye(2)-s*s*p*p.T
    av=sp.Matrix([r,-s*p1,-s*p2])
    passive=sp.expand(gram+av*av.T-sp.eye(3)).subs(r*r,1-s*s)
    assert passive == sp.zeros(3)
    q,t=sp.symbols('q t')
    A=(1-q)**2*t**2+(q*q+2*q-1)*(t-q)
    L=(1+q)*(t-q)-(1-q)*t*t
    identity=sp.expand(A*A-4*q**3*(t-q)*(1-2*q+q*t)
        -(1-q)**2*(L*L-4*q*(1-t)*(t-q)**2))
    assert identity == 0
    return {'status':'Exact symbolic assertions passed.', 'b3':'1/7',
            'response_in_Q_metric':str(response),
            'checks':['lossless Q Stein identity','Q-normalized extremal vector',
                      'Q-selfadjoint response matrix','exact normal-response coefficient',
                      'abstract passive chart identity','F1 quartic/radical polynomial identity']}

def stein(letters,source,side):
    n=source.rows
    system=mp.eye(n*n);rhs=mp.matrix(n*n,1)
    for i in range(n):
        for j in range(n):
            ij=i*n+j;rhs[ij]=source[i,j]
            for k in range(n):
                for l in range(n):
                    if side=='row':
                        c=sum(a[i,k]*mp.conj(a[j,l]) for a in letters)
                    elif side=='column':
                        c=sum(mp.conj(a[k,i])*a[l,j] for a in letters)
                    else: raise ValueError(side)
                    system[ij,k*n+l]-=c
    sol=mp.lu_solve(system,rhs);raw=mp.matrix(n,n)
    for i in range(n):
        for j in range(n):raw[i,j]=sol[i*n+j]
    symmetry=maxabs(raw-adj(raw));out=herm(raw)
    terms=(a*out*adj(a) for a in letters) if side=='row' else (adj(a)*out*a for a in letters)
    res=out-source-sum(terms,mp.zeros(n))
    relative=maxabs(res)/max(mp.mpf(1),maxabs(out),maxabs(source))
    return out,relative,symmetry

def block_diag(a,b):
    out=mp.zeros(a.rows+b.rows,a.cols+b.cols)
    for i in range(a.rows):
        for j in range(a.cols):out[i,j]=a[i,j]
    for i in range(b.rows):
        for j in range(b.cols):out[a.rows+i,a.cols+j]=b[i,j]
    return out

def k3_case(q_text,tol):
    q=mp.mpf(q_text);s=mp.sqrt(q);r=mp.sqrt(1-q)
    D=mp.matrix([[-mp.mpf(4)/7,mp.sqrt(5)/7],[-5*mp.sqrt(5)/14,-mp.mpf(3)/7]])
    ell=mp.matrix([[-mp.sqrt(7)/14,mp.sqrt(35)/7]])
    lam=-mp.mpf(6)/7;w=mp.matrix([[-mp.mpf(3)/14,mp.sqrt(5)/35]])
    ww=(w*adj(w))[0]
    F=mp.eye(2)-(1-mp.sqrt(1-q*ww))/ww*(adj(w)*w)
    N=mp.matrix(4,2)
    N[0,:]=s*w;N[1:3,:]=D*F;N[3,:]=ell*F
    z=mp.matrix([s*lam,mp.sqrt(1-q*lam*lam),0,0])
    W=mp.eye(4)-(1-r)*z*adj(z)
    V=mp.matrix(4,3);V[:,0]=s*z;V[:,1:3]=W*N
    B=V[:3,:];ls=V[3,:]
    av=mp.matrix(3,1);av[0]=r;av[1:3,:]=-s*adj(N)*z
    A=mp.matrix([[av[j] for j in range(3)],[0,0,0],[0,0,0]])
    resolvent=(mp.eye(3)-r*A)**-1
    omega=s*ls*(B**3)*resolvent
    boundary=(mp.eye(3)-r*A)*adj(omega)/s
    overlap=mp.re((omega*adj(omega))[0]);projection=1-overlap
    passive_res=maxabs(adj(A)*A+adj(B)*B+adj(ls)*ls-mp.eye(3))
    n_isometry=maxabs(adj(N)*N-mp.eye(2))
    z_unit=abs((adj(z)*z)[0]-1)
    interpolation=maxabs(s*resolvent*boundary-adj(omega))
    observation,ores,osym=stein((A,B),adj(ls)*ls,'column')
    observation_identity=maxabs(observation-mp.eye(3))
    Af=mp.diag([r,0,0,0]);Bf=mp.zeros(4)
    for i in range(3):Bf[i,i+1]=1
    letters=(block_diag(Af,adj(A)),block_diag(Bf,adj(B)))
    beta=mp.matrix(7,1);beta[3]=1;beta[4:7,:]=adj(ls)
    alpha=mp.matrix(1,7);alpha[0]=q;alpha[:,4:7]=-adj(boundary)
    P,pres,psym=stein(letters,beta*adj(beta),'row')
    Q,qres,qsym=stein(letters,adj(alpha)*alpha,'column')
    L=mp.cholesky(P)
    roots,_=mp.eighe(herm(adj(L)*Q*L))
    error2=mp.re(roots[roots.rows-1])
    gram_minP=mp.eighe(P,eigvals_only=True)[0]
    gram_minQ=mp.eighe(Q,eigvals_only=True)[0]
    # Numerical stability check on the 3-state candidate's complete CP map.
    cp=mp.matrix(9)
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    cp[3*i+j,3*k+l]=sum(a[i,k]*mp.conj(a[j,l]) for a in (A,B))
    eigs,_=mp.eig(cp);rho=max(abs(v) for v in eigs)
    remainder=(error2-(mp.mpf(9)/16-q/7))/(q*q)
    projection_gap=(error2-projection)/(q*q)
    overlap_remainder=(overlap-(mp.mpf(7)/16+q/7))/(q*q)
    residuals=max(passive_res,n_isometry,z_unit,interpolation,observation_identity,pres,qres,ores,psym,qsym,osym)
    passed=(residuals<tol and gram_minP>0 and gram_minQ>0 and rho<1
            and error2>=projection-tol and abs(remainder)<100 and abs(projection_gap)<100)
    return {'q':q_text,'squared_infinite_hankel_error':num(error2),
            'selected_root_of_error_Gramian':num(error2),
            'second_largest_squared_singular_value':num(roots[roots.rows-2]),
            'projection_residual_squared':num(projection),'overlap_squared':num(overlap),
            'error_first_correction_quotient':num((error2-mp.mpf(9)/16)/q),
            'error_remainder_div_q_squared':num(remainder),
            'error_minus_projection_div_q_squared':num(projection_gap),
            'overlap_remainder_div_q_squared':num(overlap_remainder),
            'passive_identity_residual':num(passive_res),
            'N_isometry_residual':num(n_isometry),'z_unit_residual':num(z_unit),
            'exact_interpolation_residual':num(interpolation),
            'observation_Gramian_minus_identity':num(observation_identity),
            'maximum_relative_Stein_residual':num(max(pres,qres,ores)),
            'maximum_Gramian_Hermitian_residual':num(max(psym,qsym,osym)),
            'reachability_Gramian_minimum_eigenvalue':num(gram_minP),
            'observability_Gramian_minimum_eigenvalue':num(gram_minQ),
            'candidate_joint_spectral_radius':num(rho),
            'joint_stability_gap_div_q':num((1-rho)/q),'passed':bool(passed)}

def bisect(lo,hi,f,increasing):
    for _ in range(5*mp.mp.dps+65):
        mid=(lo+hi)/2
        if mid==lo or mid==hi:break
        if (f(mid)<0)==increasing:lo=mid
        else:hi=mid
    return (lo+hi)/2

def f1_case(q_text,tol):
    q=mp.mpf(q_text)
    def udiff(u):
        return q*(u-q)**2/((1-q)**2*u*u)-(1-q)*(q-u*u)/(q*(1-u*u))
    u=bisect(q,mp.sqrt(q),udiff,True)
    quartic_value=(q*q+(1-2*q)*u*u)/(q*(1-u*u))
    def radical(t):return 1-t-q/(1+mp.sqrt(q*t/(t-q)))**2
    radical_value=bisect(q,mp.mpf(1),radical,False)
    t=radical_value
    selector=(1-q)**2*t*t+(q*q+2*q-1)*(t-q)
    rhs=4*q**3*(t-q)*(1-2*q+q*t)
    quartic_relative=abs(selector**2-rhs)/max(abs(selector**2),abs(rhs))
    diff=abs(t-quartic_value)
    passed=(diff<tol and abs(radical(t))<tol and quartic_relative<tol
            and 2*q/(1+q)<t<1 and selector>0 and q<u<mp.sqrt(q))
    return {'q':q_text,'u_selected_root':num(u),
            'value_from_u_selected_quartic':num(quartic_value),
            'value_from_monotone_radical_equation':num(radical_value),
            'value_difference':num(diff),'radical_equation_residual':num(abs(radical(t))),
            'direct_quartic_relative_residual':num(quartic_relative),
            'direct_quartic_positive_selector':num(selector),
            'dominant_vector_lower_bound':num(1-q/(1+mp.sqrt(q))**2),
            'passed':bool(passed)}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dps',type=int,default=90)
    parser.add_argument('--output',default=str(Path(__file__).with_name('generalization_checks_results.json')))
    args=parser.parse_args()
    if args.dps<70:parser.error('--dps must be at least 70')
    mp.mp.dps=args.dps;tol=mp.mpf(10)**(-args.dps+25)
    exact=exact_algebra_checks()
    k3=[]
    for q in K3_CASES:
        record=k3_case(q,tol);k3.append(record)
        print('k3',q,'slope',record['error_first_correction_quotient'],
              'remainder/q²',record['error_remainder_div_q_squared'],flush=True)
    f1=[f1_case(q,tol) for q in F1_CASES]
    summary={'all_passed':all(x['passed'] for x in k3+f1),
             'k3_infinite_Hankel_norm_cases':len(k3),'F1_radical_quartic_cases':len(f1),
             'maximum_k3_relative_Stein_residual':num(max(mp.mpf(x['maximum_relative_Stein_residual']) for x in k3)),
             'maximum_F1_value_difference':num(max(mp.mpf(x['value_difference']) for x in f1)),
             'maximum_absolute_k3_remainder_div_q_squared':num(max(abs(mp.mpf(x['error_remainder_div_q_squared'])) for x in k3)),
             'maximum_absolute_projection_gap_div_q_squared':num(max(abs(mp.mpf(x['error_minus_projection_div_q_squared'])) for x in k3))}
    payload={'status':'Numeric diagnostics, not interval certificates or proof of global optimality.',
             'environment':{'python':platform.python_version(),'mpmath':mp.__version__,'sympy':sp.__version__},
             'decimal_precision':args.dps,'comparison_tolerance':num(tol),
             'branch':{'lambda':'-6/7','w':['-3/14','sqrt(5)/35'],
                       'construction':'exact normal-response Stiefel branch'},
             'method':'Full 7-state infinite-error Stein Gramians and selected largest Hermitian eigenvalue; no word-horizon truncation.',
             'exact_algebra':exact,'summary':summary,'three_state_cases':k3,'rank_one_benchmark_cases':f1}
    Path(args.output).write_text(json.dumps(payload,indent=2)+'\n')
    print(json.dumps(summary,indent=2),flush=True)
    if not summary['all_passed']:raise SystemExit(1)
if __name__=='__main__':main()
