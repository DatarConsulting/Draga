#!/usr/bin/env python3
"""Reproducible finite-Stein diagnostics of the binary boundary theorem.

Dependencies: numpy, scipy. Run:
    python hankel_boundary_checks.py --output hankel_boundary_checks_results.json

All recorded computations are floating-point diagnostics, not proof certificates,
global optimum computations, or finite-section approximations of Hankel norms.
The complex lower-bound theorem and Hardy bridge require the analytic proofs.
"""
from __future__ import annotations
import argparse
import json
import platform
from pathlib import Path

import numpy as np
import scipy
from scipy.linalg import block_diag, solve


def adj(x):
    return x.conj().T


def herm(x):
    return (x+adj(x))/2


def psd_sqrt(x):
    eig,u=np.linalg.eigh(herm(x))
    scale=max(1.,float(np.max(np.abs(eig))))
    if eig.min() < -2e-9*scale:
        raise AssertionError(f"Gramian is not positive within tolerance: {eig.min()}")
    return (u*np.sqrt(np.maximum(eig,0.))) @ adj(u)


def stein(letters,source,side):
    n=source.shape[0]
    dtype=np.result_type(source.dtype,*[x.dtype for x in letters])
    kernel=np.zeros((n*n,n*n),dtype=dtype)
    for a in letters:
        kernel += np.kron(a.conj(),a) if side=='row' else np.kron(a.T,adj(a))
    vec=solve(np.eye(n*n,dtype=dtype)-kernel,
              source.reshape(-1,order='F'),assume_a='gen')
    result=herm(vec.reshape((n,n),order='F'))
    update=sum((a@result@adj(a) if side=='row' else adj(a)@result@a)
               for a in letters)
    residual=float(np.linalg.norm(result-source-update,'fro')/
                   max(1.,np.linalg.norm(result,'fro'),np.linalg.norm(source,'fro')))
    return result,residual


def norms(letters,b,c):
    """Exact-in-principle infinite Hankel norm and coefficient H2 norm."""
    reach,rr=stein(letters,b@adj(b),'row')
    obs,ro=stein(letters,adj(c)@c,'column')
    root=psd_sqrt(obs)
    top=float(np.linalg.eigvalsh(herm(root@reach@root)).max())
    energy=float(np.real((c@reach@adj(c))[0,0]))
    return np.sqrt(max(0.,top)),max(0.,energy),max(rr,ro)


def monomial_realization(m):
    """Standard-order alpha A_w beta for the one-letter series z^m."""
    shift=np.zeros((m+1,m+1))
    for j in range(m): shift[j,j+1]=1.
    b=np.zeros((m+1,1)); b[m,0]=1.
    c=np.zeros((1,m+1)); c[0,0]=1.
    return shift,b,c


def rational_realization(numerator,denominator,m):
    """h=N/P, deg N<m, denominator constant nonzero; h_j=c D^j b.

    The state is (h_j,...,h_{j+m-1}) and D is its recurrence companion.
    """
    if m==0:
        return np.zeros((0,0)),np.zeros((0,1)),np.zeros((1,0))
    den=np.asarray(denominator,dtype=float)/denominator[0]
    num=np.asarray(numerator,dtype=float)/denominator[0]
    initial=np.zeros(m)
    for j in range(m):
        initial[j]=((num[j] if j<len(num) else 0.)-
                    sum(den[t]*initial[j-t] for t in range(1,min(j,m)+1)))
    D=np.zeros((m,m))
    for j in range(m-1): D[j,j+1]=1.
    D[-1,:]=-den[m:0:-1]
    b=initial.reshape(m,1)
    c=np.zeros((1,m)); c[0,0]=1.
    return D,b,c


def polynomial_case(m):
    sqrt33=np.sqrt(33.)
    raw=[np.array([1.]),np.array([1.,1/np.sqrt(2.)]),np.array([1.,1.,.5]),
         np.array([1.,np.sqrt(3+sqrt33/3)/2,(1+sqrt33)/8,
                   np.sqrt(15-sqrt33)/8])][m]
    C=[1.,2/np.sqrt(3.),4/3,np.sqrt(2*(1103+33*sqrt33)/1153)][m]
    p=raw/np.linalg.norm(raw)
    reverse=p[::-1].copy()
    square=np.convolve(p,p)
    el_error=float(np.max(np.abs(square[:m+1]-C*reverse)))
    cubic=np.convolve(square,p)
    cubic_coefficient=float(cubic[m])
    shifted=np.concatenate([np.zeros(m),p])
    numerator_full=shifted-np.convolve(reverse,reverse)/C
    high_cancellation=float(np.max(np.abs(numerator_full[m:])))
    numerator=numerator_full[:m]
    D,b,c=rational_realization(numerator,p,m)
    target,tb,tc=monomial_realization(m)
    if m:
        diff_letters=[block_diag(target,D)]
        diff_b=np.concatenate([tb,b],axis=0)
        diff_c=np.concatenate([tc,-c],axis=1)
        roots=np.roots(p[::-1])
        pole_radius=float(np.max(np.abs(np.linalg.eigvals(D))))
        denominator_min_root_modulus=float(np.min(np.abs(roots)))
    else:
        diff_letters=[target]; diff_b=tb; diff_c=tc
        pole_radius=0.; denominator_min_root_modulus=None
    hankel_error,h2_error_squared,residual=norms(diff_letters,diff_b,diff_c)
    reach,rr=stein(diff_letters,diff_b@adj(diff_b),'row')
    obs,ro=stein(diff_letters,adj(diff_c)@diff_c,'column')
    correlation_state=diff_letters[0]@reach@adj(diff_c)
    T_stein=float(np.real((adj(correlation_state)@obs@correlation_state)[0,0]))
    correlations=[float(np.dot(p[j:],p[:-j])) for j in range(1,m+1)]
    correlation_energy=sum(x*x for x in correlations)
    T_polynomial=correlation_energy/C**4
    expected_h2=1/(C*C)
    tolerance=2e-10
    record={'m':m,'k':m+1,'P_coefficients':p.tolist(),
            'P_normalization_squared':float(p@p),
            'C':C,'inverse_C':1/C,'cubic_coefficient':cubic_coefficient,
            'EL_prefix_max_abs_error':el_error,
            'h_numerator':numerator.tolist(),'h_denominator':p.tolist(),
            'high_degree_cancellation_max_abs_error':high_cancellation,
            'h_state_dimension':m,'h_transition_spectral_radius':pole_radius,
            'denominator_minimum_root_modulus':denominator_min_root_modulus,
            'H2_error_squared':h2_error_squared,'expected_H2_error_squared':expected_h2,
            'one_letter_hankel_error_norm':hankel_error,
            'positive_P_autocorrelations':correlations,
            'P_H4_norm_fourth_power':1+2*correlation_energy,
            'T_autocorrelation_energy_Stein':T_stein,
            'T_autocorrelation_energy_polynomial':T_polynomial,
            'explicit_squared_error_q2_coefficient':T_polynomial/expected_h2,
            'stein_relative_residual':residual,
            'passed':bool(el_error<tolerance and abs(cubic_coefficient-C)<tolerance
                          and high_cancellation<tolerance
                          and abs(h2_error_squared-expected_h2)<tolerance
                          and abs(T_stein-T_polynomial)<tolerance
                          and pole_radius<1.)}
    return record,(D,b,c)


def word_target(r,k):
    """f(a^n b^k)=q r^n, alpha=q e0*, beta=e_k, standard order."""
    q=1-r*r
    aa=np.zeros((k+1,k+1)); aa[0,0]=r
    ab=np.zeros((k+1,k+1))
    for j in range(k): ab[j,j+1]=1.
    b=np.zeros((k+1,1)); b[k,0]=1.
    c=np.zeros((1,k+1)); c[0,0]=q
    return [aa,ab],b,c


def word_competitor(r,h_realization):
    """g(a^n b^(j+1))=q r^n h_j with at most m+1 states."""
    D,bh,ch=h_realization
    m=D.shape[0]; q=1-r*r
    aa=np.zeros((m+1,m+1)); aa[0,0]=r
    ab=np.zeros((m+1,m+1))
    if m:
        ab[0,1:]=ch[0]
        ab[1:,1:]=D
    b=np.concatenate([np.zeros((1,1)),bh],axis=0)
    c=np.zeros((1,m+1)); c[0,0]=q
    return [aa,ab],b,c


def word_check(q,poly,h_realization):
    k=poly['k']; r=np.sqrt(1-q)
    fa,fb,fc=word_target(r,k)
    ga,gb,gc=word_competitor(r,h_realization)
    ea=[block_diag(a,b) for a,b in zip(fa,ga)]
    eb=np.concatenate([fb,gb],axis=0)
    ec=np.concatenate([fc,-gc],axis=1)
    error,_,residual=norms(ea,eb,ec)
    E2=poly['expected_H2_error_squared']
    upper_squared=E2+q*poly['one_letter_hankel_error_norm']**2
    record={'k':k,'m':k-1,'q':q,'r':r,
            'target_state_dimension':k+1,'competitor_state_dimension':k,
            'actual_infinite_hankel_error':float(error),
            'actual_error_squared':float(error*error),
            'boundary_limit':poly['inverse_C'],
            'actual_minus_boundary_limit':float(error-poly['inverse_C']),
            'analytic_squared_upper_bound':float(upper_squared),
            'stein_relative_residual':residual,
            'passed':bool(error*error<=upper_squared+2e-9)}
    Hnorm2=poly['one_letter_hankel_error_norm']**2
    T=poly['T_autocorrelation_energy_polynomial']
    if q*Hnorm2<=E2/2:
        quadratic_lower=E2+q*q*T/E2
        quadratic_upper=quadratic_lower+2*q**3*Hnorm2*T/E2**2
        quadratic_passed=bool(quadratic_lower-2e-10<=error*error<=quadratic_upper+2e-10)
        record.update({'all_k_quadratic_rate_squared_lower':quadratic_lower,
                       'all_k_quadratic_rate_squared_upper':quadratic_upper,
                       'all_k_quadratic_rate_passed':quadratic_passed,
                       'observed_squared_error_q2_ratio':float((error*error-E2)/q**2)})
        record['passed']=record['passed'] and quadratic_passed
    if k==2 and q<=1/8:
        sharp_upper=.75+q*q/6+4*q**3/3
        competitor_rayleigh_lower=.75+q*q/6
        record.update({'k2_sharp_squared_upper_bound':sharp_upper,
                       'k2_competitor_rayleigh_squared_lower':competitor_rayleigh_lower,
                       'k2_sharp_bound_passed':bool(competitor_rayleigh_lower-2e-10<=error*error
                                                  <=sharp_upper+2e-10)})
        record['passed']=record['passed'] and record['k2_sharp_bound_passed']
    return record


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',default='hankel_boundary_checks_results.json')
    args=parser.parse_args()
    polys=[]; words=[]
    for m in range(4):
        rec,h=polynomial_case(m)
        polys.append(rec)
        for q in (.2,.05,.01,.003,.001): words.append(word_check(q,rec,h))
    summary={'all_passed':all(x['passed'] for x in polys+words),
             'polynomial_cases':len(polys),'infinite_word_hankel_checks':len(words),
             'k2_sharp_checks':sum('k2_sharp_bound_passed' in x for x in words),
             'all_k_quadratic_rate_checks':sum('all_k_quadratic_rate_passed' in x for x in words),
             'max_EL_prefix_error':max(x['EL_prefix_max_abs_error'] for x in polys),
             'max_h_numerator_cancellation_error':max(x['high_degree_cancellation_max_abs_error'] for x in polys),
             'max_stein_relative_residual':max(x['stein_relative_residual'] for x in polys+words),
             'boundary_constants':{str(x['k']):x['inverse_C'] for x in polys}}
    data={'status':'floating-point diagnostics only; not global optimum certificates or proofs',
          'environment':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__},
          'method':'finite stable realizations and finite Stein Gramians; no word-horizon truncation',
          'convention':'standard realization f(w)=c A_w b; counting-measure word-Hankel norm',
          'interpretation':'actual errors belong to explicit admissible competitors; they are not claimed equal to delta_k at finite r',
          'summary':summary,'polynomial_cases':polys,'word_cases':words}
    Path(args.output).write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(summary,indent=2))
    if not summary['all_passed']: raise SystemExit('Diagnostic failure; inspect JSON.')


if __name__=='__main__': main()
