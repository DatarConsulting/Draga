#!/usr/bin/env python3
import sympy as S, json
from pathlib import Path
x,s,a=S.symbols('x s a',real=True)
F0=s+(1-2*s)*(1-a)*x/(a*(1-x)+(1-a)*x)
F1=s+(1-2*s)*a*x/((1-a)*(1-x)+a*x)
T=S.Matrix([[1-s,s],[s,1-s]])
A0=S.diag(a,1-a)*T/2
alpha=S.Matrix([[S.Rational(1,2),S.Rational(1,2)]]); beta=alpha.T
p0=(alpha*A0*beta)[0]; p00=(alpha*A0*A0*beta)[0]
minor=S.factor(S.Rational(1,2)*p00-p0**2)
assert S.simplify(minor-(2*a-1)**2*(1-2*s)/64)==0
assert S.simplify(S.diff(F0,a)+(1-2*s)*x*(1-x)/(a*(1-x)+(1-a)*x)**2)==0
assert S.simplify(S.diff(F1,a)-(1-2*s)*x*(1-x)/((1-a)*(1-x)+a*x)**2)==0
special=[S.factor(z.subs({s:S.Rational(1,3),a:S.Rational(3,4)})) for z in (F0,F1)]
assert S.simplify(S.diff(special[0],x)-1/(3-2*x)**2)==0
assert S.simplify(S.diff(special[1],x)-1/(1+2*x)**2)==0
source={'passed':True,'general_rank_minor':str(minor),'baseline_maps':list(map(str,special)),
 'parameter_derivative_bound':'18/25','parameter_range_s':['1/4','5/12'],
 'parameter_range_a':['2/3','5/6'],'belief_interval':['3/8','5/8'],
 'robust_prediction_factor_stream':'229/132','hankel_hs_squared_family_bound':'144/361'}
(Path(__file__).resolve().parents[1]/'results/ANALYTIC.json').write_text(json.dumps(source,indent=2)+'\n')
print(json.dumps(source))
