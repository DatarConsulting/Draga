"""Polynomial transfer oracle: floating candidates, exact interval acceptance.

No auxiliary words or state tree are enumerated by this module. Acceptance
certifies the stored polynomial itself, not an unjustified exact jet of L^n 1.
Candidate search is intentionally finite; universal compiler termination is
NOT asserted. All runtime decisions use integers and rational numbers.
"""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'reference_code'))
from fractions import Fraction as F
from math import isqrt
from certificate import ceil, ceil_log2
from dyadic import Dyadic
from oracle import power, NeedPrecision, trajectory_step


def coordinates(model):
    l, h = model.cert.Jdual
    center, v = (l+h)/2, (h-l)/2
    ans = []
    for a, b, c, d in model.mats:
        den = (a+b)*center + (c+d)*(1-center)
        xi = v*(a+b-c-d)/den
        alpha = ((a*center+c*(1-center))/den-center)/v
        beta = ((a-c)*v/den-center*xi)/v
        kappa = F(a*d-b*c, den*den)
        ans.append((alpha, beta, xi, kappa))
    return center, v, ans


def disk_certificate(model):
    """Rational circle test: G(D_R) subset D_w, with 1 < w < R."""
    center, v, maps = coordinates(model)
    R = F(2)
    for _ in range(1024):
        if all(abs(xi)*R < 1 for _, _, xi, _ in maps):
            eta = max(abs((a+b*x)/(1+xi*x))
                      for a, b, xi, _ in maps for x in (-R, R))
            if eta < R:
                w = (1+R)/2
                while w <= eta:w = (w+R)/2
                A = ceil(2*model.cert.sigma_hi)
                # kappa < 1 at the center by the inherited certificate.
                assert all(0 < kap < 1 for _, _, _, kap in maps)
                W = F(ceil(sum((1-abs(xi)*R)**(-A) for _, _, xi, _ in maps)))
                return dict(R=R, w=w, eta=eta, W=W,
                            Lambda=W/(1-w/R), center=center, v=v)
        R = 1+(R-1)/2
    raise ValueError('No disk certificate found in the finite search budget')


def float_matrix(model, sigma, m):
    result = [[0.0]*(m+1) for _ in range(m+1)]
    for a, b, xi, kap in coordinates(model)[2]:
        a, b, xi, kap = map(float, (a,b,xi,kap))
        col = [kap**sigma]
        for k in range(m):
            col.append(col[-1]*(-xi)*(2*sigma+k)/(k+1))
        for j in range(m+1):
            for k in range(m+1): result[k][j] += col[k]
            if j < m:
                nxt = []
                for k in range(m+1):
                    nxt.append(a*col[k]+(b*col[k-1]-xi*nxt[k-1] if k else 0))
                col = nxt
    return result


def float_iterate(M, n):
    c = [1.0]+[0.0]*(len(M)-1)
    for _ in range(n):
        d = [sum(a*b for a,b in zip(row,c)) for row in M]
        c = [x/d[0] for x in d]
    ratio = sum(a*b for a,b in zip(M[0],c))
    return c, ratio


def candidate(model, n, m):
    lo, hi = map(float, (model.cert.sigma_lo,model.cert.sigma_hi))
    for _ in range(46):
        mid = (lo+hi)/2
        _, ratio = float_iterate(float_matrix(model,mid,m),n)
        if ratio > 1:lo = mid
        else:hi = mid
    sbits = 40
    sigma = F(round((lo+hi)/2*(1<<sbits)),1<<sbits)
    c, _ = float_iterate(float_matrix(model,float(sigma),m),n)
    # These are exact dyadic candidate coefficients. Their quality is verified
    # independently below, so the float search is outside the trusted core.
    coeffs = [F(round(x*(1<<48)),1<<48) for x in c]
    coeffs[0] = F(1)
    return sigma, coeffs


def interval_matrix(model, sigma, m, ctx):
    """Exact finite-section inclusion; O(m^2) interval operations.

    Col_{j+1}=(alpha+beta*x)/(1+xi*x) Col_j gives a triangular recurrence.
    The derivative weight is kappa^sigma (1+xi*x)^(-2 sigma).
    """
    result = [[(0,0) for _ in range(m+1)] for _ in range(m+1)]
    for a, b, xi, kap in coordinates(model)[2]:
        col = [power(ctx,ctx.rat(kap),sigma)]
        for k in range(m):
            col.append(ctx.mul(col[-1],ctx.rat(-xi*(2*sigma+k)/(k+1))))
        ai, bi, xi_i = ctx.rat(a),ctx.rat(b),ctx.rat(xi)
        for j in range(m+1):
            for k in range(m+1):result[k][j] = ctx.add(result[k][j],col[k])
            if j < m:
                nxt = []
                for k in range(m+1):
                    term = ctx.mul(ai,col[k])
                    if k:term = ctx.add(term,ctx.sub(ctx.mul(bi,col[k-1]),ctx.mul(xi_i,nxt[k-1])))
                    nxt.append(term)
                col = nxt
    return result


def polynomial_lower(coeffs):
    # x^(2j) in [0,1], x^(2j+1) in [-1,1].
    return coeffs[0]+sum(min(F(0),x) if k%2==0 else -abs(x)
                         for k,x in enumerate(coeffs[1:],1))


def verify_polynomial(model, sigma, coeffs, q=96):
    disk = disk_certificate(model)
    m = len(coeffs)-1
    ctx = Dyadic(q)
    try: M = interval_matrix(model,sigma,m,ctx)
    except NeedPrecision:return None
    xs = [ctx.rat(x) for x in coeffs]
    residual = F(0)
    for k,row in enumerate(M):
        y = (0,0)
        for a,b in zip(row,xs):y = ctx.add(y,ctx.mul(a,b))
        y = ctx.sub(y,xs[k])
        residual += F(max(abs(y[0]),abs(y[1])),ctx.S)
    norm = sum(abs(c)*disk['w']**j for j,c in enumerate(coeffs))
    tail = disk['W']*norm*disk['R']**(-m)/(disk['R']-1)
    lower = polynomial_lower(coeffs)
    if lower <= 0:return None
    dev = (residual+tail)/lower
    target = F(1,64*model.depth_cap)
    # Keeping the inherited log-Lipschitz envelope avoids silently changing
    # the branch interpolation and dual-coordinate precision budgets.
    deriv = sum(j*abs(c) for j,c in enumerate(coeffs))/(disk['v']*lower)
    if dev > target/8 or deriv > model.cert.K:return None
    return dict(method='rational_disk_polynomial_residual',degree=m,
                sigma=str(sigma),verification_precision=q,
                disk={k:str(v) for k,v in disk.items()},
                residual_coefficient_bound=str(residual),tail_bound=str(tail),
                polynomial_lower_bound=str(lower),log_derivative_bound=str(deriv),
                ratio_deviation_bound=str(dev),log_residual_bound=str(2*dev),
                target=str(target),word_enumeration=False,
                universal_compiler_termination_proved=False)


class PolynomialOracle:
    def __init__(self,model):
        self.model = model
        self.center,self.v,_ = coordinates(model)
        self.attempts = []
        for m in (32,48,64):
            for n in (4,8,16):
                sigma, coeffs = candidate(model,n,m)
                for q in (96,192):
                    cert = verify_polynomial(model,sigma,coeffs,q)
                    self.attempts.append(dict(degree=m,horizon=n,q=q,accepted=cert is not None))
                    if cert:
                        self.sigma,self.coeffs,self.certificate = sigma,tuple(coeffs),cert
                        cert['candidate_horizon'] = n
                        cert['candidate_generator'] = 'untrusted_float_with_exact_dyadic_output'
                        return
        raise ValueError('No accepted candidate in finite search; no false success')

    def evaluate(self,t,ctx):
        x = ctx.div(ctx.sub(t,ctx.rat(self.center)),ctx.rat(self.v))
        out = (0,0)
        for a in reversed(self.coeffs):out = ctx.add(ctx.mul(out,x),ctx.rat(a))
        return out

    def enclosure(self,t,tolerance):
        for q in (48,96,192,384,768,1536):
            ctx = Dyadic(q)
            try:
                vals=[]
                for b in (0,1):
                    x,g=trajectory_step(self.model,ctx.rat(t),b,ctx)
                    P=self.evaluate(x,ctx)
                    if P[0]<=0:raise NeedPrecision()
                    vals.append(ctx.mul(power(ctx,g,self.sigma),P))
                den=ctx.add(*vals)
                if den[0]<=0:raise NeedPrecision()
                r=ctx.div(vals[0],den)
                if F(r[1]-r[0],ctx.S)<=tolerance:
                    return F(r[0],ctx.S),F(r[1],ctx.S),q
            except NeedPrecision:pass
        raise ArithmeticError('Finite runtime precision budget exhausted')

    def branch(self,t):
        L=self.model.depth_cap;beta=self.model.cert.beta
        lo,hi,_=self.enclosure(t,beta/(512*L))
        Q=1<<ceil_log2(1024*L/beta)
        mid=(lo+hi)/2
        k=(2*mid.numerator*Q+mid.denominator)//(2*mid.denominator)
        assert beta/2<=F(k,Q)<=1-beta/2
        return k,Q

    def json(self):
        bits_per_coefficient=max(abs(int(x*(1<<48))).bit_length()+1 for x in self.coeffs)
        return dict(certificate=self.certificate,coefficients=list(map(str,self.coeffs)),
                    attempts=self.attempts,
                    compiled_polynomial_payload_bits=len(self.coeffs)*bits_per_coefficient+16,
                    coefficient_storage_fraction_bits=48,
                    payload_excludes_python_object_overhead=True)
