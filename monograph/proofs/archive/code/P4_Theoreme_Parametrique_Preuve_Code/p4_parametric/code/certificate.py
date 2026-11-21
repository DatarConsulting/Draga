"""Rational sufficient certificates for positive binary 2x2 instruments.

Every bound covers the entire coefficient box. Stochastic normalization is
required for the actual source and compiled model, not for every box point.
Failure means not certified by these sufficient conditions.
"""
from dataclasses import dataclass, asdict
from fractions import Fraction as F

def ceil(x):
    x=F(x); return -((-x.numerator)//x.denominator)

def ceil_log2(x):
    x=F(x); n=max(0,x.numerator.bit_length()-x.denominator.bit_length())
    while (1<<n)<x:n+=1
    return n

def mat(m):return tuple(F(x) for x in m)
def det(m):a,b,c,d=m;return a*d-b*c
def point(m,x):
    a,b,c,d=m;return (a*x+b*(1-x))/((a+c)*x+(b+d)*(1-x))
def transpose(m):a,b,c,d=m;return a,c,b,d

class CertificateError(ValueError):
    def __init__(self,condition,details):
        self.condition=condition;self.details=details
        super().__init__(condition+': '+str(details))

def check(test,name,details):
    if not test:raise CertificateError(name,details)

@dataclass(frozen=True)
class Certificate:
    center:tuple
    lower:tuple
    upper:tuple
    radius:F
    J:tuple
    Jdual:tuple
    m:F
    rho:F
    gm:F
    grho:F
    Ddual:F
    K:F
    tau:F
    sigma_lo:F
    sigma_hi:F
    separation:F
    observation_min:F
    observation_max:F
    invariant_margin:F
    dual_invariant_margin:F
    amin:F
    amax:F
    detmin:F
    nu:F
    R:F
    distortion:F
    cdelta:F
    beta:F
    mass_constant:F
    offset_bits:int
    h:int
    parameter_TV_constant:F
    pressure_stability_constant:F

    def json(self):
        def cv(x):
            if isinstance(x,F):return str(x)
            if isinstance(x,dict):return {k:cv(v) for k,v in x.items()}
            if isinstance(x,(list,tuple)):return [cv(v) for v in x]
            return x
        return cv(asdict(self))

def stochastic(ms):
    return all(sum(m[col]+m[col+2] for m in ms)==1 for col in (0,1))

def geometry(lower,upper,J):
    l,r=J;derivs=[];D=[];images=[];nus=[]
    for low,high in zip(lower,upper):
        a,b,c,d=low;A,B,C,E=high
        dl=a*d-B*C;du=A*E-b*c
        check(dl>0,'positive_determinants',str(dl))
        def den(m,x):return (m[0]+m[2])*x+(m[1]+m[3])*(1-x)
        denlo=min(den(low,l),den(low,r));denhi=max(den(high,l),den(high,r))
        nus.append(denlo)
        derivs.append((dl/denhi**2,du/denlo**2))
        slope=max(abs((A+C)-(b+d)),abs((a+c)-(B+E)))
        D.append(2*slope/denlo)
        image_low=(a*l+b*(1-l))/den(high,l)
        image_high=(A*r+B*(1-r))/den(low,r)
        images.append((image_low,image_high))
    return min(x[0] for x in derivs),max(x[1] for x in derivs),max(D),tuple(images),min(nus)

def certify(matrices,J=(F(3,8),F(5,8)),Jdual=(F(1,8),F(7,8)),radius=F(1,65536)):
    ms=tuple(mat(m) for m in matrices);J=tuple(map(F,J));Jdual=tuple(map(F,Jdual));radius=F(radius)
    check(len(ms)==2 and all(len(m)==4 for m in ms),'shape','two 2x2 matrices')
    check(radius>=0,'coefficient_radius',radius)
    check(stochastic(ms),'stochastic_columns','sum of each column across both matrices must be 1')
    for name,I in [('primal_interval',J),('dual_interval',Jdual)]:
        check(0<I[0]<F(1,2)<I[1]<1,name,str(I))
    lower=tuple(tuple(x-radius for x in m) for m in ms)
    upper=tuple(tuple(x+radius for x in m) for m in ms)
    amin=min(x for m in lower for x in m);amax=max(x for m in upper for x in m)
    check(amin>0,'strict_positivity',str(amin))
    detmin=min(lo[0]*lo[3]-hi[1]*hi[2] for lo,hi in zip(lower,upper))
    check(detmin>0,'positive_determinants',str(detmin))
    m,rho,D,images,nu=geometry(lower,upper,J)
    margin=min(min(im[0]-J[0],J[1]-im[1]) for im in images)
    check(margin>0,'primal_invariance',str(margin))
    check(rho<1,'primal_contraction',str(rho))
    gap=max(images[1][0]-images[0][1],images[0][0]-images[1][1])
    check(gap>0,'branch_separation',str(gap))
    sl=lower[0][0]+lower[0][2]-upper[0][1]-upper[0][3]
    su=upper[0][0]+upper[0][2]-lower[0][1]-lower[0][3]
    check(sl>0 or su<0,'observability',str((sl,su)))
    omin=min(abs(sl),abs(su));omax=max(abs(sl),abs(su))
    gm,grho,Ddual,gimages,_=geometry(tuple(map(transpose,lower)),tuple(map(transpose,upper)),Jdual)
    gmargin=min(min(im[0]-Jdual[0],Jdual[1]-im[1]) for im in gimages)
    check(gmargin>0,'dual_invariance',str(gmargin))
    check(grho<1,'dual_contraction',str(grho))
    # A finite rational search is a convenience, not a completeness claim.
    lows=[F(j,8) for j in range(1,33) if 2**8*gm**j>1]
    highs=[F(j,8) for j in range(1,33) if 2**8*grho**j<1]
    check(bool(lows and highs),'pressure_bracket','no bracket on the j/8, 1<=j<=32 search grid')
    slo=max(lows);shi=min(highs)
    check(slo<shi,'pressure_bracket',str((slo,shi)))
    K=F(65,64)*shi*Ddual/(1-grho)
    ell=Jdual[1]-Jdual[0]
    tau=grho+ell*K/2
    check(tau<1,'transfer_contraction',str(tau))
    R=amax/amin
    l,r=J;ellj=r-l;aj=min(l,1-r)
    distortion=max(r/l,(1-l)/(1-r))**2
    cdelta=min(F(1),aj**2*m/(distortion*ellj))
    ratio=(grho/gm)**ceil(shi)*3**ceil(K*ell)
    beta=1/(1+ratio)
    mass=cdelta**ceil(shi)/3**ceil(K*ell+1)
    offset=ceil_log2((9+4/beta)/mass)
    h=1
    while rho**h>F(1,2):h+=1
    Cpar=2+3*omax/(nu*(1-rho))
    Cpressure=shi*(4*amax/detmin+2/amin)/(1-grho)
    return Certificate(ms,lower,upper,radius,J,Jdual,m,rho,gm,grho,Ddual,K,tau,slo,shi,
                       gap,omin,omax,margin,gmargin,amin,amax,detmin,nu,R,distortion,cdelta,
                       beta,mass,offset,h,Cpar,Cpressure)

def hmm(s0,s1,e0,e1):
    s0,s1,e0,e1=map(F,(s0,s1,e0,e1))
    def one(e,f):return ((1-s0)*e,s1*f,s0*e,(1-s1)*f)
    return one(e0,e1),one(1-e0,1-e1)
