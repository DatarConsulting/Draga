"""Experimental integration of the independently accepted polynomial oracle."""
from spectral_oracle import PolynomialOracle, polynomial_lower
from reduced_codec import ReducedCodec, positive_power_lower
from fractions import Fraction as F
from certificate import ceil,ceil_log2


class SpectralCodec(ReducedCodec):
    def __init__(self,model):
        self.model=model
        self.oracle=PolynomialOracle(model)
        self.sigma=self.oracle.sigma
        self.n=self.oracle.certificate['candidate_horizon']
        self.certificate=self.oracle.certificate
        self.old_B=ceil(self.sigma*model.p)+model.cert.offset_bits
        self.B,self.capacity_certificate=self.compile_capacity()

    def probability(self,z):return self.oracle.branch(self.model.dual(z))

    def compile_capacity(self,grid_intervals=16):
        model=self.model;H=model.cert;L=model.depth_cap;eps=F(1,1<<model.p)
        l,r=H.J
        def z(t,x):return t*x+(1-t)*(1-x)
        geometric=min(4*z(t,l)*z(t,r) for t in H.Jdual)
        interior=min(F(1),geometric*(F(1,2)-2*model.point_error_bound/eps)/(r-l))
        assert interior>0
        pmin,pq=positive_power_lower(model,interior,self.sigma)
        root_x=(F(1,2)-self.oracle.center)/self.oracle.v
        root_P=sum(c*root_x**j for j,c in enumerate(self.oracle.coeffs))
        vmin=polynomial_lower(self.oracle.coeffs)/root_P
        mass=pmin*vmin*F(31,32)
        assert 0<mass<=1
        J0,J1=H.Jdual;raw=F(1,2);maxq=0
        for j in range(grid_intervals+1):
            t=J0+(J1-J0)*F(j,grid_intervals)
            lo,hi,q=self.oracle.enclosure(t,H.beta/4096)
            raw=min(raw,lo,1-hi);maxq=max(maxq,q)
        interpolation=H.K*(J1-J0)/(4*grid_intervals)
        exact_beta=max(H.beta,raw-interpolation)
        beta_hat=exact_beta-H.beta/(128*L)
        assert 0<beta_hat<=F(1,2)
        threshold=2/(beta_hat*mass)
        offset=ceil_log2(threshold)
        if (1<<offset)<=threshold:offset+=1
        B=ceil(self.sigma*model.p)+offset
        return B,dict(offset_bits=offset,bits=B,shape_lower=str(vmin),
            internal_mass_coefficient_lower=str(mass),rounded_branch_min_lower=str(beta_hat),
            required_root_factor=str(threshold),worst_capacity_loss=str((F(3,2)-beta_hat)/beta_hat),
            word_enumeration=False,tree_enumeration=False,
            scope='Instance certificate; universal finite-precision P4 remains open')
