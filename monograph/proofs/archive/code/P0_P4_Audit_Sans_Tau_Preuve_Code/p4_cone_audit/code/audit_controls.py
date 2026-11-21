"""Negative controls, fallback certificate and adaptive precision stress."""
from fractions import Fraction as F
from pathlib import Path
import json,time
from certificate import certify,hmm,CertificateError
from cases import CASES
from model import Model
from codec import Codec
from oracle import ratio,certify_oracle_grid,precision_budget
ROOT=Path(__file__).resolve().parents[1]

def main():
    rejections=[]
    fixtures=[('negative_determinant',hmm('3/5','3/5','1/4','3/4'),(F(3,8),F(5,8)),'positive_determinants'),
              ('overlap',hmm('1/3','1/3','1/2','1/2'),(F(3,8),F(5,8)),'branch_separation'),
              ('primal_expansion',hmm('1/20','1/20','1/4','3/4'),(F(1,100),F(99,100)),'primal_contraction')]
    for label,ms,J,expected in fixtures:
        try:certify(ms,J=J);raise AssertionError('must reject')
        except CertificateError as ex:
            assert ex.condition==expected
            rejections.append({'name':label,'condition':ex.condition})
    cert=certify(CASES[0][1],J=CASES[0][2]);model=Model.from_rationals(4,cert);codec=Codec(model)
    t=time.perf_counter();fallback=certify_oracle_grid(model,codec.n,codec.sigma)
    assert fallback is not None
    # Precision must actually double, with an independent 100-digit reference.
    import mpmath as mp
    from check_oracle_reference import V,real
    mp.mp.dps=100;s=F(3,4);tol=F(1,1<<120)
    box,c=ratio(model,F(3,8),3,s,tol)
    ref=V(model,real(F(3,8)),4,real(s))/V(model,real(F(3,8)),3,real(s))
    assert c.q>48 and mp.mpf(box[0])/c.S<=ref<=mp.mpf(box[1])/c.S
    assert c.q<=precision_budget(model,4,0,s,tol)
    # A long source representation must not survive in the compiled coefficients.
    src=[list(x) for x in cert.center];delta=F(1,1<<4096)
    src[0][0]+=delta;src[0][2]-=delta
    compact=Model.from_rationals(8,cert,src)
    assert compact.k<128
    data={'passed':True,'structural_rejections':rejections,'fallback_certificate':fallback,
          'fallback_and_stress_seconds':time.perf_counter()-t,
          'adaptive_precision':{'requested_width_bits':120,'used_fraction_bits':c.q,
                                'explicit_cap':precision_budget(model,4,0,s,tol),'mpmath_digits':100},
          'long_input':{'denominator_bits':max(v.denominator.bit_length() for m in src for v in m),
                        'compiled_k':compact.k,'scope':'Input work and memory must be counted before quantization.'}}
    (ROOT/'results/CONTROLS.json').write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps({'passed':True,'grid_points':fallback['grid_intervals']+1,'used_precision':c.q}),flush=True)
if __name__=='__main__':main()
