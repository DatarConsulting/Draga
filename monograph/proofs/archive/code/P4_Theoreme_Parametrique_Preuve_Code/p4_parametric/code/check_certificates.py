"""Sufficient-condition rejections and an admissibility boundary scan."""
from fractions import Fraction as F
from pathlib import Path
import json
from certificate import certify,hmm,CertificateError
from cases import CASES
ROOT=Path(__file__).resolve().parents[1]

def attempt(name,ms,J=(F(3,8),F(5,8)),Jdual=(F(1,8),F(7,8)),expected=None):
    try:
        c=certify(ms,J,Jdual)
        row={'name':name,'accepted':True,'tau':str(c.tau),'tau_decimal':float(c.tau)}
        assert expected is None
    except CertificateError as e:
        row={'name':name,'accepted':False,'condition':e.condition,'details':str(e.details)}
        assert e.condition==expected,(name,row,expected)
    return row

def main():
    bad=[[F(x) for x in m] for m in CASES[0][1]];bad[0][0]+=F(1,100)
    zero=[[F(x) for x in m] for m in CASES[0][1]];zero[0][0]=0;zero[0][2]+=F(1,6)
    unobs=((F(1,10),F(1,20),F(2,5),F(9,20)),(F(9,20),F(2,5),F(1,20),F(1,10)))
    rows=[
      attempt('normalization',bad,expected='stochastic_columns'),
      attempt('zero_entry',zero,expected='strict_positivity'),
      attempt('negative_determinant',hmm('3/5','3/5','1/4','3/4'),expected='positive_determinants'),
      attempt('interval_too_narrow',CASES[3][1],expected='primal_invariance'),
      attempt('weak_contraction',hmm('1/20','1/20','1/4','3/4'),J=(F(1,100),F(99,100)),expected='primal_contraction'),
      attempt('overlap',hmm('1/3','1/3','1/2','1/2'),expected='branch_separation'),
      attempt('unobservable',unobs,J=(F(1,20),F(19,20)),expected='observability'),
      attempt('dual_interval_too_narrow',CASES[0][1],Jdual=(F(49,100),F(51,100)),expected='dual_invariance'),
      attempt('transfer_not_contractive',hmm('1/4','1/4','6/25','19/25'),J=(F(1,4),F(3,4)),expected='transfer_contraction'),
    ]
    # Isolate the transfer contraction boundary for fixed interval suggestions.
    low,high=F(1,4),F(19,50)
    def status(s):
        try:return certify(hmm(s,s,F(6,25),F(19,25)),J=(F(1,4),F(3,4)))
        except CertificateError as e:
            assert e.condition=='transfer_contraction';return None
    assert status(low) is None and status(high) is not None
    for _ in range(24):
        mid=(low+high)/2
        if status(mid) is None:low=mid
        else:high=mid
    near=status(high)
    data={'passed':True,'rejection_tests':rows,'boundary':{'s_rejected':str(low),'s_accepted':str(high),
            'interval_width':str(high-low),'accepted_tau':str(near.tau),'accepted_tau_decimal':float(near.tau),
            'meaning':'Boundary of this sufficient cone certificate, not an impossibility boundary for compression.'}}
    (ROOT/'results/CERTIFICATE_TESTS.json').write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps({'passed':True,'rejections':len(rows),'boundary_tau':float(near.tau)}))
if __name__=='__main__':main()
