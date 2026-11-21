"""Replay the prerequisite and selected new certificate using exact arithmetic."""
import json,subprocess,sys
from fractions import Fraction
from pathlib import Path
import certify_relaxation
import verify_compact
ROOT=Path(__file__).resolve().parent

def main():
    manifest=json.loads((ROOT/'release_manifest.json').read_text())
    subprocess.run([sys.executable,str(ROOT/'prior_v08'/'certify_sos.py')],check=True)
    verify_compact.main()
    new=certify_relaxation.verify(ROOT/manifest['selected_certificate'])
    assert new['positive_lower_bound'] and new['improves_v08']
    lower=Fraction(new['four_decimal_endpoint']);upper=Fraction(130563,100000)
    assert str(lower)==manifest['lower']
    report={'version':'0.9.0','status':'EXACT_PASS','lower':str(lower),'upper':str(upper),
            'width':str(upper-lower),'closed_to_0_01':upper-lower<=Fraction(1,100),
            'scope':'ALL_BOUNDED_REAL_HANKEL_COMPETITORS_OF_RANK_AT_MOST_TWO',
            'new_certificate':manifest['selected_certificate'],
            'prior_certificate_replayed':True,'floating_point_proof_steps':False,
            'formal_proof_assistant':'NOT_USED','novelty':'NOT_ESTABLISHED'}
    (ROOT/'verification_report.json').write_text(json.dumps(report,indent=2))
    print('GLOBAL INTERVAL',str(lower),str(upper),'width',str(upper-lower),
          'closed_to_0.01',report['closed_to_0_01'],flush=True)

if __name__=='__main__':main()
