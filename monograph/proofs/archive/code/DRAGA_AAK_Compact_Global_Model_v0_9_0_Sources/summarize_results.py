"""Select only an exactly checked lower certificate and update display data."""
import json
from fractions import Fraction
from pathlib import Path
ROOT=Path(__file__).resolve().parent

def main():
    candidates=[]
    for p in ROOT.glob('reconstructed_certificate*_verification.json'):
        j=json.loads(p.read_text())
        if j['positive_lower_bound'] and j['improves_v08']:
            candidates.append((Fraction(j['four_decimal_endpoint']),p,j))
    lower,path,cert=max(candidates,key=lambda x:x[0]);upper=Fraction(130563,100000)
    chosen=path.name.replace('_verification.json','.json')
    (ROOT/'release_manifest.json').write_text(json.dumps({'selected_certificate':chosen,'lower':str(lower)},indent=2))
    def fr(value,digits):return f'{float(value):.{digits}f}'.replace('.','{,}')
    vals={'NewLower':fr(lower,4),'NewWidth':fr(upper-lower,5),
          'CertReference':fr(Fraction(cert['reference']),6),'CertLoss':fr(Fraction(cert['residual_budget']),12)}
    for prefix,suffix in [('Weak',''),('Strong','_stein'),('Trace','_stein_trace'),('Refined','_stein_trace_refined')]:
        p=ROOT/('relaxation_discovery'+suffix+'.json');j=json.loads(p.read_text()) if p.exists() else {}
        if 'objective' in j:
            vals.update({prefix+'Objective':f"{j['objective']:.9f}",prefix+'Iterations':str(j['iterations']),prefix+'Seconds':f"{j['wall_seconds']:.1f}"})
        else:vals.update({prefix+'Objective':'en cours',prefix+'Iterations':'--',prefix+'Seconds':'--'})
    vals['RunConclusion']='Les nombres ci-dessus servent à la découverte. Le certificat rationnel sélectionné est détaillé à la page suivante.'
    (ROOT/'run_values.tex').write_text('\n'.join('\\newcommand{\\'+k+'}{'+v+'}' for k,v in vals.items())+'\n')
    print('SELECTED',chosen,str(lower),'gap',str(upper-lower),flush=True)

if __name__=='__main__':main()
