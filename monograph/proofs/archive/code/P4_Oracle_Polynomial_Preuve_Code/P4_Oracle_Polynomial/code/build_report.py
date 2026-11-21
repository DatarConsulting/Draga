from pathlib import Path
from fractions import Fraction as F
import json
ROOT=Path(__file__).resolve().parents[1]
conf=json.loads((ROOT/'results/CONFIGURATIONS.json').read_text())
graphs=json.loads((ROOT/'results/GRAPHS.json').read_text())
bench=json.loads((ROOT/'results/BENCHMARKS.json').read_text())
lines=[r'\paragraph{Contrôles.}',
 f"Les {len(conf['rows'])} configurations S/A, pour $p=4,8,12,16,32,64,128$, passent.",
 f"Les six graphes complets à $p=4,8,12$ totalisent {graphs['states']} états et {graphs['transitions']} transitions.",
 f"Le maximum TV/$2^{{-p}}$ sur ces graphes vaut {max(r['max_TV_over_epsilon'] for r in graphs['rows']):.6f}.",
 'Les 338 entrées de matrice comparées à une formule rationnelle indépendante sont toutes encadrées.',
 'Les 14 flux de 64 observations passent également; ce sont des tests finis.',
 r'\begin{center}\begin{tabular}{crrrrrr}\toprule',
 r'Modèle & $p$ & Bits état & Degré & Horizon & $e$ certifié & TV/$\varepsilon$\\\midrule']
for x in conf['rows']:
 c=x['oracle']['certificate']
 mantissa,exponent=f"{float(F(c['ratio_deviation_bound'])):.2e}".split('e')
 sci=mantissa+r'\times10^{'+str(int(exponent))+'}'
 lines.append(f"{x['model']} & {x['p']} & {x['bits']} & {c['degree']} & {c['candidate_horizon']} & ${sci}$ & {x['max_TV_over_epsilon']:.6f}"+r'\\')
lines += [r'\bottomrule\end{tabular}\end{center}',
 'Toutes les instances acceptées ont un polynôme de degré 32. Son format compact théorique compte 1666 bits, soit 209 octets arrondis, hors autres paramètres et objets Python.',
 r'Le décalage vaut trois bits au-dessus de $\lceil\sigma p\rceil$ dans les 14 cas; il ne signifie pas trois bits au-dessus de l’optimum.',
 r'\paragraph{Mesures de temps.}',
 'CPU AMD EPYC 9V74, CPython, conteneur partagé sans affinité CPU. Médiane de trois essais de 32 observations chacun; avance et sortie incluses, vérificateur exact exclu. Chaque essai repart de la racine : il ne mesure pas la profondeur maximale pour les grands $p$.',
 r'\begin{center}\begin{tabular}{crrrrr}\toprule',
 r'Modèle & $p$ & Ancien, ms/pas & Nouveau, ms/pas & Rapport & Bits communs\\\midrule']
for x in bench['rows']:
 old=x['enumeration'];new=x['polynomial']
 assert old['bits']==new['bits']
 lines.append(f"{x['model']} & {x['p']} & {old['step_output_ms']:.2f} & {new['step_output_ms']:.2f} & {x['observed_speedup']:.2f} & {new['bits']}"+r'\\')
lines += [r'\bottomrule\end{tabular}\end{center}',
 r'\begin{center}\begin{tabular}{crrrrr}\toprule',
 r'Modèle & $p$ & Compil. ancienne, s & Compil. nouvelle, s & $n$ ancien & $n$ nouveau\\\midrule']
for x in bench['rows']:
 old=x['enumeration'];new=x['polynomial']
 lines.append(f"{x['model']} & {x['p']} & {old['compile_seconds']:.3f} & {new['compile_seconds']:.3f} & {old['horizon']} & {new['horizon']}"+r'\\')
lines += [r'\bottomrule\end{tabular}\end{center}',
 'Le nombre d’essais et les flux sont identiques entre les deux codecs pour chaque comparaison. Les temps de compilation sont des observations uniques; les temps d’exécution sont les médianes indiquées. Les horizons de génération peuvent différer.']
if not bench['passed']:
 lines += [r'\textbf{Document provisoire : campagne de temps encore en cours.}']
(ROOT/'report/results_table.tex').write_text('\n'.join(lines)+'\n')
print('Tables generated; benchmark completed:',bench['passed'])
