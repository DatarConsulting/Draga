from pathlib import Path
from fractions import Fraction as F
import json
ROOT=Path(__file__).resolve().parents[1]
def read(n):return json.loads((ROOT/f'results/{n}.json').read_text())
v=read('VERIFIED');pr=read('PRODUCTS');st=read('STREAM_MEMORY')
assert len(v['configs'])==27 and len(st['rows'])==4
def french_int(x):return f'{x:,}'.replace(',',r'\,')
macros={'NumGraphs':27,'NumStates':v['total_states'],'NumEdges':v['total_transitions'],
        'NumOutputs':v['total_outputs'],'MinIntegerWidth':min(z['minimum_integer_interval'] for z in v['configs']),
        'MaxTVRatio':f"{max(z['max_TV_over_epsilon_decimal'] for z in v['configs']):.6f}",
        'NumProductChecks':pr['total_exact_point_checks'],
        'ReducedExactBits':next(z['gcd_reduced_reference_max_bits'] for z in pr['rows'] if z['p']==256 and z['s']=='1/3'),
        'MaxStreamTVRatio':f"{max(z['max_TV_over_epsilon_decimal'] for z in st['rows']):.6f}"}
(ROOT/'report/numbers.tex').write_text(''.join('\\newcommand{\\'+k+'}{'+(french_int(val) if isinstance(val,int) else val)+'}\n' for k,val in macros.items()))
rows=[r'\begin{center}\small\begin{tabular}{llrrrrr}\toprule',
      r'$s$ & $a$ & $N_8$ & $N_{12}$ & $N_{16}$ & $B_{16}$ & TV/$\eps$ max.\\\midrule']
for s in ['8/25','1/3','17/50']:
 for a in ['37/50','3/4','19/25']:
  group=[next(z for z in v['configs'] if z['true_s']==s and z['true_a']==a and z['p']==p) for p in [8,12,16]]
  rows.append(' & '.join([s,a]+[str(z['states']) for z in group]+[str(group[-1]['bits']),f"{max(z['max_TV_over_epsilon_decimal'] for z in group):.6f}"])+r'\\')
rows.extend([r'\bottomrule\end{tabular}\end{center}'])
(ROOT/'report/graphs.tex').write_text('\n'.join(rows)+'\n')
rows=[r'\begin{center}\small\begin{tabular}{rrrrrrrr}\toprule',
      r'$p$ & $k$ & $q$ & $B$ & $n$ & Préparation (Kio) & En ligne (Kio) & Pas (ms)\\\midrule']
for z in st['rows']:
 rows.append(' & '.join([str(z[k]) for k in ['p','model_k','matrix_q','bits','horizon']]+[f"{z['preparation_tracemalloc_peak_bytes']/1024:.2f}",f"{z['model_plus_online_peak_bytes']/1024:.2f}",f"{z['median_advance_and_output_seconds']*1000:.2f}"])+r'\\')
rows.extend([r'\bottomrule\end{tabular}\end{center}'])
(ROOT/'report/memory.tex').write_text('\n'.join(rows)+'\n')
print(macros)
