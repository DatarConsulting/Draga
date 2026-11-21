from pathlib import Path
from fractions import Fraction as F
import json

R=Path(__file__).resolve().parents[1]
data={name:json.loads((R/'results'/f'{name}.json').read_text())
      for name in ('COMPARISON','GRAPHS','STREAMS','CONTROLS','HYBRID')}
assert all(x['passed'] for x in data.values())
assert len(data['COMPARISON']['rows'])==16
assert all(x['new_offset']==3 for x in data['COMPARISON']['rows'])
lookup={(x['model'],x['p']):x for x in data['COMPARISON']['rows']}
table=['\\begin{tabular}{rrrrrr}\\toprule',
       '$p$ & Ancien S & Réduit S & Ancien A & Réduit A & Grille\\\\\\midrule']
for p in (4,8,12,16,24,64,128,256):
    s,a=lookup['S',p],lookup['A',p]
    table.append(f"{p} & {s['old_bits']} & {s['new_bits']} & {a['old_bits']} & {a['new_bits']} & {s['grid_bits']}\\\\")
table.append('\\bottomrule\\end{tabular}')
(R/'report/comparison.tex').write_text('\n'.join(table)+'\n')
table=['\\begin{tabular}{crrrr}\\toprule',
       'Modèle, $p=24$ & $b$ (approx.) & $w$ (approx.) & $2/(bw)$ (approx.) & $2^C$\\\\\\midrule']
for name in ('S','A'):
    c=lookup[name,24]['capacity']
    vals=[float(F(c[k])) for k in ('rounded_branch_min_lower','internal_mass_coefficient_lower','required_root_factor')]
    # Show diagnostic decimal approximations only, not rounded certification bounds.
    table.append(name+' & '+' & '.join(f'{v:.6f}' for v in vals)+' & 8\\\\')
table.append('\\bottomrule\\end{tabular}')
(R/'report/capacity_values.tex').write_text('\n'.join(table)+'\n')
peak=max(x['max_TV_over_epsilon'] for x in data['STREAMS']['rows'])
numbers={'GraphStates':data['GRAPHS']['total_states'],
         'GraphEdges':data['GRAPHS']['total_transitions'],
         'UnusedIDs':data['GRAPHS']['total_unused_ids'],
         'StreamSteps':sum(x['steps'] for x in data['STREAMS']['rows']),
         'StreamPeak':f'{peak:.9f}',
         'SplitChecks':data['CONTROLS']['exact_split_checks'],
         'PowerPrecision':data['CONTROLS']['tiny_base_power_precision']}
(R/'report/numbers.tex').write_text(''.join('\\newcommand{\\'+k+'}{'+str(v)+'}\n' for k,v in numbers.items()))
print(numbers)
