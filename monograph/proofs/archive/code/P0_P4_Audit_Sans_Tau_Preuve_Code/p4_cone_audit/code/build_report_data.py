"""Build report tables from completed JSON, with an explicit preview mode."""
from pathlib import Path
import json,sys
ROOT=Path(__file__).resolve().parents[1]
preview='--preview' in sys.argv
def read(name):return json.loads((ROOT/'results'/name).read_text())
if not preview:
    assert read('RUN_COMPLETE.json')['passed']
    assert read('GRAPHS.json')['passed']
    assert read('CONTROLS.json')['passed']
    assert read('BOUNDARY.json')['passed'] and read('ROOT_LEAF.json')['passed']
graphs=read('GRAPHS.json');exact=read('EXACT_ORACLE_AUDIT.json')
old=json.loads((ROOT/'reference/old_project/results/GRAPHS.json').read_text())
numeric=read('NUMERIC_GLOBAL_REFERENCE.json') if (ROOT/'results/NUMERIC_GLOBAL_REFERENCE.json').exists() else {'checks':0}
commands={'NewGraphs':len(graphs['configs']),'NewStates':graphs['total_states'],
 'NewEdges':graphs['total_transitions'],'NewOutputs':graphs['total_outputs'],
 'MaxTV':f"{graphs['max_TV_over_epsilon']:.6f}",
 'WordChecks':exact['word_D_and_u_checks'],'JetChecks':exact['exact_Taylor_coefficients'],
 'PowerChecks':exact['exact_power_checks'],'NumericChecks':numeric['checks'],
 'OldGraphs':len(old['configs']),'OldStates':old['total_states'],'OldEdges':old['total_transitions']}
(ROOT/'report/numbers.tex').write_text('\n'.join('\\newcommand{\\'+k+'}{'+str(v)+'}' for k,v in commands.items())+'\n')
labels={'beyond_symmetric':'S nouveau','beyond_asymmetric':'A nouveau','asym_A':'Régression'}
table=['\\begin{center}\\small\\begin{tabular}{lrrrrrr}\\toprule',
 'Modèle & $p$ & États & Transitions & Bits ID & Horizon & TV/$\\eps$\\\\\\midrule']
for x in graphs['configs']:
 table.append(f"{labels[x['name']]} & {x['p']} & {x['states']} & {x['transitions']} & {x['bits']} & {x['oracle']['horizon']} & {x['max_TV_over_epsilon_decimal']:.5f}\\\\")
table+=['\\bottomrule\\end{tabular}\\end{center}']
(ROOT/'report/graphs.tex').write_text('\n'.join(table)+'\n')
if (ROOT/'results/STREAMS.json').exists():
 rows=read('STREAMS.json')['rows'];tab=['\\begin{center}\\small\\begin{tabular}{lrrrrrr}\\toprule',
 'Modèle & $p$ & Bits ID & Horizon & Degré Taylor & TV/$\\eps$ & ms/pas\\\\\\midrule']
 for x in rows:
  tab.append(f"{labels[x['name']]} & {x['p']} & {x['B']} & {x['oracle_horizon']} & {x['Taylor_degree']} & {x['max_TV_over_epsilon']:.5f} & {x['milliseconds_per_step']:.2f}\\\\")
 tab+=['\\bottomrule\\end{tabular}\\end{center}']
else:tab=['Aperçu de mise en page : les flux sont encore en cours.']
(ROOT/'report/streams.tex').write_text('\n'.join(tab)+'\n')
if not preview:
 status={'date':'2026-11-02','P0':{'scope_frozen':True,'targeted_literature_audit_completed':True,
   'scientific_priority_established':False,'full_text_gaps':['arXiv:2609.38524','arXiv:2610.13881'],
   'claim':'Executable uniform worst-case optimality for a fixed certified positive rank-two class; proposed combination, not established novelty.'},
   'P4_audit':{'internal_replay_completed':True,'external_independent_review':False,'formal_verification':False,
   'correction':'Preparation time must count arbitrary rational input length kappa, in addition to precision p.'},
   'extension':{'old_tau_condition_removed':True,'primal_contraction_required':True,'dual_contraction_required':True,
   'separation_required':True,'all_rank_two_models':False,'proof_status':'PROPOSED_WITH_EXECUTED_TESTS',
   'global_residual_certificate':'Taylor intervals plus proven grid fallback',
   'finite_precision':'Integer square roots and directed interval arithmetic; coefficients and temporaries counted',
   'time_scope':'Polynomial in p for H fixed, with additional polynomial input-kappa preparation cost',
   'memory_scope':'O_H(p)+O(|H|); initial exact source adds O(kappa)',
   'novelty_of_Birkhoff_method':False},'totals':commands,
   'controls':read('CONTROLS.json'),'completed_run':read('RUN_COMPLETE.json')}
 (ROOT/'results/AUDIT_STATUS.json').write_text(json.dumps(status,indent=2)+'\n')
print(json.dumps({'preview':preview,'numbers':commands}))
