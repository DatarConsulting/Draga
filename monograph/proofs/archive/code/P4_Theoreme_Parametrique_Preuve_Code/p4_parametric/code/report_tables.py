"""Build report tables directly from completed results; never invent measurements."""
from pathlib import Path
from fractions import Fraction as F
import json
from certificate import ceil,ceil_log2
ROOT=Path(__file__).resolve().parents[1]
def read(name):
    obj=json.loads((ROOT/'results'/name).read_text())
    assert obj['passed'],name
    return obj
def table(file,columns,header,rows):
    text='\\begin{center}\\small\\begin{tabular}{'+columns+'}\\toprule\n'
    text+=' & '.join(header)+' \\\\\\midrule\n'
    text+='\n'.join(' & '.join(map(str,row))+' \\\\' for row in rows)
    text+='\n\\bottomrule\\end{tabular}\\end{center}\n'
    (ROOT/'report'/file).write_text(text)
def main():
    graphs=read('GRAPHS.json')['configs']
    extra=read('GENERAL_INSTRUMENT.json')['config']
    neighbors=read('NEIGHBORHOOD.json')['configs']
    allrows=graphs+[extra]+neighbors+[read('ROOT_LEAF.json')['config']]
    products=read('PRODUCTS.json')
    refs=read('ORACLE_REFERENCE.json')
    macros={'TotalGraphs':len(allrows),'TotalStates':sum(x['states'] for x in allrows),
            'TotalEdges':sum(x['transitions'] for x in allrows),'TotalOutputs':sum(x['outputs'] for x in allrows),
            'MaxGraphTV':format(max(x['max_TV_over_epsilon_decimal'] for x in allrows),'.6f'),
            'PointChecks':products['exact_point_checks'],'RelativeChecks':products['relative_error_component_checks'],
            'OracleChecks':refs['checks']}
    (ROOT/'report/numbers.tex').write_text('\n'.join('\\newcommand{\\'+k+'}{'+str(v)+'}' for k,v in macros.items())+'\n')
    labels={'central':'Central','asym_A':'A','asym_B':'B','asym_C':'C','asym_D':'D',
            'instrument_E':'E','box_minus':'A, perturbation $-$','box_plus':'A, perturbation $+$',
            'root_leaf':'Racine feuille'}
    rows=[]
    for name in labels:
        group=[r for r in allrows if r['name']==name]
        counts=[next((str(r['states']) for r in group if r['p']==p),'--') for p in (4,8,12)]
        rows.append([labels[name]]+counts+[format(max(r['max_TV_over_epsilon_decimal'] for r in group),'.6f')])
    table('graphs.tex','lrrrr',['Modèle','$N_4$','$N_8$','$N_{12}$','TV$/\\eps$ max.'],rows)
    constants=[];machine=[]
    for name in ['central','asym_A','asym_B','asym_C','asym_D','instrument_E']:
        c=json.loads((ROOT/'results'/('certificate_'+name+'.json')).read_text())
        get=lambda k:F(c[k])
        minus=ceil_log2(4*get('distortion')**4/(get('observation_min')*get('separation')))
        plus=c['offset_bits']+2+ceil(1/(64*c['h']*(1-get('grho'))))
        degree=ceil(1/(1-get('tau')))+6
        constants.append([labels[name],format(float(get('rho')),'.4f'),format(float(get('tau')),'.4f'),minus,plus,degree])
        machine.append({'name':name,'C_minus':minus,'C_plus':plus,'coarse_polynomial_degree':degree})
    table('constants.tex','lrrrrr',['Certificat','$\\rho$','$\\tau$','$C_-$','$C_+$','Degré majorant'],constants)
    streams=read('STREAMS.json')['rows']
    table('streams.tex','lrrrrr',['Modèle','$p$','$B$','TV$/\\eps$ max.','Pic (Kio)','Pas (ms)'],
          [[labels[r['name']],r['p'],r['bits'],format(r['max_TV_over_epsilon_decimal'],'.6f'),
            format(r['model_plus_online_peak_bytes']/1024,'.2f'),format(r['median_update_and_output_seconds']*1000,'.2f')] for r in streams])
    result={'passed':True,'aggregates':macros,'theorem_constants':machine,
            'scope':'Constants computed from certificates; measured totals aggregated from completed runs.'}
    (ROOT/'results/SUMMARY.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))
if __name__=='__main__':main()
