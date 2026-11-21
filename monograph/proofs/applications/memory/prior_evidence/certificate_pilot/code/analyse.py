"""Aggregate the frozen 50-run experiment; never train or select by test."""
import base64
import csv
import hashlib
import html
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SEEDS = list(range(101, 106))
ARMS = ['main_only', 'single', 'single_param', 'multi_path', 'multi_cert']
CONTROLS = ['main_only', 'single', 'single_param']
LABELS = {
    'main_only': 'Principale seule',
    'single': 'Résidu : une transition',
    'single_param': 'Une transition, poids ajustés',
    'multi_path': 'Multilettre, contraction individuelle',
    'multi_cert': 'Multilettre, stabilité de comptage',
}
SHORT = ['Principale\nseule', 'Résidu\nunique', 'Unique, poids\najustés', 'Multilettre\nindividuelle', 'Multilettre\ncomptage']
COLORS = ['#7d8899', '#819fba', '#b4a580', '#d48a44', '#087f8c']


def read(path):
    return json.loads((ROOT / path).read_text())


def stat(values):
    a = np.asarray(values, dtype=float)
    return {'mean': float(a.mean()), 'sample_sd': float(a.std(ddof=1)),
            'values': a.tolist(), 'n': len(a)}


def fmt(s, digits=3, multiplier=1):
    return f"{s['mean'] * multiplier:.{digits}f} ± {s['sample_sd'] * multiplier:.{digits}f}"


def table(headers, rows):
    return '<div class="tablewrap"><table><thead><tr>' + ''.join(
        '<th>' + html.escape(str(x)) + '</th>' for x in headers
    ) + '</tr></thead><tbody>' + ''.join(
        '<tr>' + ''.join('<td>' + html.escape(str(x)) + '</td>' for x in row) + '</tr>' for row in rows
    ) + '</tbody></table></div>'


def aggregate():
    manifest = read('results/run_manifest.json')
    for filename, expected in manifest['hashes'].items():
        actual = hashlib.sha256((ROOT / filename).read_bytes()).hexdigest()
        assert actual == expected, (filename, 'changed since launch')
    outcomes = read('results/run_outcomes.json')
    assert len(outcomes['runs']) == 50 and all(x['exit_code'] == 0 for x in outcomes['runs'])
    runs = {}
    rows = []
    for task in ['text', 'order']:
        for arm in ARMS:
            for seed in SEEDS:
                stem = f'{task}_{arm}_seed{seed}'
                r = read(f'results/{stem}.json')
                assert (r['task'], r['arm'], r['seed']) == (task, arm, seed)
                assert r['state_bytes_fp32'] == 256
                metric = 'bpc' if task == 'text' else 'loss'
                best = min(r['trace'], key=lambda z: z['validation'][metric])
                assert best['step'] == r['selected_step']
                assert abs(best['validation'][metric] - r['validation_best']) < 1e-10
                assert (ROOT / 'checkpoints' / f'{stem}.pt').is_file()
                runs[(task, arm, seed)] = r
                row = {'task': task, 'arm': arm, 'seed': seed,
                       'parameters': r['parameters'], 'parameter_bytes_fp32': r['parameter_bytes_fp32'],
                       'state_bytes_fp32': 256, 'selected_step': r['selected_step'],
                       'training_seconds': r['training_seconds'],
                       'cpu_us_per_token_window64': r['diagnostics']['cpu_us_per_token_batch1_window64']}
                if task == 'text':
                    row.update(bpc=r['test']['bpc'], character_perplexity=r['test']['perplexity'])
                else:
                    row.update({f'accuracy_{L}': r['test'][str(L)]['accuracy'] for L in [32, 64, 128]})
                    row.update({f'pair_both_correct_{L}': r['test'][f'pairs{L}']['pair_both_correct'] for L in [32, 96]})
                    row.update({f'pair_predictions_different_{L}': r['test'][f'pairs{L}']['pair_predictions_different'] for L in [32, 96]})
                rows.append(row)
    summary = {'seeds': SEEDS, 'number_of_runs': 50, 'state_bytes_fp32': 256,
               'aggregation': 'Mean and sample SD across five seeds; no significance test.',
               'arms': {}, 'paired_comparisons': {}, 'run_manifest': manifest,
               'wall_seconds': outcomes['wall_seconds']}
    for task in ['text', 'order']:
        summary['arms'][task] = {}
        for arm in ARMS:
            rr = [runs[(task, arm, seed)] for seed in SEEDS]
            a = {'parameters': rr[0]['parameters'], 'state_bytes_fp32': 256,
                 'parameter_bytes_fp32': rr[0]['parameter_bytes_fp32'],
                 'transition_modes': rr[0]['transition_modes'],
                 'selected_steps': [r['selected_step'] for r in rr],
                 'training_seconds': stat([r['training_seconds'] for r in rr]),
                 'cpu_us_per_token_window64': stat([r['diagnostics']['cpu_us_per_token_batch1_window64'] for r in rr])}
            if task == 'text':
                a.update(bpc=stat([r['test']['bpc'] for r in rr]),
                         character_perplexity=stat([r['test']['perplexity'] for r in rr]),
                         test_predictions=rr[0]['test']['tokens'])
            else:
                a.update({f'accuracy_{L}': stat([r['test'][str(L)]['accuracy'] for r in rr]) for L in [32, 64, 128]})
                a.update({f'pair_both_correct_{L}': stat([r['test'][f'pairs{L}']['pair_both_correct'] for r in rr]) for L in [32, 96]})
                a.update({f'pair_predictions_different_{L}': stat([r['test'][f'pairs{L}']['pair_predictions_different'] for r in rr]) for L in [32, 96]})
            if arm != 'main_only':
                rho = [x for r in rr for x in r['diagnostics']['counting_cp_spectral_radius_per_head']]
                a.update(cp_rho_min=min(rho), cp_rho_max=max(rho),
                         numerical_cp_stable_runs=sum(r['diagnostics']['counting_stability_sufficient_condition'] for r in rr),
                         commutator_frobenius=stat([r['diagnostics']['mean_commutator_frobenius'] for r in rr]))
            summary['arms'][task][arm] = a
    for arm in CONTROLS + ['multi_path']:
        text_gains = [runs[('text', arm, s)]['test']['bpc'] - runs[('text', 'multi_cert', s)]['test']['bpc'] for s in SEEDS]
        order_gains = [runs[('order', 'multi_cert', s)]['test']['pairs96']['pair_both_correct'] - runs[('order', arm, s)]['test']['pairs96']['pair_both_correct'] for s in SEEDS]
        summary['paired_comparisons'][arm] = {
            'text_bpc_control_minus_cert': stat(text_gains),
            'order_pair96_cert_minus_control': stat(order_gains),
            'text_direction_positive_seeds': sum(x > 0 for x in text_gains),
            'order_direction_positive_seeds': sum(x > 0 for x in order_gains),
            'text_criterion_pass': float(np.mean(text_gains)) >= .02 and all(x > 0 for x in text_gains),
            'order_criterion_pass': float(np.mean(order_gains)) >= .10 and all(x > 0 for x in order_gains),
        }
    summary['text_success'] = all(summary['paired_comparisons'][a]['text_criterion_pass'] for a in CONTROLS)
    summary['order_success'] = all(summary['paired_comparisons'][a]['order_criterion_pass'] for a in CONTROLS)
    summary['joint_empirical_success'] = summary['text_success'] and summary['order_success']
    cert = read('certificates/summary.json')
    replay = read('results/compression_replay.json')
    verification = read('results/certificate_verification.json')
    assert verification['all_verified'] and verification['tampered_bound_rejected']
    audit = read('results/audit.json')
    assert audit['all_checks_passed']
    assert audit['checks']['validation_only_checkpoint_selection']['present_runs'] == 50
    assert len(cert['runs']) == len(replay['runs']) == 10
    assert {(r['task'], r['seed']) for r in cert['runs']} == {(t, s) for t in ['text','order'] for s in SEEDS}
    assert {(r['task'], r['seed']) for r in replay['runs']} == {(t, s) for t in ['text','order'] for s in SEEDS}
    for r in cert['runs']:
        full = read('certificates/' + r['certificate_file'])
        assert full['checkpoint_sha256'] == hashlib.sha256((ROOT / 'checkpoints' / r['checkpoint']).read_bytes()).hexdigest()
        assert len(full['heads']) == 4
        assert full['all_stability_checks_passed'] and full['all_error_bounds_exactly_rounded_up']
    summary['compatibility_exact_checks_passed'] = cert['all_passed']
    summary['joint_success'] = summary['joint_empirical_success'] and summary['compatibility_exact_checks_passed']
    summary['certificates'] = cert
    summary['certificate_verification'] = verification
    summary['compression'] = {'source_state_bytes':256, 'projected_state_bytes':192,
                              'state_reduction_fraction':.25, 'by_task':{}, 'runs':replay['runs'],
                              'scope':replay['scope']}
    for task in ['text', 'order']:
        rr = sorted([r for r in replay['runs'] if r['task'] == task], key=lambda r:r['seed'])
        c = {k:stat([r['numeric_replay'][k] for r in rr]) for k in [
            'tv_mean','tv_max','argmax_flip_fraction','max_logit_error','computed_TV_upper_mean',
            'fraction_computed_TV_upper_below_one','max_error_to_computed_bound_ratio']}
        c['maximum_tv_over_all_positions_and_seeds'] = max(r['numeric_replay']['tv_max'] for r in rr)
        c['numerical_bound_consistency_violations'] = sum(r['numeric_replay']['numeric_bound_violations_with_1e_minus5_tolerance'] for r in rr)
        if task == 'text':
            c.update(original_bpc=stat([r['original']['bpc'] for r in rr]),
                     projected_bpc=stat([r['projected']['bpc'] for r in rr]),
                     projected_minus_original_bpc=stat([r['projected']['bpc']-r['original']['bpc'] for r in rr]))
        else:
            c.update({f'{kind}_accuracy_{L}':stat([r[kind][str(L)]['accuracy'] for r in rr])
                      for kind in ['original','projected'] for L in [32,64,128]})
            c.update({f'{kind}_pair_both_correct_{L}':stat([r[kind][f'pairs{L}']['pair_both_correct'] for r in rr])
                      for kind in ['original','projected'] for L in [32,96]})
        summary['compression']['by_task'][task] = c
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with (ROOT / 'results/summary.csv').open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields); w.writeheader(); w.writerows(rows)
    return summary, runs, rows


def plot(summary):
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8), constrained_layout=True)
    fig.suptitle('256 octets d’état récurrent · cinq initialisations · moyenne ± écart-type', fontsize=14)
    x = np.arange(len(ARMS))
    for i, arm in enumerate(ARMS):
        s = summary['arms']['text'][arm]['bpc']
        axes[0].errorbar(i, s['mean'], yerr=s['sample_sd'], fmt='o', color=COLORS[i], capsize=4)
        axes[0].scatter(i + np.linspace(-.13, .13, 5), s['values'], s=18, color=COLORS[i], alpha=.6)
        p = summary['arms']['order'][arm]['pair_both_correct_96']
        axes[1].errorbar(i, p['mean'] * 100, yerr=p['sample_sd'] * 100, fmt='o', color=COLORS[i], capsize=4)
        axes[1].scatter(i + np.linspace(-.13, .13, 5), np.array(p['values']) * 100, s=18, color=COLORS[i], alpha=.6)
        vals = [summary['arms']['order'][arm][f'accuracy_{L}'] for L in [32, 64, 128]]
        axes[2].errorbar([32, 64, 128], [s['mean'] * 100 for s in vals],
                         yerr=[s['sample_sd'] * 100 for s in vals], marker='o', color=COLORS[i],
                         capsize=3, label=LABELS[arm], linewidth=1.5)
    axes[0].set(title='Texte neuf : Sophie', ylabel='Bits par caractère ↓', xticks=x, xticklabels=SHORT)
    axes[1].set(title='Paires ab / ba, longueur 96', ylabel='Deux réponses correctes (%) ↑', xticks=x, xticklabels=SHORT)
    axes[2].set(title='Suivi de permutations S₃', ylabel='État final correct (%) ↑', xlabel='Longueur (entraînement : 32)', xticks=[32, 64, 128])
    axes[2].axhline(100 / 6, color='#999999', linestyle=':', label='Hasard uniforme : 1/6')
    axes[2].legend(fontsize=7.5, loc='best')
    for ax in axes:
        ax.spines[['right', 'top']].set_visible(False)
        ax.grid(axis='y', alpha=.16)
    axes[0].tick_params(axis='x', labelsize=8)
    axes[1].tick_params(axis='x', labelsize=8)
    pair_limits = axes[1].get_ylim()
    axes[1].set_ylim(bottom=min(-.5, pair_limits[0]), top=max(5., pair_limits[1]))
    axes[2].set_ylim(bottom=0)
    fig.savefig(ROOT / 'results/comparison.png', dpi=180)
    plt.close(fig)


def compression_sections(summary):
    cert = summary['certificates']
    rows = []
    for task in ['text','order']:
        rr = [r for r in cert['runs'] if r['task']==task]
        eta = [x for r in rr for x in r['eta_per_head']]
        rows.append(['Texte' if task=='text' else 'Ordre',
                     f"{sum(r['all_stability_checks_passed'] for r in rr)}/5",
                     '20 sources + 20 projections',
                     f'{min(eta):.5g} – {max(eta):.5g}',
                     f"{min(r['uniform_TV_bounds_float']['64'] for r in rr)*100:.2f} – {max(r['uniform_TV_bounds_float']['64'] for r in rr)*100:.2f}"])
    comp = summary['compression']['by_task']
    replay = summary['compression']['runs']
    individual = []
    for task in ['text','order']:
        for r in sorted([r for r in replay if r['task']==task], key=lambda r:r['seed']):
            n = r['numeric_replay']
            metric = (f"{r['original']['bpc']:.4f} → {r['projected']['bpc']:.4f}" if task=='text'
                      else f"{r['original']['32']['accuracy']*100:.2f} → {r['projected']['32']['accuracy']*100:.2f}")
            individual.append(['Texte (BPC)' if task=='text' else 'Ordre (L32 %)',r['seed'],metric,
                               f"{n['tv_mean']*100:.3f}",f"{n['tv_max']*100:.3f}",
                               f"{n['computed_TV_upper_mean']*100:.3f}",
                               f"{n['fraction_computed_TV_upper_below_one']*100:.2f}"])
    fully_vacuous_text = all(r['numeric_replay']['fraction_computed_TV_upper_below_one']==0 for r in replay if r['task']=='text')
    utility = ('<strong>Sur le texte, toutes les bornes de TV évaluées valent 100 %.</strong> Elles sont valides mais ne contraignent pas utilement la différence de prédiction. '
               if fully_vacuous_text else '')
    utility += (f"Sur l’ordre, la borne moyenne calculée vaut {comp['order']['computed_TV_upper_mean']['mean']*100:.2f} % et seulement {comp['order']['fraction_computed_TV_upper_below_one']['mean']*100:.2f} % des positions ont une borne strictement inférieure à 100 %. "
                'La validité mathématique ne suffit donc pas à établir une garantie de compression utile.')
    return f'''<p><strong>Compatibilité exacte : {'vérifiée' if summary['compatibility_exact_checks_passed'] else 'non vérifiée'}.</strong>
Les 40 têtes sources et leurs 40 projections ont chacune une vérification rationnelle de stabilité ; les dix modèles sélectionnés sont couverts. Les fichiers de preuve conservent les coefficients effectifs, pivots LDL, candidats de Stein, résidus exacts et bornes rationnelles arrondies vers le haut. Le vérificateur indépendant <code>code/verify_certificate.py</code> utilise des fractions exactes, retrouve les coefficients des checkpoints, recalcule les inégalités pour les dix modèles (10/10) et rejette un certificat volontairement altéré. Son compte rendu est <code>results/certificate_verification.json</code>.</p>
{table(['Tâche','Modèles vérifiés','Têtes/projections vérifiées','η par tête, min–max','Borne uniforme TV à t64 (%)'],rows)}
<p>Pour chaque réduction, η majore la norme opérateur du Hankel d’erreur via sa norme Hilbert–Schmidt. Les équations de Stein sont résolues numériquement, puis leurs résidus sont recalculés exactement sur des rationnels. La correction rigoureuse des résidus rend la borne indépendante de la précision de la résolution numérique. Le théorème de transport donne <code>||Δlogitsₜ||∞ ≤ Σₕ ||vₕ||∞ ηₕ ||uₕ,1:t||₂</code>, puis <code>TV ≤ min(1, δₜ/2)</code>.</p>
<h2>Compression prévue : 256 → 192 octets d’état</h2>
<p>Chaque tête passe de8 à4 coordonnées en gardant les quatre premières, sans réentraînement et sans choix de coordonnées fondé sur le test. La mémoire récurrente passe de64 à48 scalaires, soit <strong>25 % d’octets d’état en moins</strong>. C’est une approximation réalisable avec une borne d’erreur, <strong>pas une solution optimale du problème AAK</strong>. Le nombre total d’octets du modèle exécuté n’est pas mesuré et les checkpoints restent ceux des modèles complets.</p>
{table(['Mesure','256 octets','192 octets'],[
 ['Texte : BPC ↓',fmt(comp['text']['original_bpc']),fmt(comp['text']['projected_bpc'])],
 ['Ordre L32 : exactitude (%) ↑',fmt(comp['order']['original_accuracy_32'],1,100),fmt(comp['order']['projected_accuracy_32'],1,100)],
 ['Ordre L128 : exactitude (%) ↑',fmt(comp['order']['original_accuracy_128'],1,100),fmt(comp['order']['projected_accuracy_128'],1,100)],
 ['Ordre paires L96 : réussite (%) ↑',fmt(comp['order']['original_pair_both_correct_96'],1,100),fmt(comp['order']['projected_pair_both_correct_96'],1,100)]])}
<p>Variation textuelle projetée moins originale : <strong>{fmt(comp['text']['projected_minus_original_bpc'],4)} BPC</strong>. La TV observée compare les distributions de sortie des deux réalisations sur les mêmes préfixes ; elle ne mesure pas l’exactitude par rapport à la cible. Pour l’ordre, le replay de TV porte sur tous les préfixes de1 024 séquences L128 ; les métriques d’exactitude portent sur l’état final.</p>
{table(['Tâche / métrique','Graine','Original → projeté','TV moyenne (%)','TV maximale (%)','Borne TV moyenne (%)','Bornes <100 % (%)'],individual)}
<p class="note">{utility}</p>
<p><strong>Portée exacte :</strong> les preuves concernent l’arithmétique réelle avec les coefficients dyadiques sauvegardés, des états initiaux nuls et les mêmes routes/injections. Les statistiques de replay et les évaluations des bornes dans le tableau sont numériques. L’absence de violation numérique avec une tolérance de10⁻⁵ ne constitue pas une preuve des arrondis FP32. Les continuations générées librement, où les tokens divergent, ne sont pas couvertes.</p>'''


def report(summary, runs, rows):
    yes = lambda v: 'atteint' if v else 'non atteint'
    img = base64.b64encode((ROOT / 'results/comparison.png').read_bytes()).decode()
    text_rows = [[LABELS[a], fmt(summary['arms']['text'][a]['bpc']),
                  fmt(summary['arms']['text'][a]['character_perplexity'], 2),
                  summary['arms']['text'][a]['parameters'], 256,
                  fmt(summary['arms']['text'][a]['cpu_us_per_token_window64'], 1)] for a in ARMS]
    order_rows = [[LABELS[a]] + [fmt(summary['arms']['order'][a][f'accuracy_{L}'], 1, 100) for L in [32, 64, 128]]
                  + [fmt(summary['arms']['order'][a][f'pair_both_correct_{L}'], 1, 100) for L in [32, 96]] for a in ARMS]
    criterion_rows = []
    for a in CONTROLS:
        c = summary['paired_comparisons'][a]
        criterion_rows.append([LABELS[a], fmt(c['text_bpc_control_minus_cert'], 4),
                               f"{c['text_direction_positive_seeds']}/5", yes(c['text_criterion_pass']),
                               fmt(c['order_pair96_cert_minus_control'], 1, 100),
                               f"{c['order_direction_positive_seeds']}/5", yes(c['order_criterion_pass'])])
    stability_rows = []
    for task in ['text', 'order']:
        for a in ARMS[1:]:
            s = summary['arms'][task][a]
            stability_rows.append(['Texte' if task == 'text' else 'Ordre', LABELS[a],
                                   f"{s['cp_rho_min']:.6f} – {s['cp_rho_max']:.6f}",
                                   f"{s['numerical_cp_stable_runs']}/5", fmt(s['commutator_frobenius'], 5)])
    t_individual = [[LABELS[r['arm']], r['seed'], f"{r['bpc']:.6f}", f"{r['character_perplexity']:.6f}",
                    r['parameters'], r['parameter_bytes_fp32'], r['selected_step']] for r in rows if r['task'] == 'text']
    o_individual = [[LABELS[r['arm']], r['seed']] + [f"{r[f'accuracy_{L}'] * 100:.4f}" for L in [32, 64, 128]]
                    + [f"{r[f'pair_both_correct_{L}'] * 100:.4f}" for L in [32, 96]]
                    + [r['parameters'], r['selected_step']] for r in rows if r['task'] == 'order']
    pair_values = summary['arms']['order']['multi_cert']['pair_both_correct_96']['values']
    c = summary['paired_comparisons']['multi_path']
    constraint_text = fmt(c['text_bpc_control_minus_cert'], 4)
    constraint_order = fmt(c['order_pair96_cert_minus_control'], 1, 100)
    audit = read('results/audit.json')
    protocol = html.escape((ROOT / 'PROTOCOL.md').read_text())
    doc = f'''<!doctype html><html lang="fr"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Mémoire multilettre : gain sous contrainte de stabilité</title><style>
:root{{--ink:#1e303d;--muted:#526775;--accent:#087f8c;--line:#d9e2e6}}*{{box-sizing:border-box}}
body{{margin:0;background:#f2f5f6;color:var(--ink);font:16px/1.55 system-ui,-apple-system,Segoe UI,sans-serif}}main{{max-width:1200px;margin:28px auto;padding:36px;background:white;border:1px solid var(--line);border-radius:10px}}
h1{{font-size:2rem;line-height:1.18;margin:.4rem 0 1rem}}h2{{font-size:1.25rem;margin:2rem 0 .6rem}}h3{{font-size:1rem}}p{{max-width:105ch}}
.eyebrow{{font-size:.75rem;letter-spacing:.08em;text-transform:uppercase;color:var(--accent);font-weight:700}}.lede{{font-size:1.08rem}}.meta{{display:flex;gap:10px;flex-wrap:wrap}}.pill{{background:#edf2f4;padding:6px 10px;border-radius:5px;font-size:.85rem}}
.verdict{{border-left:5px solid var(--accent);padding:16px 20px;background:#edf7f7;margin:22px 0}}.note{{padding:14px 18px;background:#fff8eb;border:1px solid #e6d3aa}}.small,figcaption{{font-size:.85rem;color:var(--muted)}}
.tablewrap{{overflow-x:auto;margin:14px 0}}table{{border-collapse:collapse;width:100%;font-size:.85rem}}th,td{{text-align:left;padding:9px 10px;border-bottom:1px solid var(--line);vertical-align:top}}th{{background:#eef3f5}}tbody tr:nth-child(even){{background:#f8fafb}}
figure{{margin:22px 0}}img{{width:100%;height:auto}}a{{color:#006b79}}code,pre{{font-family:ui-monospace,monospace}}pre{{white-space:pre-wrap;font-size:.8rem;background:#f4f6f7;border:1px solid var(--line);padding:16px}}li{{margin:.4rem 0}}
@media(max-width:700px){{main{{padding:20px;margin:0;border:0;border-radius:0}}h1{{font-size:1.6rem}}}}@media print{{body{{background:white}}main{{border:0;margin:0;padding:0}}h2{{break-after:avoid}}tr{{break-inside:avoid}}details{{display:none}}}}
</style></head><body><main>
<div class="eyebrow">Expérience exécutée · 1 novembre 2026 · CPU</div>
<h1>Mémoire multilettre : le gain subsiste-t-il sous contrainte de stabilité ?</h1>
<p class="lede">Un nouveau corpus français, une lecture scalaire additive et une condition de stabilité vérifiable exactement permettent de tester le passage de l’architecture empirique au composant compatible avec le certificat.</p>
<div class="meta"><span class="pill">50 entraînements</span><span class="pill">5 graines : 101–105</span><span class="pill">256 octets d’état FP32</span><span class="pill">5 variantes</span></div>
<div class="verdict"><strong>Critère empirique conjoint : {yes(summary['joint_empirical_success'])}.</strong>
<p>Critère texte : <strong>{yes(summary['text_success'])}</strong>. Critère ordre : <strong>{yes(summary['order_success'])}</strong>. Chaque critère exigeait un gain face aux trois témoins, selon les seuils fixés avant l’entraînement. La compatibilité mathématique et l’utilité des bornes de compression sont évaluées séparément ci-dessous.</p>
<p>Multilettre sous stabilité de comptage : <strong>{fmt(summary['arms']['text']['multi_cert']['bpc'])} BPC</strong> et <strong>{fmt(summary['arms']['order']['multi_cert']['pair_both_correct_96'], 1, 100)} %</strong> de paires L96 réussies. Résultats individuels sur ces paires : {', '.join(f'{x*100:.2f} %' for x in pair_values)}.</p></div>
<p>Sur ce corpus, la mémoire principale seule obtient {summary['arms']['text']['main_only']['bpc']['mean']:.3f} BPC, contre {summary['arms']['text']['multi_cert']['bpc']['mean']:.3f} pour le multilettre sous stabilité de comptage. Le témoin multilettre à contraction individuelle obtient lui aussi une moyenne moins bonne que la mémoire principale ({summary['arms']['text']['multi_path']['bpc']['mean']:.3f} BPC). L’absence de gain textuel ne peut donc pas être imputée uniquement à la nouvelle contrainte de stabilité.</p>
<p class="note">La contrainte Σₐ RₐᵀRₐ ≺ 0,99 I est une condition suffisante, plus restrictive que la stabilité spectrale générale. Une contre-performance de ce modèle ne réfute donc pas toutes les architectures compatibles avec le théorème. Le test ne démontre pas non plus une nouveauté architecturale ou une supériorité sur les LLM.</p>
<figure><img alt="Résultats par variante sur le texte, les paires sensibles à l’ordre et les longueurs de permutations" src="data:image/png;base64,{img}"><figcaption>Points : moyennes et graines individuelles. Barres d’erreur : écart-type échantillonnal, sans interprétation comme intervalle de confiance. Sur les exactitudes, moyenne ± écart-type peut dépasser [0,100] ; les valeurs individuelles restent des proportions. L’axe texte est resserré pour rendre les différences visibles.</figcaption></figure>

<h2>Comparaison et rôle des témoins</h2>
{table(['Variante','État principal','État résiduel','Transitions','Lecture du résidu'],[
 [LABELS['main_only'],'64','0','—','—'],
 [LABELS['single'],'32','4 × 8','1 par tête','Scalaire, additive'],
 [LABELS['single_param'],'32','4 × 8','1 par tête','Scalaire, additive ; MLP principal élargi'],
 [LABELS['multi_path'],'32','4 × 8','4 sur texte / 3 sur ordre','Scalaire, additive'],
 [LABELS['multi_cert'],'32','4 × 8','4 sur texte / 3 sur ordre','Scalaire, additive']])}
<p>Toutes les variantes conservent <strong>64 coordonnées FP32, soit 256 octets par exemple</strong>. Cela exclut les poids, les activations, les buffers et le framework. Le code traite des segments complets et conserve les états intermédiaires ; il ne mesure pas une RAM totale de 256 octets. Le quatrième et le cinquième modèles ont les mêmes paramètres bruts initiaux et la même lecture ; seule la normalisation des transitions change.</p>
<p>Le témoin à poids ajustés est le plus proche parmi les largeurs entières du MLP : il compte 40 paramètres de moins que le multilettre sur texte (environ 0,40 %) et5 de moins sur l’ordre (environ 0,11 %). Le nombre de paramètres est donc presque égal, sans égalité exacte.</p>
<p>Chaque tête suit <code>rₜ = Rₐ rₜ₋₁ + b uₜ</code> et fournit <code>v(cᵀrₜ)</code> aux logits. Le MLP lit uniquement l’état principal. Le résidu ne modifie ni le routeur ni la mémoire principale. Sur texte, le symbole de transition est appris à partir du caractère courant, avec argmax en avant et gradient straight-through. Sur S₃, il est donné par l’opération.</p>

<h2>Prédiction du texte neuf</h2>
{table(['Variante','BPC ↓','Perplexité caractère ↓','Paramètres','État (octets)','CPU µs/caractère*'],text_rows)}
<p class="small">Moyenne ± écart-type sur cinq graines. La perplexité est moyennée après calcul par graine. *Temps indicatif : lot1, fenêtre64, médiane de sept passes ; les autres entraînements concurrents peuvent affecter cette mesure. Aucun gain de débit de décodage optimisé n’est revendiqué.</p>
<p>Corpus : <em>Les malheurs de Sophie</em>, comtesse de Ségur, <a href="https://www.gutenberg.org/ebooks/15058">Project Gutenberg n°15058</a>. Les 180 000 premiers caractères du corps normalisé en minuscules sont répartis en 140 000 / 20 000 / 20 000 caractères. Le test est neuf par rapport au pilote Jules Verne. Le vocabulaire utilise seulement l’entraînement ; aucune sous-chaîne de 64 caractères de validation/test n’est retrouvée dans l’entraînement.</p>
<p>Chaque exécution présente 1 800 lots de 32 fenêtres de 64 caractères, soit 3 686 400 cibles avec rééchantillonnage. Le test comporte 19 968 prédictions. La mémoire est remise à zéro au début de chaque fenêtre ; aucune performance en contexte textuel continu long n’en découle.</p>

<h2>Ordre des opérations et extrapolation en longueur</h2>
{table(['Variante','L32 (%)','L64 (%)','L128 (%)','Paires L32 (%)','Paires L96 (%)'],order_rows)}
<p>La tâche suit la permutation complète de trois objets, parmi six classes. Les opérations identité, échange 0↔1 et échange 1↔2 sont composées dans l’ordre. L’entraînement supervise tous les préfixes de longueur32, pour 1 500 lots de64, soit 3 072 000 cibles. Les tests ordinaires ont 1 024 séquences par longueur.</p>
<p>Les 512 paires L96 partagent un préfixe de8 opérations, puis <code>ab</code> versus <code>ba</code>, puis un suffixe commun de86 opérations. Les comptes de symboles sont identiques et les cibles finales diffèrent toujours. Une paire réussit uniquement si <strong>les deux membres sont corrects</strong>. Le hasard uniforme indépendant donne 1/6 par séquence et 1/36 par paire ; ce dernier nombre n’est pas un plancher pour des prédictions corrélées. Répéter la même réponse sur les deux membres donne zéro réussite de paire.</p>

<h2>Vérification des critères fixés à l’avance</h2>
<p>Face à chacun des trois témoins : gain moyen d’au moins <strong>0,02 BPC</strong> et gain moyen d’au moins <strong>10 points</strong> sur les paires L96, avec différence favorable sur les cinq graines pour chaque métrique. Les chiffres ci-dessous sont appariés par graine ; un gain positif favorise le modèle sous stabilité de comptage.</p>
{table(['Témoin','Gain BPC','Graines +','Critère texte','Gain paires (points)','Graines +','Critère ordre'],criterion_rows)}
<p>Comparaison d’isolation <strong>multi_cert contre multi_path</strong> : gain BPC {constraint_text} ; gain paires L96 {constraint_order} points. Cette comparaison porte sur le coût de la contrainte dans la paramétrisation choisie, avec la même lecture scalaire. Elle ne permet pas d’attribuer à un seul changement les différences avec le pilote précédent, qui utilisait un MLP sur toutes les coordonnées et un autre corpus.</p>

<h2>Stabilité : diagnostic numérique et preuve exacte</h2>
{table(['Tâche','Variante','ρ(Σ R⊗R) min–max','Exécutions ρ < 1','Commutateur moyen'],stability_rows)}
<p>Ces rayons spectraux sont des diagnostics numériques sur les coefficients effectifs. La vérification exacte séparée utilise les coefficients FP32 stockés comme nombres rationnels dyadiques et vérifie les pivots LDL de <code>(99/100)I − Σ RᵀR</code>. La normalisation numérique vise 0,9801 I, laissant une marge. Une contraction de chaque transition, seule, ne suffit pas à la stabilité sur la mesure de comptage.</p>
<p><strong>Coût théorique sur les suffixes uniformes.</strong> Pour deux résidus recevant le même suffixe de L opérations iid uniformes parmi d modes, la contrainte implique <code>E||Δr_L||₂² ≤ (q/d)ᴸ ||Δr₀||₂²</code>. Ici d=3, q=0,99 : sur le suffixe commun de86 opérations, le facteur en moyenne quadratique est d’environ <strong>1,98 × 10⁻²¹</strong>. Cette contraction très forte de la différence résiduelle aide à interpréter la tâche d’ordre. Elle concerne le résidu seul, pas la mémoire principale ; elle ne prouve ni la fusion exacte des états FP32 ni l’échec de toute architecture stable. Le routeur textuel, non uniforme, ne satisfait pas automatiquement cette hypothèse.</p>
{compression_sections(summary)}

<h2>Toutes les graines et reproductibilité</h2>
<h3>Texte</h3>{table(['Variante','Graine','BPC','Perplexité','Paramètres','Poids FP32 (octets)','Étape choisie'],t_individual)}
<h3>Ordre</h3>{table(['Variante','Graine','L32 (%)','L64 (%)','L128 (%)','Paires32 (%)','Paires96 (%)','Paramètres','Étape choisie'],o_individual)}
<p>AdamW, taux initial0,003 puis0,0006 pour les 20 % dernières mises à jour, régularisation0,01, clipping1. Les exemples sont identiques par tâche/graine entre variantes. Les budgets égalisent mises à jour et exemples, pas les FLOPs. Les checkpoints sont choisis exclusivement par validation toutes les100 étapes, incluant0. Aucun réglage n’a été effectué après inspection du test.</p>
<p>Le code d’analyse exige les 50 résultats, les 50 checkpoints, l’absence de modification des fichiers figés et la cohérence de la sélection par validation. Audit indépendant : <strong>{'tous les contrôles présents passent' if audit['all_checks_passed'] else 'échec — voir audit.json'}</strong>. Environnement : Python {html.escape(runs[('text','main_only',101)]['python'])}, PyTorch {html.escape(runs[('text','main_only',101)]['torch_version'])}, CPU, un thread par entraînement, trois travailleurs. Durée totale du lancement : {summary['wall_seconds']:.1f} s.</p>
<p class="small">Code, source brute et licence, partitions, checkpoints, courbes de validation, preuves rationnelles, replay et résultats individuels sont inclus dans le dossier reproductible. Cinq graines et un seul petit corpus fournissent une expérience contrôlée limitée ; aucun test de significativité n’est présenté.</p>
<details><summary>Protocole figé avant entraînement</summary><pre>{protocol}</pre></details>
</main></body></html>'''
    (ROOT / 'RAPPORT.html').write_text(doc)


if __name__ == '__main__':
    s, runs, rows = aggregate()
    plot(s)
    report(s, runs, rows)
    (ROOT / 'results/summary.json').write_text(json.dumps(s, indent=2, ensure_ascii=False))
    print(json.dumps({k: s[k] for k in ['number_of_runs','text_success','order_success','joint_empirical_success']}, indent=2))
