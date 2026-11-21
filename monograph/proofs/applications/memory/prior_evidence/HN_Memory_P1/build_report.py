#!/usr/bin/env python3
"""Self-contained French report, with a compact theorem and measured scope."""
import base64
import html
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parent;OUT=ROOT/'results'
load=lambda name:json.loads((OUT/(name+'.json')).read_text())
rows=load('results');summary=load('summary');struct=load('structural_only')
checks=load('independent_checks');validation=load('validation')


def fr(x,d=2):return f'{x:.{d}f}'.replace('.',',')
def sci(x):return f'{x:.4g}'.replace('.',',')
def table(headers,body):
    return '<div class="table"><table><thead><tr>'+''.join('<th>'+h+'</th>' for h in headers)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+str(c)+'</td>' for c in line)+'</tr>' for line in body)+'</tbody></table></div>'


plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,
                     'axes.spines.right':False,'figure.facecolor':'white'})
fig,ax=plt.subplots(1,2,figsize=(12,4.8))
nonzero=struct[1:];ep=[float(r['epsilon']) for r in nonzero]
ax[0].loglog(ep,[r['max_state_error_observed'] for r in nonzero],'o-',c='#187b85',label='Maximum observé')
ax[0].loglog(ep,[r['uniform_state_error_bound'] for r in nonzero],'s--',c='#a9572c',label='Borne uniforme prouvée')
ax[0].set_title('Couplages retirés, sans quantification',loc='left',fontweight='bold')
ax[0].set_ylabel('Erreur sur l’état — norme infinie')
focus=[next(r for r in rows if r['epsilon']==s['epsilon'] and r['budget_bits']==16 and r['transition']=='triangular' and r['policy']=='calibrated' and r['case']==0) for s in nonzero]
for key,label,color,marker in [('observed_rms_error','RMS observée','#187b85','o'),
                              ('rms_aposteriori_bound','Majorant a posteriori RMS','#496b9c','^'),
                              ('uniform_rms_output_bound','Borne uniforme','#a9572c','s')]:
    ax[1].loglog(ep,[r[key] for r in focus],marker+'-',c=color,label=label)
ax[1].set_title('Codec 2 octets, allocation calibrée',loc='left',fontweight='bold')
ax[1].set_ylabel('Erreur sur les sorties / majorant RMS')
for a in ax:
    a.set_xlabel('Couplage ε');a.grid(alpha=.2);a.legend(frameon=False,fontsize=8,loc='upper left')
fig.suptitle('Une fuite faible se contrôle ; la quantification impose sa propre erreur',fontweight='bold',fontsize=13,x=.06,ha='left')
fig.tight_layout(rect=[0,0,1,.92])
fig.savefig(OUT/'bounds.png',dpi=180,bbox_inches='tight');fig.savefig(OUT/'bounds.svg',bbox_inches='tight');plt.close(fig)
img=base64.b64encode((OUT/'bounds.png').read_bytes()).decode()
st=table(['ε','Défaut maximal du drapeau','Erreur max d’état observée','Borne uniforme d’état'],[
    [s['epsilon'],sci(2*float(s['epsilon'])),sci(s['max_state_error_observed']),sci(s['uniform_state_error_bound'])] for s in struct])
sr=[r for r in summary if r['transition']=='triangular' and r['budget_bits']==16]
qt=table(['ε','Allocation','Cas améliorés / 24','Réduction médiane MSE','Pire réduction MSE'],[
    [s['epsilon'],'Calibration' if s['policy']=='calibrated' else 'Borne globale',s['better_cases'],fr(s['median_reduction_pct'])+' %',fr(s['worst_reduction_pct'])+' %'] for s in sr])
primary=[r for r in rows if r['epsilon']=='0.000001' and r['budget_bits']==16 and r['transition']=='triangular' and r['case']==0]
bt=table(['Allocation','Bits','RMS observée','Majorant a posteriori RMS','Borne uniforme'],[
    [r['policy'],str(r['bits']),fr(r['observed_rms_error'],4),fr(r['rms_aposteriori_bound'],4),fr(r['uniform_rms_output_bound'],4)] for r in primary])
count=sum(r[s]['scalar_coordinate_checks'] for r in validation for s in ('calibration','test'))
allchecks=table(['Contrôle','Résultat'],[
    ['Témoins rationnels de boîte invariante et contraction','5 couplages certifiés, de 0 à 10⁻²'],
    ['Rejeux indépendants en fractions exactes','12 traces × 96 pas'],
    ['Inégalités exactes par coordonnée et par pas',checks['exact_pointwise_coordinate_checks']],
    ['Inégalités futures exactes, horizons 1, 8 et 32',checks['exact_future_coordinate_checks_horizons_1_8_32']],
    ['Comparaisons numériques d’erreur / majorants',f'{count:,}'.replace(',',' ')+' ; tolérance 10⁻¹⁰'],
    ['Violations dépassant cette tolérance','0'],
    ['Excès numérique signé maximal sur borne a posteriori','1,67 × 10⁻¹⁶'],
    ['Résultats publiés recalculés depuis les preuves par séquence','2 160'],
    ['Contrôles négatifs : oublier couplage / quantification','Échecs détectés'],
    ['Franchissement local de grille à ε = 10⁻⁶ et 10⁻⁹','Un saut de grille vérifié exactement'],
    ['ε = 0,05','Pas de certificat par cette matrice de comparaison ; aucune conclusion d’instabilité']])

report=f'''<!doctype html><html lang="fr"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>P1 — Filtration approchée et erreur future certifiée</title><style>
:root{{--ink:#182f40;--teal:#187b85}}*{{box-sizing:border-box}}body{{margin:0;background:#edf2f4;color:var(--ink);font:16px/1.62 system-ui,sans-serif}}
main{{max-width:1120px;margin:30px auto;padding:48px 58px;background:white;border-radius:10px}}h1{{font-size:36px;line-height:1.17;margin:14px 0 22px}}
h2{{margin-top:38px;padding-top:24px;border-top:1px solid #d9e4e8;font-size:24px}}h3{{font-size:18px}}.eyebrow{{font-size:12px;color:var(--teal);font-weight:700;letter-spacing:.1em}}
.lead{{font-size:21px}}.box{{padding:18px 22px;border-left:4px solid var(--teal);background:#f0f5f7;margin:22px 0}}.warning{{border-left-color:#ad5e30}}
.equation{{white-space:pre-wrap;background:#f3f6f8;padding:16px 18px;font:14px/1.7 ui-monospace,monospace;overflow-wrap:anywhere}}.table{{overflow-x:auto}}
table{{width:100%;border-collapse:collapse;font-size:14px;margin:18px 0 26px}}th{{background:var(--ink);color:white;padding:10px;text-align:left}}td{{padding:10px;border-bottom:1px solid #dce5e8;vertical-align:top}}
tr:nth-child(even) td{{background:#f5f8f9}}a{{color:#126780}}.muted{{color:#536975}}img{{width:100%;height:auto}}code{{background:#f1f5f6}}
@media(max-width:760px){{main{{margin:0;padding:24px 16px}}h1{{font-size:28px}}table{{font-size:12px}}}}
@media print{{main{{margin:0;padding:0;max-width:none}}body{{background:white}}h2{{break-after:avoid}}table,img,.box{{break-inside:avoid}}}}
</style><main><div class="eyebrow">HN MEMORY · P1 v1.0 · 1 NOVEMBRE 2026</div>
<h1>De l’invariance exacte<br>à une erreur future bornée</h1>
<p class="lead">Le couplage 10⁻⁶ détruit les sous-espaces communs exacts de P0. Le même drapeau reste approximativement invariant, et une borne explicite couvre la fuite entre étages ainsi que la quantification récurrente.</p>
<div class="box"><b>Résultat central.</b> À ε = 10⁻⁶, retirer les six couplages qui cassent le drapeau produit, sans quantification, une erreur maximale d’état mesurée de <b>2,7044 × 10⁻⁶</b>, sous une borne uniforme de <b>5,0003 × 10⁻⁵</b>. Avec quantification, la preuve ajoute le résidu local d’arrondi. Les bornes sont valides dans le cadre déclaré, mais peuvent être trop larges pour choisir efficacement l’allocation.</div>

<h2>1 · Construction et sens précis de la filtration approchée</h2>
<p>On reprend les quatre composantes de P0, sans entraînement sur texte. L’entrée u est bornée par 1 en valeur absolue ; les symboles a appartiennent à {{0,1}}. Les états et sorties de référence sont x<sub>t+1</sub> = A<sub>aₜ</sub>x<sub>t</sub> + bu<sub>t</sub>, y<sub>t</sub> = Cx<sub>t</sub>.</p>
<div class="equation">A₀ = diag(0.98, 0.90, 0.65, 0.25)
A₁ = diag(0.95, 0.75, 0.30, −0.20) + ε (11ᵀ − I)
b = (0.10, 0.25, 0.50, 0.80)ᵀ ; C = diag(g)
Pₖ = diag(1,…,1,0,…,0), avec k ∈ {{1,2,3}}.</div>
<p>Le défaut de fermeture est L<sub>a,k</sub> = (I−P<sub>k</sub>)A<sub>a</sub>P<sub>k</sub>. Ici, L<sub>0,k</sub> = 0 et L<sub>1,k</sub> est une matrice de rang un à coefficients ε dans le bloc inférieur gauche. Par conséquent :</p>
<div class="equation">‖L₁,ₖ‖₂ = ε √[k(4−k)] ; maxₐ,ₖ ‖Lₐ,ₖ‖₂ = 2ε.</div>
<p>Les espaces restent non triviaux comme espaces vectoriels. Pour ε ≠ 0, ils cessent d’être des sous-représentations exactes. On les appelle ici <b>drapeau approché hérité de P0</b> : aucune nouvelle filtration HN canonique du modèle perturbé n’est revendiquée. La norme spectrale du défaut reste la même après la rotation orthogonale Q utilisée dans P0.</p>
<p>On définit T₀ = A₀ et T₁ = triu(A₁), puis Eₐ = Aₐ−Tₐ. Les matrices T préservent exactement le drapeau. Elles retirent six coefficients situés sous la diagonale. <b>Les quatre états mémoire sont conservés</b> ; cette étape simplifie les transitions et n’effectue pas une réduction de dimension. Les couplages au-dessus de la diagonale restent présents.</p>

<h2>2 · Certificat de propagation, valable pour tout mot futur</h2>
<p>Le codec exécute q<sub>t+1</sub> = Q<sub>b</sub>(T<sub>aₜ</sub>q<sub>t</sub> + bu<sub>t</sub>). Notons η<sub>t</sub> l’erreur exacte introduite par cette quantification, et e<sub>t</sub> = x<sub>t</sub>−q<sub>t</sub>. On obtient l’identité :</p>
<div class="equation">e[t+1] = A[aₜ] e[t] + E[aₜ] q[t] − η[t].</div>
<p>Avec v₀ ≥ |e₀|, les valeurs absolues et inégalités étant prises coordonnée par coordonnée, la récurrence suivante donne un majorant :</p>
<div class="equation">v[t+1] = |A[aₜ]| v[t] + |E[aₜ] q[t] − η[t]|
|e[t]| ≤ v[t] ; |y[t] − ŷ[t]| ≤ |C| v[t].</div>
<p><b>Preuve :</b> substituer l’identité d’erreur, appliquer l’inégalité triangulaire et procéder par récurrence. Cette borne se calcule à partir de l’état du codec, de ses matrices et du résidu d’arrondi ; elle ne demande pas l’état natif x. Elle est a posteriori, car elle utilise la transition et le résidu réellement exécutés.</p>
<p>Pour certifier le futur avant de connaître les symboles et entrées, prendre K ≥ |Aₐ|, D ≥ |Eₐ|, une boîte invariante |q| ≤ R et |η| ≤ δ. Poser r = DR + δ. Alors, pour toute continuation admissible et tout horizon h :</p>
<div class="equation">|e[t+h]| ≤ Kʰ v[t] + Σⱼ₌₀ʰ⁻¹ Kʲ r.
|y[t+h] − ŷ[t+h]| ≤ |C| [Kʰ v[t] + Σⱼ₌₀ʰ⁻¹ Kʲ r].</div>
<p>La majoration tient compte des <b>nouvelles quantifications à chaque étape future</b>. Les deux systèmes reçoivent les mêmes entrées et symboles externes. Le certificat ne couvre pas une situation où leurs propres sorties engendrent des entrées futures différentes.</p>
<p>Si K est contractante dans une norme adaptée et e₀ = 0, le majorant uniforme est v̄ = (I−K)⁻¹r. Comme x et q appartiennent à la même boîte, on peut aussi utiliser min(v̄,2R) coordonnée par coordonnée. Dans notre exemple, K est la diagonale des maxima absolus, complétée par ε hors diagonale. Le témoin rationnel R = (I−K)⁻¹b vérifie exactement KR + b = R et R &gt; 0. Il donne une contraction pondérée de facteur maxᵢ (KR)ᵢ/Rᵢ &lt; 1.</p>
<p>La grille de bᵢ bits utilise 2<sup>bᵢ</sup>−1 niveaux dans [−Rᵢ,Rᵢ], dont zéro ; δᵢ = Rᵢ/(2<sup>bᵢ</sup>−2). Un code est inutilisé. La fermeture de la boîte vaut pour la dynamique native et la dynamique triangulaire. La preuve est une majoration classique de perturbations, adaptée à ce codec ; ce n’est pas un nouveau théorème HN.</p>

<h2>3 · Effet des seuls couplages supprimés</h2>
<p>Ce contrôle isole E en désactivant la quantification (δ = 0). Les observations ci-dessous portent sur 128 séquences de test de longueur 512. La borne mathématique couvre tous les temps et toutes les entrées admissibles, à partir de l’état nul.</p>
{st}
<p>Pour les gains de sortie initiaux g = (8 ; 2 ; 0,5 ; 0,125), à ε = 10⁻⁶, le maximum de sortie mesuré vaut 5,4088 × 10⁻⁶, et la borne uniforme vaut 1,00006 × 10⁻⁴. Ces deux chiffres sont <b>sans quantification</b> ; ils ne sont pas les erreurs totales d’un codec de deux octets.</p>

<h2>4 · Mémoire quantifiée : comparaison à budget égal</h2>
<p>On compare transitions complètes et triangulaires, pour cinq couplages, trois budgets (2, 3, 4 octets) et les 24 permutations de gains de sortie de P0. Calibration : 64 séquences × 256 pas ; test : 128 séquences × 512 pas ; 64 pas d’échauffement. Les entrées sont des entiers uniformes dans [−1024,1024] divisés par 1024. Les trajectoires sont partagées entre configurations de gains.</p>
<p>La mémoire est réellement emballée en octets et décodée à chaque pas. Les politiques sont : précision uniforme ; minimisation de la MSE sur calibration ; minimisation de la borne globale analytique, sans données de calibration ni de test. Les tableaux complets contiennent 2 160 enregistrements. Le tableau suivant concerne les transitions triangulaires à <b>2 octets</b>, comparées à l’uniforme sur les mêmes transitions :</p>
{qt}
<p>À ε = 10⁻⁶, la calibration améliore les 24 configurations, avec 91,28 % de réduction médiane de MSE. Le choix par borne globale améliore 23 cas mais en dégrade un. Une garantie valide peut donc rester un mauvais objectif d’optimisation de l’erreur moyenne. À ε = 10⁻², cette difficulté devient plus marquée.</p>
<div class="box warning"><b>Comparaison avec P0 :</b> P1 mesure l’erreur de sortie après chaque stockage quantifié, sur la récurrence déployée entière. P0 mesurait principalement une continuation à un pas sans quantification supplémentaire. Les graines et entrées ont aussi changé. Les pourcentages des deux phases ne doivent pas être comparés comme une progression de performance.</div>
<img src="data:image/png;base64,{img}" alt="Erreur structurelle et bornes avec quantification selon le couplage">
<p class="muted">Panneau droit : gains initiaux, allocation calibrée, 2 octets. Le majorant a posteriori RMS est la racine de la moyenne du carré des bornes par coordonnée ; la borne uniforme couvre chaque instant.</p>
<h3>Exemple quantifié, ε = 10⁻⁶ et gains initiaux</h3>
{bt}
<p>Pour l’allocation calibrée (8,4,2,2), l’erreur RMS mesurée est 0,2626 tandis que la borne uniforme vaut 4,3543. Le majorant a posteriori RMS vaut 1,2808. La validité est établie, mais la précision de la borne demeure une limite pratique. Sur les 24 gains à 2 octets et ε = 10⁻⁶, la borne globale est environ 6,85 à 44,83 fois la RMS observée.</p>

<h2>5 · Un petit couplage peut franchir une cellule de quantification</h2>
<p>Un contrôle local construit un état admissible de la grille et une entrée |u| ≤ 1 pour lesquels retirer un couplage ε = 10⁻⁹ change de côté d’un seuil d’arrondi. Avant quantification, la différence est proportionnelle à ε ; après quantification, elle vaut un pas de grille, environ <b>0,35714</b> sur la deuxième coordonnée. Le calcul est entièrement rationnel.</p>
<p>Ce contrôle est local : l’accessibilité de cet état à partir de zéro avec une entrée commune n’est pas affirmée. Il établit néanmoins pourquoi une borne fondée uniquement sur le défaut de fermeture ne suffit pas pour le codec. Le terme η, ou sa majoration δ, est indispensable. La présence d’une petite fuite n’implique pas l’identité des codes quantifiés.</p>

<h2>6 · Vérifications et coût du certificat</h2>
{allchecks}
<p>Les contrôles numériques massifs réutilisent les mêmes données entre allocations ; leur nombre n’est pas un nombre d’échantillons statistiques indépendants. Les intervalles bootstrap rééchantillonnent des séquences entières, avec 500 tirages appariés et sans correction de multiplicité.</p>
<p>La preuve générale concerne l’arithmétique réelle. Les témoins de contraction, défauts de fermeture et rejeux de 96 pas sont vérifiés en fractions exactes. Les simulations de performance utilisent float64 et une tolérance déclarée de 10⁻¹⁰. Aucun certificat exhaustif de toutes les décisions d’arrondi de toutes les implémentations flottantes n’est revendiqué.</p>
<div class="box"><b>Comptabilité mémoire.</b> La borne globale est précomputée et partagée par codec : elle n’ajoute aucun état par séquence. Suivre le vecteur a posteriori v en ligne demanderait <b>32 octets supplémentaires par séquence</b> en float64. Dans ce test, il est conservé pour audit hors ligne. Les 2/3/4 octets annoncés décrivent uniquement la mémoire persistante du codec, pas ce suivi d’audit ni les buffers temporaires.</div>
<p>Les plages et largeurs de quantification sont des métadonnées statiques ; matrices, buffers de décodage, entrées et sorties sont comptés séparément. Les transitions triangulaires suppriment six coefficients, mais aucune mesure de latence, de pic RAM ou de gain total de stockage de poids n’a été réalisée.</p>

<h2>7 · Ce que P1 établit et ce qui reste ouvert</h2>
<p><b>Établi :</b> le drapeau de P0 admet un défaut de fermeture continu en ε ; une simplification des transitions est exécutable ; l’erreur future dispose d’une borne qui inclut le couplage et chaque quantification ; les témoins exacts et les essais indépendants concordent.</p>
<p><b>Limite :</b> ces résultats n’établissent pas un avantage propre à HN, ni une filtration apprise, ni une compression utile d’un modèle de langage. Les choix d’étages viennent du modèle construit. La borne globale devient très prudente lorsque la contraction s’affaiblit ; minimiser cette borne ne minimise pas nécessairement l’erreur observée.</p>
<p>La prochaine étape utile serait de resserrer le certificat en exploitant les cellules de quantification effectivement accessibles et la sensibilité future des sorties, puis de tester une base apprise sur un modèle couplé. Cela permettrait de mesurer si l’organisation par étages apporte un bénéfice au-delà d’une allocation prédictive ordinaire.</p>

<h2>Reproduction et sources</h2>
<div class="equation">python -m pip install -r requirements.txt
python experiment.py
python check_exact.py
python build_report.py</div>
<p>Le dossier inclut code, données, erreurs par séquence, résultats CSV/JSON, témoins rationnels, rejeux exacts, contrôle de frontière de grille et une note LaTeX autonome exposant la preuve. Les hachages de distribution sont dans <code>SHA256SUMS</code>.</p>
<p><b>Contexte :</b> <a href="https://www.sciencedirect.com/science/article/pii/S1751570X13000150">Petreczky, Wisniewski et Leth, Balanced truncation for linear switched systems (2013)</a> étudient des bornes de réduction pour systèmes commutés. <a href="https://link.springer.com/article/10.1007/s10444-018-9610-z">Gosea et al. (2018)</a> développent une autre approche par Gramiens couplés. P1 n’implémente pas leur réduction équilibrée : la preuve explicite ci-dessus utilise une matrice de comparaison positive et la récurrence exacte de l’erreur.</p>
<p class="muted">Statuts : PROUVÉ sous hypothèses explicites pour la propagation et la boîte invariante ; VÉRIFIÉ EXACTEMENT pour les témoins et les rejeux rationnels ; MESURÉ pour la compression synthétique ; NON ÉTABLI pour une contribution propre à HN ou pour les modèles entraînés sur texte.</p>
</main></html>'''
(ROOT/'HN_Memory_P1_Report.html').write_text(report)
print(ROOT/'HN_Memory_P1_Report.html')
