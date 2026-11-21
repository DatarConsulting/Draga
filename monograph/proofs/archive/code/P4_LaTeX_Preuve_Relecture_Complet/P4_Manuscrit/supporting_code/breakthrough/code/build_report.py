"""Build the French execution report from recorded results, without rerunning tests."""
import json,math,statistics,html,hashlib,sys,platform,subprocess
from pathlib import Path
from fractions import Fraction as F
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle,PageBreak,Image,KeepTogether
from reportlab.lib.styles import getSampleStyleSheet,ParagraphStyle
from reportlab.lib import colors
from reportlab.lib.units import mm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
ROOT=Path(__file__).resolve().parents[1];D=ROOT/'results';OUT=ROOT/'report';OUT.mkdir(exist_ok=True)
def load(name):return json.loads((D/name).read_text())
pilot=load('pilot_summary.json');replays=[json.loads(p.read_text()) for p in sorted(D.glob('replay_*_uniform.json'))+sorted(D.glob('replay_*_mixed.json'))]
assert len(replays)==72
prob=load('probability.json');pr=prob['synthetic']+prob['trained'];tab=load('output_and_table_audit.json')['records'];higher=load('higher_rank.json')['records']
meanref=statistics.mean(c['test']['reference_BPC'] for c in prob['trained']);meanq=statistics.mean(c['test']['quantized_BPC'] for c in prob['trained'])
metrics={'rank2_instances':36,'rank2_C_python_steps':sum(r['exact_C_python_steps'] for r in replays),
 'mixed_state_savings_3_to_2_instances':sum(r['byte_saving_percent']>=25 for r in pilot['records']),
 'max_mixed_certified_error':max(r['mixed_bound'] for r in pilot['records']),
 'max_observed_error_all_policies':max(r['max_observed_abs_error'] for r in replays),
 'probabilistic_instances':len(pr),'probabilistic_C_python_steps':sum(c['checks']['C_python_steps'] for c in pr),
 'higher_rank_instances':len(higher),'higher_rank_C_python_steps':sum(t['C_python_steps'] for c in higher for t in c['trials']),
 'reference_BPC_mean':meanref,'quantized_BPC_mean':meanq,'iid_BPC':prob['trained'][0]['test']['iid_BPC'],
 'max_probability_TV_bound':max(float(F(c['total_word_TV_bound'])) for c in pr),
 'scope':'Constructed examples and small character task; no breakthrough or intrinsic optimal memory established.'}
(D/'consolidated_metrics.json').write_text(json.dumps(metrics,indent=2))
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9})
rs=pilot['records'];fig,ax=plt.subplots(figsize=(9,2.8));xx=np.arange(36)
obs={r['id']:r['max_observed_abs_error'] for r in replays if r['policy']=='mixed'}
ax.plot(xx,[r['mixed_bound'] for r in rs],'o-',ms=3,label='Borne certifiée, tous les mots',color='#0f766e')
ax.plot(xx,[obs[r['id']] for r in rs],'o-',ms=3,label='Maximum observé, référence float64',color='#ce7a22')
ax.axhline(.01,color='#b33845',ls='--',lw=1,label='Budget fixé : 0,01');ax.set_ylim(0,.0107)
ax.set_xticks([5.5,17.5,29.5],['Proche rang 2','Changement de base non normal','Frontière de stabilité'])
ax.set_ylabel('Erreur absolue scalaire');ax.spines[['top','right']].set_visible(False);ax.legend(loc='upper center',bbox_to_anchor=(.5,1.22),ncol=2,frameon=False,fontsize=8);fig.tight_layout();fig.savefig(OUT/'bounds.png',dpi=190);plt.close(fig)
pdfmetrics.registerFont(TTFont('DV','/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'))
pdfmetrics.registerFont(TTFont('DVB','/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf'))
pdfmetrics.registerFont(TTFont('DVM','/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf'))
pdfmetrics.registerFontFamily('DV',normal='DV',bold='DVB',italic='DV',boldItalic='DVB')
styles=getSampleStyleSheet()
styles.add(ParagraphStyle(name='Body',fontName='DV',fontSize=9.5,leading=14,spaceAfter=8,textColor=colors.HexColor('#183042')))
styles.add(ParagraphStyle(name='Smallx',parent=styles['Body'],fontSize=8,leading=11,spaceAfter=6))
styles.add(ParagraphStyle(name='Titlex',parent=styles['Body'],fontName='DVB',fontSize=27,leading=32,spaceAfter=16))
styles.add(ParagraphStyle(name='H1x',parent=styles['Body'],fontName='DVB',fontSize=18,leading=23,spaceAfter=14,textColor=colors.HexColor('#0f766e')))
styles.add(ParagraphStyle(name='H2x',parent=styles['Body'],fontName='DVB',fontSize=11,leading=15,spaceBefore=7,spaceAfter=7))
styles.add(ParagraphStyle(name='Equation',parent=styles['Body'],fontSize=10.2,leading=16,leftIndent=12,backColor=colors.HexColor('#eef6f4'),borderPadding=8,spaceBefore=6,spaceAfter=13))
story=[]
def p(s,style='Body'):story.append(Paragraph(s,styles[style]))
def h(s):p(s,'H2x')
def eq(s):p(s,'Equation')
def page(n,title):
    if story:story.append(PageBreak())
    p(f'EXÉCUTION • 01 NOVEMBRE 2026 • {n:02d}','Smallx');p(title,'H1x')
def table(headers,rows,widths):
    data=[[Paragraph(str(x),styles['Smallx']) for x in headers]]+[[Paragraph(str(x),styles['Smallx']) for x in row] for row in rows]
    t=Table(data,colWidths=[v*mm for v in widths],repeatRows=1,hAlign='LEFT')
    t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#dceee9')),('VALIGN',(0,0),(-1,-1),'TOP'),('LEFTPADDING',(0,0),(-1,-1),7),('RIGHTPADDING',(0,0),(-1,-1),7),('TOPPADDING',(0,0),(-1,-1),7),('BOTTOMPADDING',(0,0),(-1,-1),6),('LINEBELOW',(0,0),(-1,0),.6,colors.HexColor('#0f766e')),('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,colors.HexColor('#f4f7f8')])]))
    story.append(t);story.append(Spacer(1,8))
def link(label,url):return f'<link href="{html.escape(url,quote=True)}" color="#0f766e">{label}</link>'
page(1,'Mémoire prédictive certifiée')
p('Résultats d’exécution<br/>et bilan consolidé','Titlex')
p('<b>Résultat obtenu :</b> un dossier reproductible relie une réduction explicite de rang 4 vers 2, des certificats rationnels et une machine entière effectivement exécutée. Des variantes probabilistes et des extensions aux rangs 3 et 4 ont également été testées.')
p('<b>Le breakthrough n’est pas établi.</b> Les exemples sont construits pour être réductibles, les outils de certification ont des précédents, les bornes inférieures restent faibles et le petit essai prédictif ne bat pas son témoin sans mémoire.')
table(['Preuve et exécution','Résultat'],[
 ['Réduction signée 4 → 2','36 sources minimales de rang 4 ; 36 vérifications rationnelles réussies.'],
 ['Machine entière à deux coordonnées','52 637 040 transitions C/Python identiques ; aucun dépassement observé ; invariant valable à tout horizon.'],
 ['État récurrent','28/36 cas : 3 → 2 octets dans le format à deux coordonnées. Paquet complet : 67 → 66 octets.'],
 ['Contrôle décisif','Les 72 machines quantifiées se réencodent avec un identifiant d’un octet. Les tables occupent 69 à 699 octets au total.'],
 ['Probabilités','12 mélanges synthétiques de rang 4 et 5 modèles appris : TV sur les mots terminés ≤ 1 %, certifiée.'],
 ['Généralisation','4 cas de dimensions 5 → 3 et 6 → 4 ; 800 000 transitions comparées exactement.']],[58,116])
h('Lecture du dossier')
p('« Vérifié » désigne une vérification arithmétique ordinaire et des preuves explicites, pas une formalisation dans un assistant de preuve. Les résultats hérités sont identifiés séparément. Les tests C ne prouvent pas à eux seuls un contrat universel : celui-ci repose sur les invariants rationnels.')
p('Périmètre : plan retrouvé, manuscrit joint, six rapports antérieurs et deux archives expérimentales accessibles. Ce bilan ne prétend pas couvrir chaque échange ou chaque fichier des deux dernières semaines.','Smallx')
page(2,'Ce qui a été exécuté')
table(['Jalon','Livraison','Statut scientifique'],[
 ['P0 — Positionnement','Protocole figé, audit ciblé de précédents, métriques séparées.','Nouveauté non démontrée.'],
 ['P1 — Réduction','36 témoins rationnels, minimalité exacte, Gramiens et borne supérieure.','Réduction construite à partir d’un sous-modèle connu ; pas de solveur général.'],
 ['P2 — Coefficients','Bases native et équilibrée ; coefficients dyadiques à 24 bits fractionnaires ; recertification.','36/36 sources acceptées dans les deux bases.'],
 ['P3 — Exécution','État empaqueté, noyau C, arrondi pair, sortie dyadique exacte, invariants et replays.','Contrat scalaire à tout horizon établi sous les hypothèses documentées.'],
 ['P4 — Mémoire','Allocation mixte, uniforme, équilibrage, FP16/FP32 et tables ; témoins de séparation.','Gain de format constaté ; optimum mémoire non établi.'],
 ['P5 — Probabilités','Normalisation positive, TV globale, postérieurs finis, 5 entraînements frais ; lissage local KL exécuté.','TV ≤ 1 % et KL certifiés pour la machine lissée ; aucun gain prédictif.'],
 ['P6 — Extension','Preuve sous condition de majoration pour tout k ; essais k=3,4 ; dossier.','Pas de théorème universel de compression ni d’article à nouveauté établie.']],[29,76,69])
h('Contrat resté fixe')
p('Alphabet binaire pour la branche signée. Norme Hankel de la source ≤ 1, certifiée par sa norme de Hilbert-Schmidt. Budget scalaire 0,01 : réduction ≤ 0,005 ; coefficients ≤ 0,002 ; exécution ≤ 0,003. Aucun assouplissement de ces seuils après test.')
h('Écarts au programme initial')
p('Le sous-modèle est fourni par la construction, puis masqué par un changement de base rationnel ; sa découverte depuis une source arbitraire n’a pas été résolue. Une borne de Hilbert-Schmidt, plus conservatrice, remplace l’optimisation directe de la norme opérateur. Les valeurs singulières de sections finies sont des diagnostics flottants, pas des bornes certifiées.')
p('La famille « non normale » est similaire à la famille « proche rang 2 » et produit les mêmes séries. Elle teste le conditionnement, pas une troisième population indépendante. Les « graines » choisissent trois paramètres déterministes. Le calendrier de trente jours du plan était indicatif ; ce dossier constitue un cycle d’exécution, pas trente jours de recherche simulés.')
page(3,'Le certificat mathématique obtenu')
h('1. Réduction et coefficients stockés')
p('Pour f(w)=αA<sub>w</sub>β, poser P=ββᵀ+Σ<sub>a</sub>A<sub>a</sub>PA<sub>a</sub>ᵀ et Q=αᵀα+Σ<sub>a</sub>A<sub>a</sub>ᵀQA<sub>a</sub>. La stabilité de comptage est certifiée par X−Σ<sub>a</sub>A<sub>a</sub>XA<sub>a</sub>ᵀ=I avec X strictement positif. Les identités sont vérifiées en nombres rationnels.')
eq('‖H<sub>f</sub>‖² ≤ ‖H<sub>f</sub>‖²<sub>HS</sub> = tr(PQ).')
p('La même construction est appliquée au bloc résiduel et à la réalisation jointe de la différence avant/après quantification. La racine carrée est majorée par un dyadique calculé par arithmétique entière. Le vérificateur séparé contrôle les identités, les signes et les déterminants de minimalité.')
h('2. Écritures récurrentes arrondies')
p('Convention ligne : x<sub>t+1</sub>=x<sub>t</sub>A<sub>a</sub>. Pour les matrices stockées, M<sub>ij</sub>=max<sub>a</sub>|A<sub>a,ij</sub>|. Une inverse (I−Mᵀ)<sup>−1</sup> positive fournit un témoin de contraction. Avec pas h<sub>i</sub> et η<sub>i</sub>=h<sub>i</sub>/2 :')
eq('E = (I−Mᵀ)<sup>−1</sup>η ;<br/>R = (I−Mᵀ)<sup>−1</sup>(|x̂<sub>0</sub>|+η).')
p('L’arrondi initial vérifie |x̂<sub>0</sub>−α|≤η. Par récurrence, |x̂<sub>t</sub>−x<sub>t</sub>|ᵀ≤E, car MᵀE+η=E. De même, |x̂<sub>t</sub>|ᵀ≤R. Vérifier R dans l’intervalle représentable interdit le débordement d’état ; un calcul distinct borne chaque accumulateur entier.')
eq('sup<sub>w</sub> |f(w)−g<sub>machine</sub>(w)|<br/>≤ U<sub>réduction</sub> + U<sub>coefficients</sub> + |β̂|ᵀE ≤ 0,01.')
h('Portée exacte')
p('Cette dernière borne est ponctuelle, uniforme sur les mots et sans limite d’horizon. Elle ne certifie pas la norme Hankel de la machine arrondie, qui n’est plus une réalisation linéaire de rang 2. Elle ne donne aucune garantie TV/KL pour une série signée. Les coefficients, l’ordre des opérations, l’arrondi et le format de sortie font partie du contrat.')
page(4,'Construction explicite et code exécuté')
h('Un témoin complet : near_rank2_m0_s101')
eq('A<sub>0</sub> = [[7/25, 7/200], [7/640, 21/100]]<br/>A<sub>1</sub> = [[7/40, 7/100], [7/320, 7/25]]<br/>α = [1/2, 1/4] ; β = [1/2, 1/16]ᵀ.')
p('La source à quatre dimensions ajoute le bloc R<sub>0</sub>=[[1/5,1/11],[0,2/7]], R<sub>1</sub>=[[1/7,0],[1/13,1/4]], α<sub>R</sub>=[1/2,1/3], β<sub>R</sub>=[3/10000,1/13000]ᵀ. Un changement de base rationnel dense est appliqué à la source. Des mineurs non nuls prouvent son rang minimal 4.')
example=json.loads((ROOT/'models/near_rank2_m0_s101.json').read_text());mix=example['best']['mixed'];f=example['frames'][mix['basis']]
table(['Quantité','Valeur du témoin'],[
 ['Erreur de réduction majorée',f"{float(F(example['model_error_bound'])):.10f}"],
 ['Erreur des coefficients majorée',f"{float(F(f['coefficient_error_bound'])):.10f}"],
 ['Erreur récurrente majorée',f"{float(F(mix['execution_bound'])):.10f}"],
 ['Borne totale',f"{float(F(mix['total_bound'])):.10f} < 0,01"],
 ['État uniforme / mixte','18 bits / 14 bits ; 3 octets / 2 octets'],
 ['Allocation mixte en base équilibrée','4 + 10 bits ; pas 1/32 et 1/256 ; état initial entier [0, −132]']],[80,94])
h('Exécution concrète')
p('Le noyau C lit et réécrit uniquement un tableau d’octets empaqueté. Les coefficients sont des int32 dyadiques ; les produits utilisent des entiers signés de 128 bits ; l’arrondi est au plus proche, égalités vers le pair. La sortie est un numérateur int64 et un exposant dyadique, vérifiés exactement sur chaque état atteignable.')
p('Pour chacune des 72 configurations : tous les préfixes binaires jusqu’à longueur 16, puis six séquences de 100 000 symboles. Les 731 070 comparaisons par configuration totalisent 52 637 040 égalités C/Python. La fermeture complète des états atteignables est également calculée. Trois certificats volontairement altérés sont rejetés.')
p('Mesure de compilation locale : pile statique annoncée de 208 octets pour step et 56 octets pour output_value. Ce n’est ni la mémoire du processus Python ni un gain RSS. Les compilateurs et le matériel ne sont pas formellement vérifiés.','Smallx')
page(5,'Mémoire : gain de format, limite intrinsèque')
table(['Famille','Cas 3 → 2 octets','Interprétation'],[
 ['Proche rang 2','8 / 12','4 cas restent à 3 octets.'],
 ['Non normale','8 / 12','Même série, autre base ; pas un échantillon indépendant.'],
 ['Frontière de stabilité','12 / 12','Gain constaté sur les 12 paramètres.'],
 ['Total','28 / 36','−33,3 % d’état ; −1,49 % du paquet complet dans ces cas.']],[60,40,74])
story.append(Image(str(OUT/'bounds.png'),width=174*mm,height=54*mm))
p(f"Maximum de la borne mixte : {metrics['max_mixed_certified_error']:.8f}. Maximum observé toutes configurations : {metrics['max_observed_error_all_policies']:.8f}. Ce second nombre utilise une référence float64 et ne remplace pas la preuve.",'Smallx')
h('Le témoin par table empêche une conclusion d’optimalité')
p('La fermeture des machines mixtes contient 6 à 38 états comportementalement distincts ; les témoins uniformes en contiennent 16 à 69. La minimisation de Moore et le contrôle de toutes les arêtes permettent un identifiant persistant d’un octet dans les deux cas. La sortie et toutes les transitions restent exactement identiques à celles de la machine quantifiée.')
p('Le coût total augmente : 69 à 389 octets pour les tables mixtes et 169 à 699 pour les tables uniformes, contre 66 ou 67 octets pour les réalisations à deux coordonnées. Le choix dépend donc du compromis entre mémoire par flux et paramètres partagés. Le gain mixte ne subsiste pas en octets d’état face à ce témoin admissible.')
h('Borne inférieure et contrôles')
p('Des préfixes dont les continuations diffèrent de plus de 2ε imposent des états distincts à toute machine déterministe autonome de précision ε. Les témoins rationnels trouvés séparent 4 à 7 préfixes : seulement 2 à 3 bits nécessaires, contre 12 à 20 utilisés par l’allocation mixte. L’écart n’est pas fermé.')
p('Les comparateurs FP16/FP32 ont été exécutés dans les deux bases : 144 cas, 4/8 octets d’état. Ils restent des diagnostics NumPy sur horizon fini ; aucune preuve universelle en virgule flottante n’est annoncée.','Smallx')
page(6,'Probabilités : une extension distincte et certifiée')
p('Douze sources positives sont des mélanges de deux générateurs à deux états : p=(1−δ)p<sub>2</sub>+δr<sub>2</sub>, δ=1/1000. Leur rang 4 est vérifié exactement. Garder p<sub>2</sub> coûte au plus δ en TV sur les mots terminés. Chaque état possède une probabilité d’arrêt strictement positive.')
h('Postérieur fini sur une grille')
p('Le postérieur des deux états est représenté par x=j/2<sup>f</sup>. Pour chaque lettre, sa mise à jour est une fonction homographique. Si A<sub>a</sub> a pour sommes de lignes s<sub>0</sub>,s<sub>1</sub>, une constante vérifiable est :')
eq('q = max<sub>a</sub> |det A<sub>a</sub>| / min(s<sub>0</sub>,s<sub>1</sub>)² &lt; 1 ;<br/>e ≤ [2<sup>−f−1</sup>]/(1−q).')
p('Les coefficients sont quantifiés par ligne en préservant exactement positivité et somme 1. Soit d<sub>ligne</sub> leur erreur TV maximale, τ<sub>0</sub> l’arrêt minimal avant quantification, τ après, et C la distance TV entre les distributions de prochain événement des deux états. Un couplage pas à pas donne :')
eq('TV(p,p<sub>machine</sub>) ≤ δ + d<sub>ligne</sub>/τ<sub>0</sub> + Ce/τ.')
p('Les événements sont les lettres et l’arrêt. Les sorties C sont des numérateurs entiers avec un dénominateur commun ; leur somme est exactement le dénominateur. La mise à jour du postérieur utilise une division entière arrondie. Tous les états atteignables sont vérifiés ; 1 700 673 transitions supplémentaires sont comparées C/Python.')
table(['Lot','Borne TV','État / paquet'],[
 ['12 sources synthétiques','0,7757 % à 0,9025 %','1 ou 2 octets / 49 ou 50 octets'],
 ['5 modèles appris','0,5317 % à 0,9785 %','1 octet / 737 octets']],[53,57,64])
h('KL : réalisation par lissage local, à état inchangé')
p('Le mélange global q<sub>+</sub>=(1−d)q+d u, étudié dans la note antérieure, n’est pas déployé. Une variante locale est exécutée : chaque distribution de prochain événement devient (1−λ)q+λu, où u arrête avec probabilité 1/2 et émet chaque lettre avec probabilité 1/(2m). Les transitions d’état restent celles de q.')
eq('TV(p,q<sub>λ</sub>) ≤ d + λ/τ ;<br/>KL(p‖q<sub>λ</sub>) ≤ −log(1−d) − log(1−λ)/τ<sub>0</sub><br/>+ d log(2/λ) + Ψ<sub>p</sub>(d) log(2m/λ).')
p('λ est dyadique et choisi pour conserver TV ≤ 1 %. Les 17 sorties lissées sont vérifiées exactement en C. Un octet statique est ajouté, aucun état récurrent. Les bornes KL rationnelles vont de 2,167 à 12,601 nats : elles sont finies mais très conservatrices. La preuve et le comptage sont détaillés en annexe technique.','Smallx')
page(7,'Essai appris : précision conservée, utilité non gagnée')
h('Données supplémentaires et protocole')
p('Le pilote antérieur utilisait les 180 000 premiers caractères normalisés des Malheurs de Sophie. Ce nouveau test utilise les 12 755 caractères restants du même livre : 8 000 pour apprendre, 2 000 pour valider, 2 755 pour tester. Aucune fenêtre de test de 64 caractères ne correspond exactement à une fenêtre de l’ancien corpus ou du nouvel apprentissage.')
p('Il s’agit d’une portion supplémentaire du même ouvrage, pas d’un nouvel auteur. L’alphabet est appris sur le train : 45 symboles avec &lt;UNK&gt;. Les BPC concernent cet encodage. Le test n’a servi ni à choisir les graines, ni à régler les hyperparamètres.')
p('Cinq graines fixes, 30 itérations d’une heuristique alternée de type Baum-Welch, fenêtres de 64 symboles. Les transitions et émissions sont contraintes pour garder une contraction vérifiable. L’arrêt est fixé à 0,1. Les paramètres appris sont rationalisés et normalisés avant de définir le modèle de référence.')
table(['Graines','Référence BPC','Machine BPC','TV conditionnelle max observée'],[
 [c['id'].replace('fresh_char_seed',''),f"{c['test']['reference_BPC']:.6f}",f"{c['test']['quantized_BPC']:.6f}",f"{100*c['test']['observed_max_conditional_TV']:.4f} %"] for c in prob['trained']],[25,43,43,63])
table(['Comparaison','Valeur'],[
 ['Moyenne de référence',f'{meanref:.6f} BPC'],
 ['Moyenne quantifiée',f'{meanq:.6f} BPC'],
 ['Témoin unigramme sans mémoire',f"{prob['trained'][0]['test']['iid_BPC']:.6f} BPC"],
 ['Plus grand écart absolu de quantification','0,00007143 BPC'],
 ['Supériorité prédictive démontrée','Non : chacun des cinq modèles reste derrière le témoin unigramme.']],[95,79])
h('Interprétation')
p('La certification de la réalisation numérique est informative, mais le modèle contraint apprend très peu de mémoire utile. Les machines atteignables n’ont que 1 à 3 états de postérieur ; l’une se comporte sans mémoire après quantification. Ce résultat explique pourquoi la petite taille ne suffit pas à démontrer une compression prédictive intéressante.')
p('La TV certifiée compare les lois des mots avec arrêt du modèle de référence et de sa machine. Elle ne mesure pas l’écart à la vraie distribution du français. Les TV conditionnelles du tableau sont uniquement observées sur le test.','Smallx')
page(8,'Généralisation : ce qui passe, ce qui échoue')
h('L’invariant s’étend à tout k sous les mêmes hypothèses')
p('La preuve de la page 3 ne dépend pas de k=2. Pour k coordonnées, il suffit de vérifier M≥|A<sub>a</sub>|, une inverse positive (I−Mᵀ)<sup>−1</sup>, l’arrondi initial et les marges représentables. En posant wᵀ=|β|ᵀ(I−Mᵀ)<sup>−1</sup> et h<sub>i</sub>=2<sup>eᵢ−pᵢ+1</sup> :')
eq('ε<sub>exécution</sub> ≤ Σ<sub>i=1…k</sub> w<sub>i</sub> 2<sup>eᵢ−pᵢ</sup> ;<br/>B<sub>état</sub> = ⌈(Σ<sub>i</sub>p<sub>i</sub>)/8⌉ octets.')
p('L’amplification dépend de la résolvante, de β et des plages e<sub>i</sub>. Elle peut diverger à l’approche de la frontière. Les coefficients et les multiplications croissent avec |Σ|k² ; la preuve n’implique pas un coût polynomial uniforme de découverte du meilleur modèle.')
table(['Réduction','Échelle','Bits uniformes / mixtes','Octets uniformes / mixtes'],[
 [f"{r['n']} → {r['k']}",r['scale'],f"{r['trials'][0]['allocation']['state_bits']} / {r['trials'][1]['allocation']['state_bits']}",f"{r['trials'][0]['allocation']['state_bytes']} / {r['trials'][1]['allocation']['state_bytes']}"] for r in higher],[36,30,54,54])
p('Quatre exemples dyadiques avec bloc faible connu : minimalité, Gramiens et invariants vérifiés ; 100 000 transitions C/Python par politique et par exemple. Le choix mixte utilise une programmation dynamique sur une borne séparable, puis une vérification de faisabilité ; aucun optimum universel n’est revendiqué.')
h('Un contre-exemple à la nécessité de notre condition')
p('Prendre A<sub>0</sub>=(3/5)I et A<sub>1</sub>=(3/5)S, S permutation des deux coordonnées. La matrice X=(25/7)I satisfait exactement X−ΣA<sub>a</sub>XA<sub>a</sub>ᵀ=I : le système est stable au sens du comptage. Pourtant M=(3/5)·1 possède un rayon spectral 6/5 ; notre certificat par maximum entrée par entrée rejette le système.')
p('C’est un échec de la condition suffisante, pas une impossibilité de précision finie. La norme infinie de chaque transition contracte ici d’un facteur 3/5. Pour l’échelle 1 ou 4/5, au contraire, l’identité de trace (1−2c²)tr(X)=2 interdit un témoin de Stein positif. Les cas de frontière restent dans le dossier.')
h('Ce qui ne suit pas du rang')
p('Un rang linéaire k ne garantit ni une petite réalisation positive, ni une petite machine entière, ni une queue de longueur favorable à KL. Les séries signées, les lois positives et les machines arrondies conservent des contrats distincts.')
page(9,'Résultats antérieurs : expérience et compression')
p('Les chiffres ci-dessous proviennent des rapports accessibles. Seuls les contrôles explicitement indiqués ont été relancés dans cette exécution.')
table(['Axe / provenance','Résultat rapporté','Statut ici'],[
 ['Pilote multilettre 256 B [D1]','24 entraînements, 3 graines. BPC 2,820 contre 2,884 ; variance élevée sur l’ordre. Test réutilisé et stabilité de comptage non acquise.','Relu ; pas de nouveau gain validé.'],
 ['Pilote compatible certificat [D2]','50 entraînements, 5 graines. Témoin 2,486 ± 0,004 BPC ; modèle certifié 2,499 ± 0,018. Critère conjoint non atteint. Ordre long : 0 % pour multi_cert.','10 certificats / 40 têtes revérifiés exactement ; contrôle falsifié rejeté.'],
 ['Réduction 256 → 192 B [D2]','Perte de 0,2334 BPC ; TV observée environ 16–23 %. Les bornes TV texte atteignent 100 %.','Aucune compression prédictive certifiée utile établie.'],
 ['Filtration approximative HN [D3]','2 160 configurations publiées. Allocation calibrée souvent ≈90–95 % de MSE en moins à budget égal ; allocation fondée sur la borne parfois moins bonne.','2 160 lignes recalculées ; 4 608 contrôles exacts ponctuels et 720 contrôles futurs.'],
 ['Codec GDN sans perte [D4]','18,646 % du budget nominal hors ligne ; ≈12,07 % pour le codec en ligne plus simple. 1 380 sorties logits avec empreintes identiques.','Relu seulement ; cache natif conservé, pas de gain RSS ou de vitesse prouvé.'],
 ['Codec KV Qwen [D5]','Deux textes de 160 tokens, une couche. TV moyenne 0,3925 %, maximum 6,6294 % ; gain cache 19,76 % face au FP32.','Relu seulement ; certificat inutilisable, cache natif conservé, pas de comparaison favorable FP16.']],[44,88,42])
h('Ce que la nouvelle exécution ajoute')
p('La réalisation entière et son invariant d’erreur ont maintenant été exercés sur un banc contrôlé. Cela ne rétrocertifie ni le cache Qwen, ni les anciens modèles multilettres. La preuve du codec GDN vient de la réversibilité binaire, indépendamment de la théorie de Hankel.')
p('Dans HN, les quatre coordonnées sont conservées et six coefficients de transition sont supprimés. Le drapeau est hérité d’une expérience précédente ; il n’est pas démontré comme filtration HN du modèle perturbé. Le suivi en ligne de la borne a posteriori ajouterait 32 octets par flux. Aucun bénéfice attribuable spécifiquement à HN n’est établi.','Smallx')
page(10,'Résultats antérieurs : mathématiques et nouveauté')
table(['Document','Apport du texte','Limite conservée'],[
 ['Exact Certificates for Multiletter Hankel [D6]','Th. 3.1 : existence d’une approximation optimale de rang borné. Th. 4.1 : certificats par inégalités de Stein. Cor. 5.3 : synthèse algébrique par tranches fixes.','Ne fournit pas automatiquement un exécutable en précision finie, une réalisation positive ou TV/KL.'],
 ['Sharp Hankel-to-KL [D7]','Note de rang 1 : module asymptotique sous queue exponentielle, Θ(εᵞ log²(1/ε)) ; obstruction avec seule moyenne de longueur bornée.','Résultat hérité de la note ; preuve complète et priorité non réexpertisées ici.'],
 ['Occupation Flows [D8]','Dictionnaire positif fixé ; programme linéaire de flux ; TV ≤ J/2 ; arrêt propre ; contrôle KL après lissage et hypothèse de queue.','Qualité et découverte du dictionnaire restent ouvertes ; rang linéaire et nombre d’états positifs distincts.']],[48,81,45])
h('Positionnement par rapport aux précédents')
p('La vérification des erreurs d’arrondi, des débordements et des cycles limites de filtres quantifiés, y compris à horizon illimité, est documentée chez Cox, Sankaranarayanan et Chang (2014). Le seul ajout d’un certificat de précision finie ne suffit donc pas à établir la nouveauté [W1].')
p('Erbas, Intes et Pandey (2026) étudient déjà l’effet de l’écriture récurrente quantifiée et des mémoires de correction. Notre contrôle de cette écriture doit être positionné face à ces travaux [W2]. Glaude et Pietquin (2016) donnent aussi des méthodes d’apprentissage d’automates probabilistes propres avec garanties PAC [W3].')
h('Verdict scientifique')
p('<b>Résultat identifiable :</b> une chaîne explicite et reproductible de témoins rationnels jusqu’aux octets exécutés, avec distinction des trois erreurs et des différentes notions de mémoire. Les contrôles par tables et les essais appris rendent la limite visible.')
p('<b>Résultat non établi :</b> une séparation nouvelle taille–erreur, une borne inférieure proche du coût atteint, une réduction générale découverte depuis des données, ou un avantage prédictif / système sur un modèle utile. Le manuscrit joint et les notes antérieures restent des contributions candidates, sans validation externe déduite de cette session.')
p('L’audit bibliographique est ciblé, pas exhaustif. Il ne permet ni de certifier une priorité ni d’affirmer qu’aucun travail concurrent n’existe. Les calculs de cette session ne sont pas un processus de revue par les pairs.','Smallx')
page(11,'Reproduire, vérifier, puis décider')
h('Contenu du dossier reproductible')
p('<b>code/</b> contient les générateurs, les deux vérificateurs rationnels indépendants, les noyaux C, les replays, les contrôles flottants et les audits par tables. <b>models/</b> contient les témoins JSON et les paquets binaires. <b>results/</b> contient les résultats complets, dont les échecs et les bornes. <b>prior_evidence/</b> conserve les archives relancées, les extraits lus et les journaux. Le manuscrit joint et le plan figé sont inclus.')
p('Vérification courte : python code/verify.py ; python code/verify_probability.py. Reconstruction complète : suivre RUN_ORDER.txt, sans Python -O, avec gcc, NumPy, SciPy et SymPy. Le manifeste SHA256SUMS couvre les fichiers livrés. Les bibliothèques .so sont recompilables et dépendent de la plateforme.','Smallx')
h('Décision proposée à partir des résultats')
p('Conserver le banc de preuve comme outil. Pour viser un résultat nouveau, la priorité est une famille difficile dont la réduction n’est pas donnée, accompagnée d’une borne inférieure en bits qui tienne compte des tables et des paramètres partagés. Le test appris doit ensuite montrer une mémoire utile avant d’annoncer une compression prédictive.')
p('Restent ouverts : découverte du sous-modèle, gap d’optimalité mémoire, optimisation positive générale, garantie KL plus informative, transfert à un modèle entraîné utile, mesure de latence et de mémoire résidente. Ces points ne sont pas remplacés par des promesses de breakthrough.')
h('Sources et provenance')
refs=[
 ('D1 — Pilote multilettre 256 B','libfile_54a80b4ef3d88191a4534c18393d2aad'),
 ('D2 — Pilote compatible certificat','libfile_ee3d16d856848191bf3be8e7b0134380'),
 ('D3 — Archive HN P1','libfile_5334700fafb88191a4568f8b8f8efdad'),
 ('D4 — Rapport GDN sans perte','libfile_d12804e781f4819190337638d0a5ccda'),
 ('D5 — Rapport KV Qwen','libfile_79de102a081c8191b7f83b50fde7f22b'),
 ('D6 — Manuscrit Exact Certificates joint','libfile_81a98554ca648191a18454ae83f4b98c'),
 ('D7 — Sharp Hankel-to-KL','libfile_896a205d9c50819198e40adfd7050fa7'),
 ('D8 — Occupation Flows','libfile_8d912961cec48191bf38aa7890854d79')]
for label,ref in refs:p(link(label,'https://chatgpt.com/library?file_id='+ref+'&account_id=personal'),'Smallx')
p(link('W1 — Cox et al. (2014), A Bit Too Precise?','https://plv.colorado.edu/bec/papers/filters-sttt13.html'),'Smallx')
p(link('W2 — Erbas et al. (2026), When Quantization Breaks Memory','https://arxiv.org/abs/2609.04490'),'Smallx')
p(link('W3 — Glaude et Pietquin (2016), PAC learning of Probabilistic Automaton','https://proceedings.mlr.press/v48/glaude16.html'),'Smallx')
p('Sources web consultées le 1 novembre 2026. Les nombres des expériences nouvelles sont produits par les JSON du dossier ; les chiffres historiques sont attribués aux documents D1–D8.','Smallx')
page(12,'Annexe : pourquoi le lissage KL est certifié')
p('Soient p la loi de référence et q la machine finie non lissée, avec TV(p,q)≤d&lt;1. La machine q<sub>λ</sub> garde exactement le même état et la même mise à jour après chaque lettre ; seules ses probabilités de prochain événement sont lissées. Toutes les lois sont propres grâce à l’arrêt positif.')
h('Comparaison des probabilités de mots')
p('Pour un mot de longueur L, le produit des L émissions et de l’arrêt donne les deux minorations suivantes, valables sur chaque chemin déterministe de la machine :')
eq('q<sub>λ</sub>(w) ≥ (1−λ)<sup>L+1</sup> q(w) ;<br/>q<sub>λ</sub>(w) ≥ (λ/(2m))<sup>L</sup>·λ/2.')
h('Décomposition de la masse commune')
p('La masse commune de p et q vaut au moins 1−d. On peut donc écrire p=(1−d)p<sub>0</sub>+dp<sub>1</sub> avec (1−d)p<sub>0</sub>≤q et dp<sub>1</sub>≤p. Par convexité de KL dans son premier argument, on borne séparément ces deux parties.')
p('Pour la partie commune, log(p<sub>0</sub>/q<sub>λ</sub>) ≤ −log(1−d) − (L+1)log(1−λ). Sa contribution pondérée est au plus −log(1−d) − E<sub>p</sub>(L+1)log(1−λ). Pour la partie résiduelle, on majore KL par la cross-entropie et on utilise la seconde minoration : sa contribution est au plus d log(2/λ)+Ψ<sub>p</sub>(d)log(2m/λ).')
p('Si l’arrêt de la source vaut au moins τ<sub>0</sub>, E(L+1)≤1/τ<sub>0</sub>. Pour une sous-mesure η≤p de masse d, E<sub>η</sub>L=Σ<sub>j≥1</sub>η(L≥j)≤Σ min(d,r<sup>j</sup>), r=1−τ<sub>0</sub>. Pour tout entier N≥0 :')
eq('Ψ<sub>p</sub>(d) ≤ dN + r<sup>N+1</sup>/τ<sub>0</sub>.')
h('Arithmétique et coût réellement comptés')
p('Le calcul prend d comme majorant dyadique de la TV certifiée. Les logarithmes sont majorés par la série rationnelle de 2 atanh((z−1)/(z+1)), après réduction z=2ᵏt avec 1≤t&lt;2. Le reste de série est majoré explicitement, puis la borne finale est arrondie vers le haut. final_audit.py recalcule indépendamment ces inégalités.')
p('Le noyau C émet des numérateurs et un dénominateur entiers communs. Chaque valeur intermédiaire tient dans uint64 sur les cas acceptés. Le paquet lissé ajoute un octet pour l’exposant de λ : 50/51 octets en synthétique et 738 octets pour le texte, dont le même 1 ou 2 octets d’état. Les tableaux de sortie et les multiplications sont des temporaires, pas une réduction de mémoire totale du processus.')
p('Cette preuve contrôle KL entre le modèle de référence et sa machine lissée. Elle ne certifie ni la qualité du modèle de référence face aux données ni une perte de test maximale. Le contrôle TV du lissage vient du couplage : probabilité de divergence par visite ≤λ, nombre moyen de visites de q ≤1/τ.','Smallx')
def footer(canvas,doc):
    canvas.setStrokeColor(colors.HexColor('#c8d8d8'));canvas.line(18*mm,17*mm,192*mm,17*mm)
    canvas.setFont('DV',7.5);canvas.setFillColor(colors.HexColor('#567078'));canvas.drawString(18*mm,12*mm,'Mémoire prédictive certifiée • Rapport d’exécution • 01/11/2026');canvas.drawRightString(192*mm,12*mm,str(doc.page))
pdf=OUT/'Bilan_Execution_Memoire_Certifiee.pdf'
doc=SimpleDocTemplate(str(pdf),pagesize=(210*mm,297*mm),rightMargin=18*mm,leftMargin=18*mm,topMargin=16*mm,bottomMargin=24*mm,title='Mémoire prédictive certifiée — Exécution et bilan consolidé',author='Dossier de recherche')
doc.build(story,onFirstPage=footer,onLaterPages=footer)
print(pdf)
