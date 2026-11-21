from pathlib import Path
import re, json, shutil, hashlib

ROOT=Path(__file__).resolve().parent
WORK=Path('/workspace/scratch/8e6d2ecb2b28')
SPECS=[
 ('p4_resolution','theorem.tex','Optimum exact pour une série amortie','P4_Optimum_Exact_Memoire.pdf','Optimum en bits sur une suite de précisions. Source scalaire amortie, distincte du flux conditionnel des chapitres suivants.'),
 ('p4_conditional','P4_Prediction_Conditionnelle.tex','Prédiction conditionnelle et huit optima finis','P4_Prediction_Conditionnelle.pdf','Modèle HMM central; tables admises et comptées. Huit budgets fermés, pas une formule à toute précision.'),
 ('p4_pressure','P4_Loi_Memoire_Sans_Table.tex','Loi de mémoire par pression géométrique','P4_Loi_Memoire_Sans_Table.pdf','Première preuve de la loi asymptotique pour le modèle central. Les ressources et le codage sont affinés dans les chapitres suivants.'),
 ('p4_direct','P4_Identifiant_Direct.tex','Calcul d’identifiant sans table des états','P4_Identifiant_Direct.pdf','Étape historique du calcul sans table. Le terme direct ne signifie pas une mise à jour en temps constant ni sans parcours de chemin.'),
 ('p4_uniform','P4_Codage_Uniforme.tex','Codage uniforme et poids de transfert','P4_Codage_Uniforme.pdf','La preuve historique utilise une condition suffisante auxiliaire, remplacée ensuite par le cône projectif. Les constantes de codage sont réduites au dernier chapitre.'),
 ('p4_arithmetic','P4_Memoire_Lineaire.tex','Arithmétique bornée et mémoire totale linéaire','P4_Memoire_Lineaire.pdf','Étape du contrôle des temporaires. Les budgets numériques et leurs obligations de validation sont à lire avec la relecture actuelle.'),
 ('p4_robust','P4_Compilation_Robuste.tex','Compilation robuste des paramètres','P4_Compilation_Robuste.pdf','Famille structurée et perturbations. La classe paramétrique générale et le terme de coût d’entrée sont précisés ultérieurement.'),
 ('p4_parametric','P4_Theoreme_Parametrique.tex','Théorème paramétrique sous certificat','P4_Theoreme_Parametrique.pdf','Preuve générale proposée. L’ancienne condition sur tau est retirée au chapitre suivant; l’optimalité composée reste soumise à validation indépendante.'),
 ('p4_cone_audit','P0_P4_Audit_Sans_Tau.tex','Audit P0 et suppression de la condition auxiliaire','P0_P4_Audit_Sans_Tau.pdf','Version sans tau. Les contractions primale et duale subsistent. Les mentions d’accès bibliographique sont mises à jour par la comparaison directe.'),
 ('p4_next_step','P0_P4_Comparaison_Directe.tex','Grille directe, comparaison et positionnement','P0_P4_Comparaison_Directe.pdf','La grille donne une formule sans chemin avec coefficient de mémoire un. Les coûts de l’ancien codec et la prochaine cible décrivent l’état avant la réduction du surcoût.'),
 ('p4_reduced','P4_Reduction_Surcout_Certifie.tex','Réduction certifiée du surcoût','P4_Reduction_Surcout_Certifie.pdf','Dernier état de l’implémentation: offset de trois bits pour seize instances. Les garanties globales héritées conservent le statut défini dans la relecture.'),
]

def flatten_inputs(text, folder):
 def expand(m):
  name=m.group(1)
  chunks=Path(name).parts
  archived=ROOT/'historical_tex'/chunks[0]/chunks[-1] if chunks and chunks[0].startswith('p4_') else folder/chunks[-1]
  choices=[archived,folder/name,folder/chunks[-1]]
  path=next((p for p in choices if p.exists()),None)
  if path is None:raise FileNotFoundError((name,choices))
  return flatten_inputs(path.read_text(),path.parent)
 return re.sub(r'\\input\{([^}]+)\}',expand,text)

def macros(pre):
 found=[]
 for match in re.finditer(r'\\newcommand',pre):
  pos=match.start();i=match.end()
  while i<len(pre) and pre[i].isspace():i+=1
  if pre[i]!='{':raise ValueError(pre[pos:pos+60])
  # command-name argument
  depth=1;i+=1
  while depth:
   if pre[i]=='{':depth+=1
   elif pre[i]=='}':depth-=1
   i+=1
  while i<len(pre) and pre[i].isspace():i+=1
  if pre[i]=='[':
   i=pre.index(']',i)+1
  while i<len(pre) and pre[i].isspace():i+=1
  if pre[i]!='{':raise ValueError(pre[pos:pos+60])
  depth=1;i+=1
  while depth:
   if pre[i]=='{' and pre[i-1]!='\\':depth+=1
   elif pre[i]=='}' and pre[i-1]!='\\':depth-=1
   i+=1
  found.append(pre[pos:i].replace('\\newcommand','\\providecommand',1))
 return '\n'.join(found)

PREAMBLE=r'''\documentclass[11pt,a4paper,oneside]{report}
\usepackage[T1]{fontenc}
\usepackage[utf8]{inputenc}
\usepackage[provide=*,french]{babel}
\usepackage{lmodern,amsmath,amssymb,amsthm,mathtools,booktabs,array,longtable}
\usepackage[margin=20mm,headheight=15pt]{geometry}
\usepackage{microtype,xcolor,fancyhdr,graphicx,pdfpages,etoolbox}
\usepackage[colorlinks=true,linkcolor=teal,urlcolor=teal,citecolor=teal,bookmarksnumbered=true,hypertexnames=false]{hyperref}
\definecolor{teal}{RGB}{15,102,102}
\definecolor{pale}{RGB}{238,246,246}
\hypersetup{pdftitle={P4 - Résultats, preuves et relecture critique},pdfauthor={Dossier de recherche consolidé},pdfsubject={Mémoire prédictive certifiée; rang deux; précision finie}}
\pagestyle{fancy}\fancyhf{}
\fancyhead[L]{\small Mémoire prédictive certifiée}
\fancyhead[R]{\small P0--P6 / P4 • 3 novembre 2026}
\fancyfoot[L]{\small Édition critique • preuve proposée et résultats attribués}
\fancyfoot[R]{\thepage}
\fancypagestyle{plain}{\fancyhf{}\fancyfoot[L]{\small Mémoire prédictive certifiée}\fancyfoot[R]{\thepage}}
\setlength{\parindent}{0pt}\setlength{\parskip}{5pt}
\setlength{\emergencystretch}{2em}
\setcounter{tocdepth}{0}\setcounter{secnumdepth}{2}
\patchcmd{\thebibliography}{\chapter*}{\section*}{}{}
\newtheorem{theorem}{Théorème}[chapter]
\newtheorem{lemma}[theorem]{Lemme}
\newtheorem{proposition}[theorem]{Proposition}
\newtheorem{corollary}[theorem]{Corollaire}
\newcommand{\eps}{\varepsilon}
\newcommand{\TV}{\operatorname{TV}}
\newcommand{\HH}{\mathcal H}
\newcommand{\Lip}{\operatorname{Lip}}
\newcommand{\ceil}[1]{\left\lceil #1\right\rceil}
\newcommand{\floor}[1]{\left\lfloor #1\right\rfloor}
\newcommand{\norm}[1]{\left\lVert#1\right\rVert}
\newcommand{\sourcebox}[2]{\par\medskip\noindent\colorbox{pale}{\parbox{\dimexpr\linewidth-2\fboxsep\relax}{\small\textbf{#1.} #2}}\par\medskip}
\newif\ifincludeoriginals
\ifdefined\WITHORIGINALS\includeoriginalstrue\else\includeoriginalsfalse\fi
\begin{document}
\begin{titlepage}
\vspace*{18mm}
{\color{teal}\Large DOSSIER MATHÉMATIQUE CONSOLIDÉ\par}
\vspace{12mm}
{\Huge\bfseries Mémoire prédictive\\[5pt] certifiée\par}
\vspace{8mm}
{\LARGE Résultats, preuves et problèmes ouverts\par}
\vspace{7mm}
{\Large Programme P0--P6 et développements P4\par}
\vspace{15mm}
\rule{\linewidth}{.6pt}
\vspace{5mm}
{\large Sources mathématiques éditables en \LaTeX\\[4pt]
Douze PDF historiques conservés intégralement\\[4pt]
Relecture contradictoire par un agent IA distinct\par}
\vspace{10mm}
\sourcebox{Statut}{Les lemmes, les résultats expérimentaux et les dépendances de preuve sont distingués. L’optimalité générale et la priorité scientifique ne sont pas présentées comme validées extérieurement. Aucune médaille ou qualité de rapporteur humain n’est attribuée à l’agent de relecture.}
\vfill
{\large 3 novembre 2026\par}
\ifincludeoriginals{\small Version avec les 105 pages des douze PDF sources en annexe.\par}\else{\small Version mathématique consolidée. Les PDF sources sont dans l’archive et dans la version longue.\par}\fi
\end{titlepage}
\tableofcontents
\clearpage
\part{Synthèse et relecture critique}
'''

parts=[PREAMBLE,(ROOT/'frontmatter.tex').read_text(),(ROOT/'review_addendum.tex').read_text()]
reviewfile=ROOT/'review/referee_final.tex'
if reviewfile.exists():parts.append(reviewfile.read_text())
else:parts.append(r'\section{Rapport final de seconde lecture}'+'\n'+r'La seconde lecture est en cours de consolidation; cette version est un brouillon de travail.'+'\n')
parts.append(r'\part{Résultats et preuves du corpus}'+'\n')
parts.append((ROOT/'bilan_chapter.tex').read_text())
coverage=[]
for no,(folder,name,title,pdf,note) in enumerate(SPECS,1):
 path=ROOT/'historical_tex'/folder/name
 source=flatten_inputs(path.read_text(),path.parent)
 pre,body=source.split(r'\begin{document}',1)
 body=body.split(r'\end{document}',1)[0]
 original_body=body
 body=body.replace(r'\maketitle','').strip()
 if body.startswith(r'\begin{center}'):
  end=body.index(r'\end{center}')+len(r'\end{center}')
  body=body[end:].lstrip()
 body=body.replace(r'\begin{abstract}',r'\begin{quote}\small\textbf{Résumé original.} ').replace(r'\end{abstract}',r'\end{quote}')
 body=body.replace(r'\newpage',r'\par\medskip')
 # Preserve mathematical text; isolate all labels and bibliography keys.
 prefix=f'h{no:02d}:'
 body=re.sub(r'\\(label|ref|eqref|pageref|autoref|cite|bibitem)\{([^}]+)\}',lambda m:'\\'+m.group(1)+'{'+','.join(prefix+x.strip() for x in m.group(2).split(','))+'}',body)
 # Report paths and historical page references remain source-relative by design.
 localmacros=macros(pre)
 chapter='\\chapter{'+title+'}\n\\label{ch:hist'+str(no)+'}\n'
 chapter+='\\sourcebox{Texte historique édité}{'+note+' Les renvois à des pages et chemins dans le texte original se rapportent au document source et à son archive, sauf renvoi explicitement consolidé. Source: \\path{'+pdf+'}.}\n'
 chapter+='\\begingroup\n'+localmacros+'\n'+body+'\n\\endgroup\n'
 parts.append(chapter)
 coverage.append({'source_tex':str(path.relative_to(ROOT)),'source_pdf':pdf,'flattened_body_sha256':hashlib.sha256(original_body.encode()).hexdigest(),'body_characters':len(original_body),'math_environments_before':len(re.findall(r'\\begin\{(?:equation|align|proof|theorem|lemma|corollary)\}',original_body)),'math_environments_after':len(re.findall(r'\\begin\{(?:equation|align|proof|theorem|lemma|corollary)\}',body))})

inventory=json.loads((ROOT/'source_inventory.json').read_text())
chronology=['Bilan_Execution_Memoire_Certifiee.pdf']+[s[3] for s in SPECS]
byname={i['name']:i for i in inventory}
parts.append(r'''\appendix
\chapter{Inventaire, provenance et reproductibilité}
Les douze PDF suivants forment le corpus de cette édition. Ils sont conservés dans \path{source_pdfs/} et identifiés par des empreintes SHA256 dans \path{source_inventory.json}. Leurs pages n’ont pas été corrigées silencieusement: les compléments et réserves de la première partie priment pour l’interprétation actuelle.
\begin{longtable}{rp{.73\linewidth}r}\toprule
No & PDF source & Pages\\\midrule\endhead
''')
for i,name in enumerate(chronology,1):parts.append(f'{i} & \\path{{{name}}} & {byname[name]["pages"]}\\\\\n')
parts.append(r'''\bottomrule\end{longtable}
\section*{Contenu du projet}
\path{P4_Manuscrit_Complet.tex} contient les preuves en mathématiques éditables et compile de façon autonome. \path{historical_tex/} conserve les sources et tableaux historiques; le bilan initial est retranscrit. La version longue ajoute les douze PDF de \path{source_pdfs/}.

\path{review/} contient les rapports et contrôles nouveaux de l’agent. \path{supporting_code/} conserve les codes et résultats historiques, dont toutes les dépendances externes ne sont pas embarquées. Le codec courant utilise la bibliothèque standard Python. Le \path{README.txt} détaille la reproduction.
\begin{verbatim}
latexmk -pdf P4_Manuscrit_Complet.tex
latexmk -pdf P4_Manuscrit_Avec_12_PDF.tex
python3 review/check_local_lemmas.py
\end{verbatim}
La compilation n’exige aucun accès réseau. Elle vérifie la syntaxe et la mise en page, pas la vérité des théorèmes. Les assertions et les tests ne sont pas une formalisation dans Lean ou Coq.
''')
parts.append(r'\ifincludeoriginals'+'\n'+r'\part{Les douze PDF sources intégraux}'+'\n')
for i,name in enumerate(chronology,1):
 title=name[:-4].replace('_',' ')
 parts.append('\\includepdf[pages=-,pagecommand={\\thispagestyle{empty}},addtotoc={1,chapter,0,{Source '+str(i)+' : '+title+'},srcpdf'+str(i)+'}]{source_pdfs/'+name+'}\n')
parts.append('\\fi\n\\end{document}\n')
main='\n'.join(parts).replace('\u202f','~').replace('\u00a0','~').replace('’',"'").replace('–','--').replace('—','---').replace('‑','-')
main=re.sub(r'\{\\rm ([^{}]+)\}',lambda m: '{\\text{'+m.group(1).replace('\\ ',' ')+'}}' if any(ord(c)>127 for c in m.group(1)) else m.group(0),main)
(ROOT/'P4_Manuscrit_Complet.tex').write_text(main)
(ROOT/'P4_Manuscrit_Avec_12_PDF.tex').write_text('\\def\\WITHORIGINALS{1}\n\\input{P4_Manuscrit_Complet.tex}\n')
(ROOT/'coverage.json').write_text(json.dumps({'source_pdf_count':12,'source_pdf_pages':sum(x['pages'] for x in inventory),'historical_latex_documents':coverage,'bilan_transcribed':True},indent=2,ensure_ascii=False))
assert all(x['math_environments_before']==x['math_environments_after'] for x in coverage)
print(json.dumps({'latex_characters':len(main),'historical_sources':len(coverage),'source_pages':sum(x['pages'] for x in inventory)},indent=2))
