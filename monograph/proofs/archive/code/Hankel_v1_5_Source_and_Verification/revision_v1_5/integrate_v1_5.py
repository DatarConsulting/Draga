#!/usr/bin/env python3
"""Rebuild standalone v1.5 from its preserved v1.4 baseline."""
from pathlib import Path
HERE=Path(__file__).resolve().parent
text=(HERE/'baseline_v1_4.tex').read_text()
def replace(old,new):
    global text
    assert text.count(old)==1,(old[:100],text.count(old))
    text=text.replace(old,new)
def before(marker,name):
    replace(marker,(HERE/name).read_text()+'\n'+marker)
replace('Research manuscript v1.4','Research manuscript v1.5')
replace(r'Version 1.4\par 7 November 2026',r'Version 1.5\par 7 November 2026')
replace('''Version 1.4 adds conditional positive-basis certificates, a corrected
relative-anchor construction, selection from a supplied finite pool,
and a complete sampling theorem at a fixed rational basis.
Exact counterexamples and six open theorem targets separate these
results from general basis discovery and executable-memory claims.''','''Version 1.5 adds quantitative truncation and dyadic precision,
estimated-basis learning on a bounded-depth promise class, normalized
transfer convergence, an optimal homogeneous address executable,
and reversible coding bounds. Proved restricted results, rejected
universal claims and residual open targets are explicitly separated.''')
replace('''empirical spectral recovery remain open.
\\end{abstract}''',r'''empirical spectral recovery remain open in that version.
The version-1.5 continuation proves conditional square-root structured
truncation and explicit dyadic precision; frequent-prefix discovery and
estimated-basis learning under a bounded-depth source promise; a rank-two
parametric construction on a supplied valid interval; normalized transfer
convergence and polynomial approximation; an optimal table-free executable
on homogeneous separated similarity attractors; and exact integer-lifting
code-length bounds against specified baselines. The revised six targets
identify the remaining generality, certification and execution gaps.
\end{abstract}''')
before(r'\section{Specialized solutions beyond the obstruction}\label{sec:special}','quantitative.inc')
before(r'\section{From exact algebra to a finite executable}\label{sec:execution}','discovery.inc')
before(r'\section{Polynomial transfer computation: the remaining completeness gap}\label{sec:oracle}','memory.inc')
before(r'\section{The bridge to GDN: exact identities and an unsuccessful predictor}\label{sec:gdn}','transfer.inc')
before(r'\section{Audit of the original P0--P6 research plan}\label{sec:audit}','codec.inc')
start=text.index(r'\section{Proved scope and open theorem targets in version 1.4}')
end=text.index('\\clearpage\n\\appendix',start)
text=text[:start]+(HERE/'open.inc').read_text()+'\n\n'+text[end:]
before('\\clearpage\n\\addcontentsline{toc}{section}{References}','audit.inc')
replace('''discovery claims and exact counterexamples.
Publication novelty''',r'''discovery claims and exact counterexamples. The version-1.5 results
are proved in Sections~\ref{sec:quantitative-truncation},
\ref{sec:learned-basis}, \ref{sec:parametric-rank-two},
\ref{sec:address-memory}, \ref{sec:normalized-transfer} and
\ref{sec:integer-lifting}; Appendix~\ref{app:six-closure-audit}
maps each proposed closure to its precise established scope.
Publication novelty''')
replace('''not empirical discovery of the interval hypotheses.
Theorems''',r'''not empirical discovery of the interval hypotheses.
Theorems''') if False else None
replace('''results generates a separating pool or certifies a learned basis.
For the general target,''',r'''results alone generates a separating pool or certifies a learned basis.
Theorem~\ref{thm:bounded-pool-learning} now supplies both steps for a
bounded-depth relative-Doeblin promise class with a supplied core.
Its latent class assumptions are not inferred from finite observations.
For the general target,''')
replace('''integration with empirical spectral recovery and estimated bases
remains open.''',r'''Theorem~\ref{thm:bounded-pool-learning} adds the estimated-basis
frequency route. Integration with empirical spectral recovery and
discovery beyond its source promise remain open.''')
replace('''For the positive separated contracting class, prove that the R6 candidate
construction terminates with accepted rational certificates.''',r'''For the positive separated contracting class, prove that the R6 candidate
construction terminates with accepted rational certificates.
The normalized eigenfunction and its polynomial approximation are now
provided by Theorem~\ref{thm:normalized-transfer} and
Corollary~\ref{cor:normalized-polynomial}; downstream acceptance
and update conditions are still required.''')
replace('''The remaining objective is a computable net and compact
update with stated compilation, static storage, arithmetic width, and
runtime costs.''',r'''Theorem~\ref{thm:address-memory} now provides a compact address
executable for homogeneous separated similarity systems. The remaining
objective is an optimal effective update on broader regular classes,
with stated compilation, static storage, arithmetic width and runtime costs.''')
replace('''not the missing exact code-length identity.''',r'''not a universal exact code-length identity. The new lifting and coding
theorems provide a reversible mathematical contract with all descriptions
and model mismatch charged; useful learned-codec gains remain to be shown.''')
replace('''next target is an observable certificate that produces a suitable
relative-anchor pool, followed by perturbation control for the
estimated basis itself. The fixed-pool SPA and fixed-basis sampling
theorems isolate these two missing steps.''',r'''next target beyond Theorem~\ref{thm:bounded-pool-learning} is an
observable justification of the source promise and automatic core
discovery. Its uniform event now closes the estimated-basis perturbation
step for the supplied-core bounded-depth class.''')
replace('''Explicit
precision and bit-complexity bounds, practical implementation, and all
stored quantities checked under the declared margins.''',r'''Theorem~\ref{thm:explicit-dyadic} now provides explicit precision
under verified witnesses. Practical implementation and obtaining useful
uniform witnesses remain to be established.''')
replace('''An effective
memory--error--update theorem with all costs charged, and R6 completeness.''',r'''The homogeneous address theorem now charges all stated memory and
update costs. Broader optimal memory and full R6 completeness remain open.''')
replace('''& Arbitrary-rank acquisition and empirical SVD recovery; independently
certified interval/cone assumptions, adaptive sampling, and useful
conditional prediction and KL guarantees.''',r'''& Beyond the new supplied-core bounded-depth learning theorem:
observable source-class/core discovery, empirical SVD recovery, and
useful conditional prediction and KL guarantees.''')
replace('''as described in Appendix~\\ref{app:discovery-audit}. Earlier source
archives retain their historical versions.''',r'''as described in Appendix~\ref{app:discovery-audit}. The version-1.5
bundle includes the six-target proposal unchanged, the added proofs,
an integration patch, and new verification results described in
Appendix~\ref{app:six-closure-audit}. Earlier source archives retain
their historical versions.''')
replace(r'\end{thebibliography}',r'''\bibitem[Calderbank et~al.(1998)]{calderbank1998}
A. R. Calderbank, Ingrid Daubechies, Wim Sweldens, and Boon-Lock Yeo.
\newblock Wavelet transforms that map integers to integers.
\newblock \emph{Applied and Computational Harmonic Analysis},
5(3):332--369, 1998.
\newblock \doi{10.1006/acha.1997.0238}.
\newblock URL \url{https://cm-bell-labs.github.io/who/wim/papers/integer.pdf}.
\end{thebibliography}''')
text=text.replace('certificate. Empirical identification of those structural assumptions\nremain open.', 'certificate. Empirical identification of those structural assumptions\nremains open.')
(HERE.parent/'Hankel_Unified_Manuscript.tex').write_text(text)
print('Wrote v1.5:',len(text.splitlines()),'lines')
