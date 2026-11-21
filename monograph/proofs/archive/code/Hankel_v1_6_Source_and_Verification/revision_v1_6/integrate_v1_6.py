#!/usr/bin/env python3
"""Build standalone v1.6 from the preserved v1.5 source and new proofs."""
from pathlib import Path
HERE=Path(__file__).resolve().parent
text=(HERE/'baseline_v1_5.tex').read_text()
def replace(old,new):
    global text
    assert text.count(old)==1,(old[:100],text.count(old))
    text=text.replace(old,new)
def before(marker,file): replace(marker,(HERE/file).read_text()+'\n'+marker)
replace('Research manuscript v1.5','Research manuscript v1.6')
replace(r'Version 1.5\par 7 November 2026',r'Version 1.6\par 7 November 2026')
replace('Version 1.5 adds quantitative truncation and dyadic precision,',
        'Version 1.5 added quantitative truncation and dyadic precision,')
replace('universal claims and residual open targets are explicitly separated.',
        'universal claims and residual open targets are explicitly separated. '
        'Version 1.6 adds a commuting-reference truncation certificate, '
        'packed startup states for effective initialization, and an audit '
        'of the extended proposal.')
replace('identify the remaining generality, certification and execution gaps.\n\\end{abstract}',
        r'''identify the remaining generality, certification and execution gaps.
Version 1.6 adds a positive-resolvent certificate based on commuting
reference Gramians and removes the homogeneous address machine's
fixed-point initialization restriction by packing its startup states.
Explicit counterexamples reject the extended proposal's stronger
coupling, projective-contraction and quantization-covariance claims.
\end{abstract}''')
before(r'\section{Specialized solutions beyond the obstruction}\label{sec:special}',
       'resolvent.inc')
replace(r'''The fixed-point initialization in Theorem~\ref{thm:address-memory}
can be replaced by an effectively
supplied length-$m$ address of $x_0$ whose nested cylinder contains it,
with its acquisition cost charged. Arbitrary initial coordinates do
not automatically supply such an address.''',
        r'''Theorem~\ref{thm:packed-startup} below removes the fixed-point
initialization restriction for effective initial coordinates by encoding
startup lengths within the state index; it needs no address of $x_0$.''')
before(r'\section{Polynomial transfer computation: the remaining completeness gap}\label{sec:oracle}',
       'initial_state.inc')
replace('Proved scope and residual open targets in version 1.5',
        'Proved scope and residual open targets in version 1.6')
replace('Version 1.5 proves concrete restricted answers to all six themes.',
        'Versions 1.5--1.6 prove concrete restricted answers to all six themes.')
replace(r'''conditional square-root operator bound, and
Theorem~\ref{thm:explicit-dyadic} supplies an a priori precision.''',
        r'''conditional square-root operator bound; Theorem~\ref{thm:reference-truncation}
adds a commuting-reference certificate with explicit resolvent and
coupling constants. Theorem~\ref{thm:explicit-dyadic} supplies an a priori precision.''')
replace(r'''nonconformal systems without such a chart, and more general initial states, with the
optimal coefficient rather than a worst-contraction upper bound.''',
        r'''and nonconformal systems without such a chart, with the optimal
coefficient rather than a worst-contraction upper bound.
Theorem~\ref{thm:packed-startup} now includes arbitrary effective initial
states in the homogeneous class under its readout-evaluation assumption.''')
replace(r'''Theorem~\ref{thm:address-memory} gives
$B_f^*(2^{-p})=sp+O(1)$ with an effective table-free executable for
homogeneous separated similarities and future-observable output.
Theorem~\ref{thm:chart-memory} settles supplied effective globally
contractive charts. Extend these results to unequal contraction ratios,
and nonconformal systems without such a chart, with the optimal
coefficient rather than a worst-contraction upper bound.
Theorem~\ref{thm:packed-startup} now includes arbitrary effective initial
states in the homogeneous class under its readout-evaluation assumption.''',
        r'''Theorems~\ref{thm:address-memory} and \ref{thm:packed-startup} give
$B_f^*(2^{-p})=sp+O(1)$ with table-free homogeneous executables,
including effective initial states under the stated readout assumptions.
Theorem~\ref{thm:chart-memory} settles supplied contractive global charts.
Extend these results to unequal contraction ratios and nonconformal
systems without such a chart, preserving the optimal memory coefficient.''')
replace(r'''Appendix~\ref{app:six-closure-audit}. Earlier source archives retain
their historical versions.''',
        r'''Appendix~\ref{app:six-closure-audit}. The version-1.6 bundle adds the
extended conversation export, the commuting-reference and packed-startup
proofs, their integration patch and verification evidence, as audited in
Appendix~\ref{app:extended-audit}. Earlier source archives retain
their historical versions.''')
before('\\clearpage\n\\addcontentsline{toc}{section}{References}',
       'transcript_audit.inc')
replace(r'\end{thebibliography}',r'''\bibitem[Bourin and Lee(2019)]{bourin2019}
J.-C. Bourin and E.-Y. Lee.
\newblock On the Russo--Dye theorem for positive linear maps.
\newblock \emph{Linear Algebra and its Applications}, 2019.
\newblock \href{https://doi.org/10.1016/j.laa.2019.02.017}{doi:10.1016/j.laa.2019.02.017};
\href{https://arxiv.org/abs/1911.10573}{arXiv:1911.10573}.

\end{thebibliography}''')
(HERE.parent/'Hankel_Unified_Manuscript.tex').write_text(text)
print(f'Wrote v1.6: {len(text.splitlines())} lines.')
