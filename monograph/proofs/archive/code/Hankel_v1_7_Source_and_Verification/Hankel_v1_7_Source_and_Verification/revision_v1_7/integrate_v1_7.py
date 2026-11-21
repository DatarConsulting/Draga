"""Rebuild v1.7 from the preserved v1.6 baseline and reviewed fragments."""
from pathlib import Path
import difflib

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
baseline = (ROOT / 'baseline/Hankel_Unified_Manuscript_v1_6.tex').read_text()
text = baseline
text = text.replace('Research manuscript v1.6', 'Research manuscript v1.7')
text = text.replace('Version 1.6\\par 7 November 2026', 'Version 1.7\\par 8 November 2026')
text = text.replace('\\fancyfoot[R]{\\scriptsize 7 November 2026}',
                    '\\fancyfoot[R]{\\scriptsize 8 November 2026}')
needle = 'and an audit of the extended proposal.\n'
assert text.count(needle) == 1
text = text.replace(needle, needle +
    'Version 1.7 adds proved agent-memory compilation contracts, conditional\n'
    'probabilistic interfaces, explicit DHMC open obligations and the corrected\n'
    'rank-four perturbation norm estimate.\n')
needle = '\\end{abstract}'
assert text.count(needle) == 1
text = text.replace(needle,
    'Version 1.7 adds a self-contained DHMC theory: exact multiletter quotients,\n'
    'invariant simulation certificates, sufficient residuals, uncertain-event\n'
    'propagation, certified guards and base-mode replay. Conditional stochastic\n'
    'interfaces include occupation defects, validation, dyadic perturbation,\n'
    'rare-context bounds and finite-history forgetting. End-to-end accounting\n'
    'is proved given its premises; extraction, reader conformance, interactive\n'
    'feedback and useful benchmark gains remain explicitly unresolved.\n'
    + needle)
old = r'''The extra series $h$ satisfies
$\|H_h\|_{\rm HS}\le9/7$ by summing
$\sum_{t\ge0}(t+1)2^t3^{-2t}$. Since $\eta\le1/1024$,
$\|H_g\|<5/6+(9/7)\eta<1$.'''
new = r'''For $h(w)=3^{-|w|}+5^{-|w|}$, the triangle inequality gives
$\|H_h\|_{\rm HS}\le9/7+25/23$ by summing the two geometric
Hilbert--Schmidt series separately. The perturbation is $\eta h/2$,
and $(9/7+25/23)/2<9/7$. Since $\eta\le1/1024$,
$\|H_g\|<5/6+(9/7)\eta<1$.
The exact squared norm is $81/49+450/169+625/529$.'''
assert text.count(old) == 1
text = text.replace(old, new)
anchor = r'\section{Proved scope and residual open targets in version 1.6}'
assert text.count(anchor) == 1
math = '\n\n'.join((HERE / f).read_text() for f in
                   ('dhmc_deterministic.tex', 'dhmc_probability.tex'))
text = text.replace(anchor, math + '\n\n' + anchor)
anchor = '\\clearpage\n\\appendix'
assert text.count(anchor) == 1
text = text.replace(anchor, (HERE / 'dhmc_open.tex').read_text() + '\n\n' + anchor)
# Keep the historical audit appendices and bibliography verbatim.
anchor = r'\begin{thebibliography}'
pos = text.index(anchor)
text = text[:pos] + (HERE / 'dhmc_audit.tex').read_text() + '\n\n' + text[pos:]
(ROOT / 'Hankel_Unified_Manuscript.tex').write_text(text)
(HERE / 'v1_6_to_v1_7.patch').write_text(''.join(difflib.unified_diff(
    baseline.splitlines(keepends=True), text.splitlines(keepends=True),
    fromfile='Hankel_Unified_Manuscript_v1_6.tex', tofile='Hankel_Unified_Manuscript_v1_7.tex')))
print(ROOT / 'Hankel_Unified_Manuscript.tex')
