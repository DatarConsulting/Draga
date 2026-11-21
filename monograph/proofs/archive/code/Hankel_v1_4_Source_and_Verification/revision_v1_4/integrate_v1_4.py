#!/usr/bin/env python3
"""Rebuild the standalone v1.4 source from its preserved v1.3 baseline."""
from pathlib import Path

HERE = Path(__file__).resolve().parent
text = (HERE / "baseline_v1_3.tex").read_text()

def replace(old, new):
    global text
    assert text.count(old) == 1, (old[:100], text.count(old))
    text = text.replace(old, new)

replace("Research manuscript v1.3", "Research manuscript v1.4")
replace(r"Version 1.3\par 7 November 2026", r"Version 1.4\par 7 November 2026")
replace("""Version 1.3 retains the earlier patch audit and adds proved rank-two
realizations, finite-prefix clipping, endpoint acquisition for both
orientations, and explicit sampling and dyadic-rounding certificates.
Counterexamples and six open theorem targets delimit their scope.""", r"""Version 1.4 adds conditional positive-basis certificates, a corrected
relative-anchor construction, selection from a supplied finite pool,
and a complete sampling theorem at a fixed rational basis.
Exact counterexamples and six open theorem targets separate these
results from general basis discovery and executable-memory claims.""")
replace("""greedy acquisition, or a general KL guarantee.
\\end{abstract}""", r"""greedy acquisition, or a general KL guarantee.
Version 1.4 adds a projectively contracting trailing-word construction,
a corrected relative-Doeblin simplex inclusion, and a classical SPA
selection guarantee with its full noise threshold. A supplied rational
positive basis yields a complete iid sampling and dyadic-output theorem,
with explicit sample thresholds and polynomial intermediate bit cost.
Pool discovery, estimated-basis certification, and integration with
empirical spectral recovery remain open.
\end{abstract}""")
replace("""arbitrary-rank acquisition program.
Publication novelty""", r"""arbitrary-rank acquisition program. Sections~\ref{sec:observable-positive}--
\ref{sec:fixed-basis-sampling} add arbitrary-rank results with supplied
cones, pools, or bases. The SPA step is explicitly an application of a
classical theorem. Appendix~\ref{app:discovery-audit} records rejected
discovery claims and exact counterexamples.
Publication novelty""")

marker = r"\section{From exact algebra to a finite executable}\label{sec:execution}"
replace(marker, (HERE / "positive_discovery.inc").read_text() + "\n" +
        (HERE / "fixed_basis_sampling.inc").read_text() + "\n" + marker)

replace("""not empirical discovery of the interval hypotheses.
For the general target,""", r"""not empirical discovery of the interval hypotheses.
Theorems~\ref{thm:trailing-word}--\ref{thm:relative-spa} add supplied-cone
approximation and relative-anchor selection from a fixed finite pool.
Theorem~\ref{thm:fixed-basis-sampling} closes direct sampling and
dyadic output at a supplied rational positive basis. None of these
results generates a separating pool or certifies a learned basis.
For the general target,""")
replace("""and the arbitrary-rank numerical bridge remain open.""", r"""remain open. Theorem~\ref{thm:fixed-basis-sampling} now closes one
arbitrary-rank numerical route at a fixed rational certified basis;
integration with empirical spectral recovery and estimated bases
remains open.""")
replace("""compiler; neither follows from the other.

For the executable-memory route,""", r"""compiler; neither follows from the other. In arbitrary rank, the
next target is an observable certificate that produces a suitable
relative-anchor pool, followed by perturbation control for the
estimated basis itself. The fixed-pool SPA and fixed-basis sampling
theorems isolate these two missing steps.

For the executable-memory route,""")
replace("""and conditioning assumptions must also remain in the claim.
The memory contribution""", r"""and conditioning assumptions must also remain in the claim.
The new relative-Doeblin and fixed-basis results extend the conditional
positive construction to arbitrary rank. They do not justify an
unconditional spectral-to-positive learning claim.
The memory contribution""")
replace("Proved scope and open theorem targets in version 1.3",
        "Proved scope and open theorem targets in version 1.4")
replace("""Empirical SVD and arbitrary-rank acquisition &
Open: interval/cone discovery, conditioning, adaptive data, and all
numerical costs in one theorem.\\\\""", r"""Supplied cones and relative-anchor pools &
Proved with explicit hypotheses: trailing-word approximation,
corrected simplex inclusion, and conditional classical SPA selection.\\
Fixed-basis arbitrary-rank sampling &
Proved for a fixed core and supplied rational certified basis;
exact rational intermediates and exact dyadic output.\\
Empirical SVD and general basis discovery &
Open: pool generation, observable separation certificates, learned
basis perturbation, and integration of all costs.\\""")
replace("[Acquisition beyond a supplied rank-two interval model]",
        "[Discovery beyond supplied interval or separable-pool certificates]")
replace("""classes, quantitative conditioning and rare-prefix costs, and whether
a specified greedy rule terminates.""", r"""classes, quantitative conditioning and rare-prefix costs, and whether
a specified greedy rule terminates. Theorem~\ref{thm:relative-spa}
selects anchors from a supplied relatively separable finite pool; it
does not generate that pool. Theorem~\ref{thm:doeblin-sandwich}
certifies a supplied true basis, not its empirical substitute.
Observable certification of these hypotheses and a complete
estimated-basis algorithm remain required.""")
replace("""not derive these errors from empirical SVD.
For rank two,""", r"""not derive these errors from empirical SVD.
Theorem~\ref{thm:fixed-basis-sampling} now proves a direct-frequency
TV rate and dyadic serialization in arbitrary rank at a fixed rational
certified basis. The remaining numerical task is to recover and
certify an estimated basis (including noisy SPA or empirical SVD),
and compose its error with the positive-cone construction. A bound
on the indices selected by SPA alone does not control this error.
For rank two,""")
replace("""new checks described in Appendix~\\ref{app:rank-two-audit}.
Two supplied copies""", r"""new checks described in Appendix~\ref{app:rank-two-audit}.
The version-1.4 bundle adds the discovery audit, exact counterexamples,
the fixed-basis sampling construction, and its reproducible checks,
as described in Appendix~\ref{app:discovery-audit}. Earlier source
archives retain their historical versions.
Two supplied copies""")

marker = "\\clearpage\n\\addcontentsline{toc}{section}{References}"
replace(marker, (HERE / "discovery_audit.inc").read_text() + "\n" + marker)
replace(r"\end{thebibliography}", r"""\bibitem[Gillis and Vavasis(2014)]{gillis2014}
Nicolas Gillis and Stephen A. Vavasis.
\newblock Fast and robust recursive algorithms for separable nonnegative
matrix factorization.
\newblock \emph{IEEE Transactions on Pattern Analysis and Machine
Intelligence}, 36(4):698--714, 2014.
\newblock \doi{10.1109/TPAMI.2013.226}.
\newblock URL \url{https://arxiv.org/abs/1208.1237}.

\bibitem[Gaubert and Qu(2015)]{gaubert2015}
St\'ephane Gaubert and Zheng Qu.
\newblock Dobrushin's ergodicity coefficient for Markov operators on cones.
\newblock \emph{Integral Equations and Operator Theory}, 2015.
\newblock \doi{10.1007/s00020-014-2193-2}.
\newblock URL \url{https://arxiv.org/abs/1307.4649}.
\end{thebibliography}""")

(HERE.parent / "Hankel_Unified_Manuscript.tex").write_text(text)
print(f"Wrote v1.4 standalone source: {len(text.splitlines())} lines")
