\clearpage
\section{Audit of the six proposed closures and version-1.5 disposition}
\label{app:six-closure-audit}
The supplied six-page proposal labels all mathematical targets
``100\% solvable.'' That designation is not a proof. This revision
proves the valid conditional replacements below, identifies false
implications, and leaves the residual obligations in Section~\ref{sec:open-targets}.

\small
\begin{longtable}{@{}>{\raggedright\arraybackslash}p{.20\textwidth}>{\raggedright\arraybackslash}p{.73\textwidth}@{}}
\toprule Target & Proved replacement and retained boundary\\\midrule
\endfirsthead
\toprule Target & Proved replacement and retained boundary\\\midrule
\endhead
1: Truncation and precision &
Theorems~\ref{thm:quantitative-truncation} and \ref{thm:explicit-dyadic}
give an explicit square-root operator bound and a priori dyadic precision
under two full Stein witnesses and separately verified truncation
stability. The proposed leading-block inheritance is false. Uniform
constants without such control remain excluded by the earlier obstruction.\\[5pt]
2: Anchor/pool discovery &
Theorem~\ref{thm:bounded-pool-learning} discovers frequent witnesses
and returns a positive dyadic law in a bounded-depth relative-Doeblin
promise class with a supplied core. Latent separation is not certified
from finite censored samples. Proposition~\ref{prop:censored-obstruction}
explains why unrestricted observable certification is impossible.\\[5pt]
3: Numerical basis and rank two &
Lemma~\ref{lem:learned-basis} and Theorem~\ref{thm:bounded-pool-learning}
close the estimated-basis perturbation gap, including adaptive selection
on a uniform event. Theorem~\ref{thm:parametric-rank-two} gives a
parametric rate on a supplied valid interval with certified map errors.
Empirical SVD and automatic recovery of all structural margins remain
separate tasks.\\[5pt]
4: R6 &
Theorem~\ref{thm:normalized-transfer} proves normalized projective
convergence using a strict regularity cone, and
Corollary~\ref{cor:normalized-polynomial} gives a terminating polynomial
eigenfunction approximator. The proposed boundary cone and diameter
claim are not used. Downstream R6 acceptance, coordinate inversion
and update invariants are not proved merely by this approximation.\\[5pt]
5: Effective memory &
Theorem~\ref{thm:address-memory} gives optimal persistent bits and a
table-free compiler on homogeneous separated similarity attractors;
Theorem~\ref{thm:chart-memory} adds supplied contractive global charts.
Bit updates cost $O(p)$; constant word-operation cost would assume a
word of $\Theta(p)$ bits. Unequal ratios and arbitrary charted classes
are not covered by the asserted optimal coefficient.\\[5pt]
6: Exact coding &
Theorems~\ref{thm:integer-lifting} and \ref{thm:actual-lifting-code}
give reversible lifting, propagated correlated-error bounds and complete
expected length for specified models/baselines. The universal entropy
saving criterion is impossible by Proposition~\ref{prop:no-universal-code-gain}.
Learned-model compression gains and peak-memory/update gains require
separate evidence.\\
\bottomrule
\end{longtable}
\normalsize

\paragraph{Precision of the audit.}
The nilpotent matrix in Section~\ref{sec:quantitative-truncation}
refutes inheritance of an arbitrary full Lyapunov witness; it is not
a counterexample to every theorem about balanced truncation.
The new proof uses the lower Stein blocks and tree resolvents and
does not rely on that inheritance. Likewise a stronger assumption
defines a restricted positive theorem; it does not retrospectively
prove a prior unconditional statement.

The proposed normalized transfer proof confuses pointwise-cone diameter
with the diameter in a log-Lipschitz cone. At $L=c/(1-\gamma)$ the
estimated regularity constant is $L$ again, which only proves
invariance. The strict $L'<L$ argument in this revision supplies a
finite, conservative diameter in the correct order. Its complexity
bound explicitly includes $n$ matrix-vector iterations and stored
coefficients. It closes the analytic approximation obligation, not
all operations of a downstream encoder.

In the codec proposal, a valid entropy upper bound below a valid
entropy lower bound would contradict lossless coding itself. The
correct comparison is to a specified restricted coder, with model
mismatch and descriptions charged. A log-concavity assumption does
not remove that contradiction. Integer lifting is not presented as
a new transform invented in this manuscript.

\paragraph{Verification and provenance.}
The companion \texttt{verify\_closures\_v1\_5.py} and its JSON record
finite exact algebraic and numerical checks of the new formulas,
counterexamples and constructions. These checks supplement the
proofs; they are not a proof-assistant formalization, an empirical
SVD study, or a new held-out neural-codec benchmark. The original
six-page proposal is retained unchanged in the source bundle.
