\section{Audit of the extended proposal and version-1.6 additions}
\label{app:extended-audit}
The subsequently supplied 52-page conversation export repeats a proposed
LaTeX extension and Python checks, claiming complete resolution of the
six version-1.4 targets. It is distinct from this manuscript's corrected
version 1.5. Some code lines are clipped at the right edge of the PDF;
the export is retained as evidence, not represented as an executable
source file. The following audit uses the visible statements and proofs.

\paragraph{What is retained.}
Lemma~\ref{lem:positive-resolvent} records the valid classical
positive-resolvent observation. Theorem~\ref{thm:reference-truncation}
repairs its application to balanced truncation and sharpens the
commuting-reference comparison. The explicit three-state example in
Section~\ref{sec:resolvent-truncation} disproves the proposed norm-of-sum
strengthening and the product-scale coupling bound. The claimed linear
rate in the export's summary is not a consequence of its displayed
square-root estimate. Independent retained-tuple stability and all
resolvent/coupling constants remain explicit assumptions.
Theorem~\ref{thm:packed-startup} retains the useful finite-window idea,
packs all startup lengths into a single integer, and extends the
homogeneous result to effective initial states without an address oracle.
No constant fixed-word bit cost is inferred from a large-integer update.

\paragraph{Estimated bases and rank-two endpoints.}
The estimated-basis argument in the export claims
$\|\widetilde B_a\|_\infty\le1$ from an observable row-sum identity.
Signed coordinates do not satisfy that implication: the matrix
$\left(\begin{smallmatrix}2&-1\\0&1\end{smallmatrix}\right)$
has row sums one and maximum absolute row sum three.
Lemma~\ref{lem:learned-basis} instead transfers through an actual
positive basis, and Theorem~\ref{thm:bounded-pool-learning} explicitly
controls same-sample selection and finite-precision margins.
Finite endpoint iteration does not in general produce an invariant
interval, even with noiseless data. For $F(t)=(1+t)/2$ started at
$t_0=0$, the upper endpoint after $j$ hull iterations is
$r_j=1-2^{-j}$, while $F(r_j)>r_j$ for every finite $j$.
Give this map letter weight $1/4$ and the constant map $G(t)=0$
letter weight $(1+t)/8$, with STOP probability $(5-t)/8$.
These data define a proper rank-two model: the reachable states $0$
and $1/2$ have different STOP probabilities. Direct substitution
into endpoint barycentric coordinates then gives a negative weight.
The clipping construction in Theorem~\ref{thm:parametric-rank-two}
and its charged closure defect are essential. Bounded-depth source
promises and a valid coordinate interval are still not learned for
arbitrary sources.

\begin{proposition}[Failure of the proposed strict cone factor]
\label{prop:false-cone-factor}
The extended proposal's factor
$\gamma'=\gamma+(1-\gamma)\xi_{\max}/(1+\xi_{\max})$
does not follow from its stated homographic hypotheses at
$L=c/(1-\gamma)$.
\end{proposition}
\begin{proof}
Take $\sigma=s_+=1$, $\xi_0=\xi_1=1/10$, and
\[
 G_0(x)=\frac{-1/2+x/5}{1+x/10},\qquad
 G_1(x)=\frac{1/2+x/5}{1+x/10},
\]
\[
 g_0(x)=\frac{1/4}{(1+x/10)^2},\qquad
 g_1(x)=\frac{3/20}{(1+x/10)^2}.
\]
On $[-1,1]$ both maps are increasing, their images are disjoint
subintervals of $[-1,1]$, and their largest derivative is
$\gamma=25/81$. The weights lie between $15/121$ and $25/81$.
For $R=2$, $w=3/2$, the complex-disc condition follows from
$|G_b(z)|\le(1/2+2/5)/(1-1/5)=9/8<w<R$.
Thus these are admissible separated homographic data.
Here $c=2/9$, $L=9/28$, and the proposed smaller constant is
$\gamma'L=331/2772<2/9$.
But $1\in\mathcal K_L$ and
\[
 (\mathcal L1)(x)=\frac{2/5}{(1+x/10)^2},\qquad
 \operatorname{Lip}(\log\mathcal L1)=2/9.
\]
Consequently $\mathcal L1\notin\mathcal K_{\gamma'L}$.
\end{proof}
The corrected choice $L>c/(1-\gamma)$ and the regularity-cone
comparison in Theorem~\ref{thm:normalized-transfer} avoid this error.
Pointwise bounds alone do not supply the projective diameter in the
regularity-cone order. Even valid normalized convergence proves none
of the unlisted downstream R6 acceptance or inversion inequalities.
Those remain in Open problem~\ref{open:r6}.

\paragraph{Lifting and the proposed entropy comparison.}
The displayed three-shear rotation formula in the export has the
opposite orientation to its left-hand rotation; splitting angles to
keep coefficients at most one also changes its stated step count.
The corrected factorization and propagated, potentially correlated
rounding errors are in Lemma~\ref{lem:rotation-lifts} and
Theorem~\ref{thm:integer-lifting}.
The proposed covariance lower bound after quantization is itself
false. Let $X$ be uniform on $[-1/10,1/10]$ and use unit-step nearest
rounding $S$. This is a log-concave source, $S=0$ almost surely, and
$X_q=S+W$ with independent $W$ uniform on $[-1/2,1/2]$ has variance
$1/12$. Nevertheless,
\[
 \left(\sqrt{\operatorname{Var}X}-\frac12\right)^2
 =\left(\frac1{\sqrt{300}}-\frac12\right)^2>\frac1{12}.
\]
Moreover, rounded-and-dithered densities need not remain log-concave:
rounding a nondegenerate Gaussian and adding within-bin uniform noise
produces a positive piecewise-constant density with jumps.
A finite concave log density on the real line must be continuous.
No lower entropy theorem for continuous log-concave densities can
therefore be applied to that construction without further justification.
Most decisively, a valid complete lossless-code upper bound cannot
lie below a valid lower bound on the same source's entropy, by
Proposition~\ref{prop:no-universal-code-gain}. A comparison to a
specified baseline remains meaningful; beating every optimal lossless
coder does not.

\paragraph{Scope of verification.}
The visible proposed Python function for truncation checks the
reference-product cancellation and a nilpotent identity, but does not
compute the advertised general truncation error. Its basis function
checks row sums without executing the claimed positivity tests. Its
rank-two example inserts an $n^{-1/2}$ perturbation by hand, rather
than deriving a sampling guarantee. The lifting example passes
rational values through floating point when rounding. These can be
illustrative checks, but they neither validate all theorem hypotheses
nor prove the six universal completion claims.
The version-1.6 companion script separately checks the new exact
counterexamples, reference Stein defects, numerical joint-Gramian
bounds, and packed startup updates. These checks supplement the
written proofs and are not a proof-assistant certification.
