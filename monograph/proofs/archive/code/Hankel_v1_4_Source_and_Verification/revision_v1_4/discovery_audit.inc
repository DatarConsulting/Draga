\clearpage
\section{Higher-rank discovery audit and version-1.4 scope}
\label{app:discovery-audit}
The successive discovery proposals supplied after version 1.3 are
audited together here. This version retains conditional results with
proofs, gives exact counterexamples to invalid implications, and
separates a fixed-basis numerical algorithm from discovery of that
basis. It does not adopt the supplied patches as a complete solution
of Priority 2.

\small
\begin{longtable}{@{}>{\raggedright\arraybackslash}p{.25\textwidth}>{\raggedright\arraybackslash}p{.68\textwidth}@{}}
\toprule Proposed claim & Disposition in version 1.4\\\midrule
\endfirsthead
\toprule Proposed claim & Disposition in version 1.4\\\midrule
\endhead
Bounded-depth observable normalization &
Exact identities and a finite positivity certificate are proved in
Theorem~\ref{thm:observable-certificate}. Finding the basis is separate.\\[5pt]
Trailing-word contraction in arbitrary rank &
Theorem~\ref{thm:trailing-word} proves a supplied-cone approximation
with at most $N_d(L)$ states. Cone order replaces an unjustified
coordinate $\ell^1$ inequality. No rank-$R$ state bound follows.\\[5pt]
Perron-ray purification and contraction-assisted SPA &
False without additional separation hypotheses. The counterexamples
below distinguish full Hankel rank, projective contraction, and
availability of different nearly pure anchors.\\[5pt]
Unrestricted dilation from row-relative Doeblin margins &
False. Rare-letter row sums can be mixed with negative coefficients.
The relative-inner-simplex construction instead requires both
inclusions and the corrected column-sum factor.\\[5pt]
Relative Doeblin simplex and positive margins &
Proved in Theorem~\ref{thm:doeblin-sandwich}. The interval for $m$
must be nonempty; strict positivity additionally needs positive
letter masses. The actual anchor matrix is supplied.\\[5pt]
Latest patch's relative SPA constants &
Replaced by Theorem~\ref{thm:relative-spa}. The cited classical theorem
requires $\min\{1/(2\sqrt{R-1}),1/4\}$, not $1/2$ alone.
The nearly pure witness displacement needs $\sqrt2\eta_0\|G\|_2$,
and observed selected rows must be distinguished from true rows.\\[5pt]
Joint counts and floating-point normalization &
The per-pair concentration argument safely gives the $R^2b_n$ bound
in \eqref{eq:joint-count-event}. The patch's interchange of a letter
sum and row maximum is not justified. Rounded divisions and STOP
remainders need not sum exactly to one.\\[5pt]
Universal constant for partially pivoted LU &
Not retained: elimination growth, input/basis representation, and
rounding accumulation are not bounded by the asserted absolute
constant. Theorem~\ref{thm:fixed-basis-sampling} supplies an exact
rational intermediate algorithm and exact dyadic output instead.\\[5pt]
Sample/bit theorem at a certified basis &
Closed for the fixed rational basis and fixed core of
Theorem~\ref{thm:fixed-basis-sampling}, with explicit sufficient
sample thresholds and polynomial bit complexity. Estimated-basis
discovery is not included.\\[5pt]
General metastable KL jump from repeated cuts &
False as stated: whole-future conditional divergences cannot be
added over overlapping prediction epochs. The example below
disproves the claimed bound; Theorem~\ref{thm:kl-fixed-states}
remains the proved fixed-state discontinuity result.\\
\bottomrule
\end{longtable}
\normalsize

\subsection{Exact obstructions to the proposed acquisition arguments}
\begin{proposition}[Uniform Doeblin positivity prevents arbitrary hidden-state purification]
\label{prop:doeblin-purity-obstruction}
If $A_{a,ij}\ge\mu(A_a\mathbf1)_i$ for every letter, each
nonempty-prefix posterior has all coordinates at least $\mu$.
The same is true of a normalized positive left Perron vector of a
nonempty word operator. Thus a collection of such anchors with
diagonal entries at least $1-\eta$ requires $\eta\ge(R-1)\mu$.
\end{proposition}
\begin{proof}
For a posterior $\pi$, $(\pi A_a)_j\ge\mu\pi A_a\mathbf1$.
Normalize and apply this observation to the last letter. A Perron
state returns to itself after the word and obeys the same inequality.
The other $R-1$ coordinates sum to at least $(R-1)\mu$.
\end{proof}
In particular $\eta<\mu s_*/(4R)$ is incompatible with these anchors
for $R\ge2$. An empty-prefix anchor supplies at most one exception.
The relative coordinates in Theorem~\ref{thm:doeblin-sandwich} avoid
this contradiction by using the inner simplex $J_\mu$.

\begin{proposition}[Noncommuting rank three with only one Perron ray]
\label{prop:common-perron}
Let $\alpha=(1,0,0)$ and
\[
 A_0=\frac1{30}\begin{pmatrix}4&1&1\\1&3&1\\1&2&4\end{pmatrix},
 \qquad
 A_1=\frac1{30}\begin{pmatrix}2&3&1\\3&1&2\\1&2&3\end{pmatrix}.
\]
With $\tau=(3/5,19/30,17/30)^\top$, this is a proper rank-three
law with noncommuting strictly positive transitions. Every nonempty
word operator has the same normalized left Perron vector
$(1/3,1/3,1/3)$.
\end{proposition}
\begin{proof}
The rows normalize and all stopping probabilities are positive.
For the core $\{\varepsilon,0,1\}$ on both sides, direct rational
calculation gives $\det H=-11/50625000\ne0$. Hence the rank is three.
The $(1,2)$ entry of $A_0A_1-A_1A_0$ is $1/450$.
Both matrices have column sums $1/5$, so their products have constant
column sum $5^{-|w|}$. Strict positivity makes the uniform left
Perron vector unique for each nonempty word.
\end{proof}

\paragraph{Dilation without a valid outer inclusion.}
Take $R=2$, $\eta=1/100$, $\mu=1/4$, $s_*=1/2$,
\[
 A_0=\begin{pmatrix}1/10000&1/10000\\1/8&1/8\end{pmatrix},\quad
 A_1=\tfrac14\mathbf1\mathbf1^\top-A_0,\quad
 \Pi=\begin{pmatrix}1&0\\1/100&99/100\end{pmatrix}.
\]
Although $\eta<\mu s_*/(4R)=1/64$, setting $m=1/98$ and
$W=D_m\Pi$ gives
\[
 W=\frac1{9800}\begin{pmatrix}9899&-99\\-1&9801\end{pmatrix},
 \qquad (WA_0W^{-1}\mathbf1)_1=-\frac{113851}{49000000}<0.
\]
This disproves the claimed algebraic positivity implication even
when a nearly diagonal stochastic anchor matrix is supplied.

\paragraph{The missing column-sum factor is substantive.}
With $R=3$, $\mu=1/5$, $\eta=1/10$, $m=49/100$, and
\[
 \Pi=\begin{pmatrix}1&0&0\\1/10&9/10&0\\0&0&1\end{pmatrix},
 \qquad
 W=D_m\Pi J_\mu
   =\frac1{2500}\begin{pmatrix}2431&59&10\\208&2282&10\\-39&59&2480\end{pmatrix},
\]
the old interval $1/8<m<1/2$ holds but $W_{31}<0$.
The corrected upper endpoint is $5/11$. This is realizable with
actual one-letter anchors, not only abstract matrices: put $Z=\Pi J_\mu$,
let the columns of
\[
 Q=\begin{pmatrix}1/4&1/1000&1/1000\\
                  1/10000&1/4&1/1000\\
                  1/10000&1/1000&1/4\end{pmatrix}
\]
be $w_a$, and define $A_a=w_aZ_{a,:}$, $a=0,1,2$, with uniform
initial distribution and STOP as the remainder. Every one-letter
posterior is the corresponding row of $Z$, all Doeblin inequalities
hold, and $\tau=(187/250,7489/10000,7489/10000)^\top$.
The one-letter prefix/suffix core has rank three. Nevertheless
$(WA_0W^{-1}\mathbf1)_3=-94961/25000000$.
Choosing $m=1/4$ in the corrected interval makes all transformed
coefficients strictly positive for this same source.

\subsection{Verified interior two-state example}
For the earlier source
$A_0=\frac1{20}\left(\begin{smallmatrix}3&2\\3&3\end{smallmatrix}\right)$,
$A_1=\frac1{20}\left(\begin{smallmatrix}2&2\\2&3\end{smallmatrix}\right)$,
$\alpha=(1/10,9/10)$, and $\tau=(11/20,9/20)^\top$, take
\[
 U=\begin{pmatrix}30/59&29/59\\20/49&29/49\end{pmatrix},\quad
 D(4,1)=\begin{pmatrix}2&-1\\-4&5\end{pmatrix},\quad
 V=\frac1{2891}\begin{pmatrix}1760&1131\\20&2871\end{pmatrix}.
\]
Then
\[
 \alpha_V=(897/5800,4903/5800),\qquad
 \tau_V=(29539/57820,26059/57820)^\top,
\]
\[
 A_{V,0}=\frac1{100606800}
       \begin{pmatrix}24761923&2357717\\24727123&5420117\end{pmatrix},
 \quad
 A_{V,1}=\frac1{100606800}
       \begin{pmatrix}16461862&5627438\\16427062&8689838\end{pmatrix}.
\]
These exact matrices certify $\gamma=1/50$ and $s=9/20$.
The various larger or differently normalized matrices in intermediate
patches are not used. For the separate $R=2$ sandwich example with
$\mu=1/4$, $\Pi=\left(\begin{smallmatrix}1&0\\1/100&99/100\end{smallmatrix}\right)$,
and $m=1/98$, the correct values are
\[
 W=\frac1{19600}\begin{pmatrix}14799&4801\\4899&14701\end{pmatrix},
\]
\[
 WA_0W^{-1}=\frac1{970200000000}
       \begin{pmatrix}29484656149&30074228951\\
                       90086011249&91887363851\end{pmatrix}>0,
\]
where $A_0$ is the rare-letter matrix in the preceding obstruction.

\subsection{Why repeated whole-future KL cuts cannot be added}
\begin{proposition}[Counterexample to the proposed occupation-weighted KL bound]
\label{prop:repeated-kl-cuts}
The proposed lower bound $W_0\log(K/N)$ from $K>N$ disjoint future
signatures and cumulative prefix occupations is false without
nonoverlapping information increments.
\end{proposition}
\begin{proof}
On $\{0,a,b,c\}$ let
\[
 p=\tfrac13(\delta_{a0^{100}a}+\delta_{b0^{100}b}+\delta_{c0^{100}c}).
\]
For each $k\in\{a,b,c\}$ use prefixes $k0^j$, $0\le j\le100$,
and signature event $E_k=\{0^mk:m\ge0\}$. The 101 constant-length
cuts are prefix-free, the signatures are disjoint, and their
conditional probabilities equal one. Thus the proposed $N=2$, $K=3$
bound is $101\log(3/2)>40.95$. The one-state law with STOP $1/103$,
letter-$0$ probability $100/103$, and each other letter probability
$2/309$ is also representable by two identical states. Its divergence is
\[
 3\log103+100\log(103/100)+\log(3/4)<16.58.
\]
This contradicts the proposed bound. Each cut separately admits a
conditional KL lower bound, but its full future contains information
already present at the other cuts.
\end{proof}
If there are $M$ such cuts, summing their separate valid inequalities
and dividing by $M$ is legitimate. Adding them without this factor
requires signatures measurable in successive disjoint prediction
blocks and a corresponding chain rule. Neither version establishes
a general spectral KL discontinuity or directed conductance theorem.

\subsection{Reproduction and remaining obligations}
The companion \texttt{verify\_discovery\_v1\_4.py} checks the exact
matrices, the noncommuting rank witness, the counterexamples, the
corrected simplex inclusions, the probability normalization after
dyadic output, and deterministic perturbation estimates. Its JSON
records the executed cases. Tests do not replace the general proofs,
and no new LLM or codec benchmark is claimed.

Version 1.4 closes supplied-cone finite-prefix approximation,
conditional relative-anchor selection from a fixed finite pool,
the corrected simplex realization, and sampling/serialization at a
supplied rational positive basis. It leaves pool discovery, observable
certification of latent separation hypotheses, and integration of
estimated bases with the spectral--conic bridge as explicit open
obligations. A nonempty semialgebraic feasibility problem can be solved
by fixed-dimensional methods \citep{basu2006,renegar1992}; fixing
$V$ in $VW=I$ merely determines $W$, and is not a basis-discovery LP.
Uniformly polynomial dependence on varying $R$ is not established.
