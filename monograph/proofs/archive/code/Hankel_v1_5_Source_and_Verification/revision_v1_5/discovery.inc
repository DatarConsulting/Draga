\subsection{Estimated bases and discovery in a bounded-depth promise class}
\label{sec:learned-basis}
The following perturbation argument closes the estimated-basis gap
for a stated finite-dimensional source class. All matrix norms in
this subsection are maximum absolute row sums.

\begin{lemma}[A perturbed normalized basis preserves a strict positive realization]
\label{lem:learned-basis}
Suppose $V\mathbf1=\mathbf1$ realizes a proper law with
$\alpha_j,A_{a,ij}\ge\gamma>0$ and $\tau_i\ge s>0$.
Let $\widehat V=V+E$, $E\mathbf1=0$, and
$\rho=\|E\|_\infty\|V^{-1}\|_\infty<1$.
In the new basis the exact coefficients of the \emph{same} law satisfy
\begin{align*}
 \|\alpha'-\alpha\|_1&\le\rho/(1-\rho),\\
 \sum_a\|A'_a-A_a\|_\infty&\le 2R\rho/(1-\rho),\\
 \|\tau'-\tau\|_\infty&\le\rho.
\end{align*}
Thus $\rho\le\min\{1/2,\gamma/(8R),s/2\}$ suffices for
$\alpha'_j,A'_{a,ij}\ge\gamma/2$ and $\tau'_i\ge s/2$.
Also $\|\widehat V^{-1}\|_\infty\le2\|V^{-1}\|_\infty$.
\end{lemma}
\begin{proof}
Set $K=EV^{-1}$, so $\widehat V=(I+K)V$, $K\mathbf1=0$.
Then $\alpha'=\alpha(I+K)^{-1}$,
$A'_a=(I+K)A_a(I+K)^{-1}$ and $\tau'=(I+K)\tau$.
The Neumann bound gives the first and third estimates, and
\[
 A'_a-A_a=(KA_a-A_aK)(I+K)^{-1}.
\]
Nonnegativity and normalization give
$\sum_a\|A_a\|_\infty\le R$, proving the second estimate.
The row identities remain exact; the margins imply positivity and
termination. The inverse bound is the same Neumann estimate.
\end{proof}

For an estimated observable row $z$ whose true row $f$ satisfies
$f\mathbf1=1$, use the affine normalization
\[
 \mathcal N(z)=z+(1-z\mathbf1)e_1^\top.
\]
Then $\mathcal N(z)\mathbf1=1$ exactly and
$\|\mathcal N(z)-f\|_1\le2\|z-f\|_1$. This works for signed
observable coordinates and avoids an unjustified probabilistic
interpretation of those coordinates.

\begin{theorem}[Frequent-prefix discovery, estimated basis, and dyadic law]
\label{thm:bounded-pool-learning}
Fix $R\ge2$, $d$, $L_0$, a supplied core of depth at most $L_0$,
and positive rational margins. Consider the following \emph{promise
class}: the source has a minimal proper $R$-state realization with
the Doeblin hypotheses of Theorem~\ref{thm:doeblin-sandwich},
stop margin $s_*$, $w_{\min}>0$, and $0<\mu<1/R$;
$\sigma_R(H)\ge\sigma_0$, core prefix masses are at least $q_*$;
for every relative vertex there is a witness prefix of depth at most
$L_0$, mass at least $p_*$, and relative purity at least $1-\eta_0$.
The realization and witness prefixes need not be supplied.

Choose $0\le\eta<1/2$ satisfying the strict inequality
\eqref{eq:sandwich-feasibility}, and a rational $m\in(m_-,m_+)$
with $C=1+(R-1)\eta$. Set
\begin{align*}
 a_\mu&=1-R\mu,& h_0&=\sqrt R/\sigma_0,& D_0&=1+2(R-1)m,\\
 s_G&=a_\mu/\sqrt R,& S_G&=\sqrt R h_0,
 & C_0&=1+80(S_G/s_G)^2,\\
 v_-&=\frac{1+(R-2)\mu}{a_\mu(1-2\eta)},&
 v_+&=D_0h_0,&
 \gamma_0&=\frac{\min(m-m_-,m_+-m)}{1+Rm}\min(1,w_{\min}).
\end{align*}
Observe $n$ iid words censored at $2L_0+1$, including termination
flags. Put $M_0=N_d(L_0)$,
\[
 K=M_0(R+1)+R^2+2R,\qquad
 t=\sqrt{\log(4K/\delta)/(2n)},\qquad
 b_n=\sqrt{d/n}+\sqrt{2\log(2R^2/\delta)/n}.
\]
Use the frequent observed-prefix pool
$\widehat{\mathcal P}=\{w:|w|\le L_0,\widehat q_w\ge3p_*/4\}$.
Compute its rows by Proposition~\ref{prop:pool-estimation}, then apply
$\mathcal N$. The following explicit sufficient tests are imposed:
\begin{align*}
 Rh_0t&\le1/2,& t&\le\min(p_*/4,q_*/2),\\
 e_0&=16Rh_0t/p_*+(2Rh_0^2+h_0)t,& e_F&=2e_0,\\
 e&=e_F+\sqrt2\eta_0S_G,&
 e&<s_G\theta_R/C_0,\qquad
 (e_F+C_0e)/s_G\le\eta,\\
 e_V&=D_0e_F,&
 e_Vv_-&\le\min\{1/2,\gamma_0/(8R),s_*/2\}.
\end{align*}
Here $\theta_R$ is as in Theorem~\ref{thm:relative-spa}.
Run SPA, let $\widehat U$ be its selected normalized rows, and set
$\widehat V=D_m\widehat U$. For the empirical realization define
\begin{align*}
 \bar v_-&=2v_-,&\bar v_+&=v_++e_V,\\
 e_\alpha&=\bar v_-[(2R+1)h_0+2Rh_0^2]t,\\
 K_0&=2Rh_0/q_*^2+(4R^2h_0^2+4Rh_0)/q_*,\\
 e_A&=\bar v_+\bar v_-
       [2h_0R^2b_n/q_*+K_0t].
\end{align*}
Require $e_\alpha\le\gamma_0/8$,
$e_A\le\min(\gamma_0,s_*)/4$, and $h_b=2^{-b}\le\gamma_0/8$.
Use exact rational arithmetic, normalize the initial vector, verify
initial/transition entries at least $\gamma_0/4$ and STOP at least
$s_*/4$, and round downward as in Theorem~\ref{thm:fixed-basis-sampling}.
On a failed check return \texttt{UNRESOLVED}.
With probability at least $1-\delta$ all checks pass and the returned
proper law on $R$ states satisfies
\begin{equation}
 \TV(p,q_b)\le\frac{e_\alpha}{1-e_\alpha}+(R-1)h_b
                    +\frac{2(e_A+dRh_b)}{s_*}.
 \label{eq:learned-basis-tv}
\end{equation}
For fixed structural margins and sufficiently small $\eta_0$, the rate is
\[
 O\!\left(\sqrt{\frac{d+L_0\log d+\log(R/\delta)}n}+d2^{-b}\right).
\]
Counting uses a trie; the observed pool has at most $1+nL_0$ members.
With this supplied core the algorithm has polynomial bit cost in
$n,L_0,d,R,b$ and the rational parameter bit lengths.
\end{theorem}
\begin{proof}
Apply Hoeffding simultaneously to every prefix of depth at most
$L_0$ and its $R$ completed-word suffix queries, and to the core
queries. The union bound uses $K$ even though only observed frequent
prefixes are processed. Combine it with the joint-letter event in
\eqref{eq:joint-count-event}; the two failure budgets are $\delta/2$.
Every witness is retained, and every retained prefix has true mass
at least $p_*/2$. Proposition~\ref{prop:pool-estimation}, with this
lower mass, and affine normalization give row $\ell^1$ error $e_F$.

In native coordinates $\widetilde P=D_q^{-1}P$ is stochastic.
Thus $\sigma_R(J_\mu\widetilde P^{-1})\ge a_\mu/\sqrt R=s_G$.
The suffix matrix $C$ in $H=PC$ has nonnegative rows of mass at most
one, because suffixes are distinct. Hence
$\|V_*\|_\infty=\|CH^{-1}D_q\|_\infty\le h_0$ and
$\|J_\mu V_*\|_2\le S_G$.
The conservative SPA tests therefore imply that the selected
\emph{true} rows have relative purity $1-\eta$.

Their basis $V=D_mU$ has margins $\gamma_0,s_*$ by
Theorem~\ref{thm:doeblin-sandwich}. Every true observable residual
row has $\ell^1$ norm at most $h_0$, giving $\|V\|_\infty\le v_+$.
The identities
\[
 V^{-1}=\widetilde P J_\mu^{-1}\Pi^{-1}D_m^{-1},\quad
 \|\Pi^{-1}\|_\infty\le(1-2\eta)^{-1},\quad
 \|D_m^{-1}\|_\infty=1
\]
give the claimed $v_-$ bound. Moreover
$\|\widehat V-V\|_\infty\le e_V$ and both bases send
$\mathbf1$ to $\mathbf1$. Lemma~\ref{lem:learned-basis} implies
that the true law in the empirical basis has margins
$\gamma_0/2,s_*/2$ and the displayed inverse/norm bounds.

Now apply the deterministic product estimates from the proof of
Theorem~\ref{thm:fixed-basis-sampling} on the same uniform event.
They hold for every selected basis obeying these norm bounds; no
sample-independence claim for $\widehat V$ is needed. The raw initial
and transition errors are bounded by $e_\alpha,e_A$. The small-error
tests imply all stipulated margins after initial normalization.
Downward rounding preserves properness. Pathwise perturbation with
true STOP at least $s_*/2$ proves \eqref{eq:learned-basis-tv}.
A change of exact basis itself changes no word probability; it is
used to prove positivity, not counted as an additional distribution error.

All retained prefixes occur in the observed trie. Query lookup,
matrix products and SPA require polynomially many operations in its
size and the stated parameters. Squared-norm SPA comparisons can be
performed exactly on rationals; no square roots are needed to execute
the pivots. Fraction-free algebra and dyadic floor operations give
the polynomial bit bound as in Theorem~\ref{thm:fixed-basis-sampling}.
\end{proof}
This closes frequent-prefix discovery and the numerical basis bridge
\emph{within the stated promise class and with a supplied core}.
It does not certify the latent class assumptions from finite data,
discover the core in uniformly polynomial time, or use empirical SVD.
The finite exact positivity checks certify the output automaton;
the approximation guarantee additionally uses the source-class promise.

\begin{proposition}[Bounded-depth observations cannot certify an unrestricted source]
\label{prop:censored-obstruction}
For every censoring depth $D$, the laws $\delta_{a^{D+1}}$ and
$\delta_{a^{D+2}}$ induce identical censored observations, including
termination flags through $D$, but their TV distance is one.
For $\epsilon<1/2$ and $\delta<1/2$, no procedure using only these
observations can approximate both laws within $\epsilon$ with
probability at least $1-\delta$.
\end{proposition}
\begin{proof}
Every observation is the same nonterminated $D$-letter string.
The two TV balls of radius $\epsilon$ are disjoint, by the triangle
inequality. An output distribution that is identical under the two
sources cannot give both disjoint events probability greater than
$1/2$. Arbitrarily many such censored observations do not help.
\end{proof}

\begin{proposition}[Fixed-block empirical singular-space certificate]
\label{prop:fixed-block-svd}
Let $K$ be an $M\times N$ block of completed-word probabilities,
with its query words fixed before sampling, $\rank K=r$, and
$\sigma_r(K)\ge s_0>0$. From $n$ iid observations let $\widehat K$
be the frequency matrix and put
\[
 e=\sqrt{MN}\sqrt{\log(2MN/\delta)/(2n)}.
\]
With probability at least $1-\delta$,
$\|\widehat K-K\|_2\le e$. If $e<s_0/2$, the leading $r$ left
and right singular projectors obey
\[
 \|\widehat P-P\|_2,\ \|\widehat Q-Q\|_2
          \le e/(s_0-e),\qquad
 \|\widehat K_r-K\|_2\le2e.
\]
If $3e<s_0$, thresholding empirical singular values at $2e$ selects
exactly $r$. The population rank and gap remain assumptions.
\end{proposition}
\begin{proof}
Hoeffding and a union bound give entry error at most
$e/\sqrt{MN}$ and thus the operator bound via the Frobenius norm.
The variational singular-value inequalities give
$\widehat\sigma_r\ge s_0-e$ and
$\widehat\sigma_{r+1}\le e$.
Writing $\widehat K\widehat V_r=\widehat U_r\widehat\Sigma_r$
and multiplying by $I-P$ gives
$(I-P)\widehat U_r=(I-P)(\widehat K-K)\widehat V_r
\widehat\Sigma_r^{-1}$. The norm is at most $e/(s_0-e)$.
For two projectors of equal rank this is their operator-norm distance;
the transposed argument gives the right projector. Finally
$\|\widehat K_r-K\|\le\|\widehat K_r-\widehat K\|+e\le2e$.
\end{proof}
This controls a fixed finite block. It neither estimates an infinite
Hankel tail nor supplies a positive cone or justifies adaptive block
selection outside a declared uniform family. Validated numerical
SVD errors must be added when the computed decomposition is inexact.

\subsection{A bounded-coordinate rank-two parametric construction}
\label{sec:parametric-rank-two}
\begin{theorem}[Estimated maps and endpoint clipping on a supplied valid interval]
\label{thm:parametric-rank-two}
Suppose $f_t$, $t\in I$, is a supplied affine family of probability
laws on a compact interval, with $p=f_{t_0}$,
$\|f_t-f_u\|_1\le K_f|t-u|$, true STOP at least $s>0$, and
letter masses $w_a(t)$ and projective updates $F_a(t)\in I$.
Let estimated affine-fractional maps $\widehat F_a:I\to I$, with
no pole on $I$, be contractions of modulus
at most $\bar\gamma<1$, and let nonnegative affine estimated masses
$\widehat w_a$ satisfy uniformly
\[
 \max_a|\widehat F_a(t)-F_a(t)|\le e_F,
 \qquad \sum_a|\widehat w_a(t)-w_a(t)|\le\zeta<s.
\]
Starting from $[t_0,t_0]$, perform $J$ all-letter endpoint-hull
iterations with $\widehat F_a$. At the two resulting endpoints,
use masses $\widehat w_a$, STOP as the remainder, and route to the
barycentric representation of $\widehat F_a(t)$ clipped to the
current interval. If the interval is degenerate, use one state.
Round probabilities down with step $h_b$ and initial normalization
by a remainder. The proper law on at most two states satisfies
\begin{equation}
 \TV(p,q_{J,b})\le
 \frac{\zeta+\frac12K_f(e_F+\bar\gamma^J\operatorname{diam}I)}{s-\zeta}
       +h_b+\frac{2dh_b}{s-\zeta}.
 \label{eq:parametric-rank-two}
\end{equation}
Consequently fixed-depth coordinate estimates with
$e_F,\zeta=O(n^{-1/2})$ and $J=O(\log n)$ give a parametric
$O(n^{-1/2}+d2^{-b})$ TV rate. The interval and its validity are inputs.
\end{theorem}
\begin{proof}
The endpoint operator, taking minimum and maximum over $t_0$ and
all endpoint images, is $\bar\gamma$-Lipschitz in the maximum norm.
For an affine-fractional map without a pole on $I$, extrema on an
interval occur at its endpoints. Thus the interval expands monotonically and the
distance from any next image to the current interval is at most
$\bar\gamma^J\operatorname{diam}I$.
Initial mixing represents $f_{t_0}$ exactly. In a row, the STOP/letter
mass discrepancies contribute at most $2\zeta$ in $\ell^1$.
Since $\sum_a\widehat w_a\le1-s+\zeta\le1$, replacing true successors
by clipped estimated successors contributes at most
$K_f(e_F+\bar\gamma^J\operatorname{diam}I)$.
Every constructed STOP is at least $s-\zeta$. The occupation bound
and Theorem~\ref{thm:flows} give the first fraction. Initial dyadic TV
is at most $h_b$; per-row dyadic TV is at most $2dh_b$.
Apply the stopped-coupling bound for the remaining terms.
\end{proof}
For fixed-depth observable matrices, a certified nonzero denominator
on $I$ makes their fractional maps locally Lipschitz in the entries;
the inverse estimates already proved in Section~\ref{sec:fixed-basis-sampling}
then give the stated parametric errors. This does not assert that
the unknown valid interval, denominator margin or contraction bound
can always be inferred from bounded-depth data. Exact limiting rays
need not be evaluated: finite endpoint iteration and clipping suffice.
