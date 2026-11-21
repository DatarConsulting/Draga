\subsection{Observable coordinates and a finite positive certificate}
\label{sec:observable-positive}
The results below distinguish a supplied positive basis, selection from a
supplied candidate pool, and discovery of such a pool. Throughout this
continuation, matrix $\|\cdot\|_\infty$ is the maximum absolute row sum,
$\|\cdot\|_{\max}$ is the maximum absolute entry, and vector norms have
their usual meanings. All positive realizations are row-normalized:
$\alpha\mathbf1=1$ and $\sum_a A_a\mathbf1+\tau=\mathbf1$.

\begin{theorem}[Observable normalization and finite positivity certificate]
\label{thm:observable-certificate}
Let $p(w)=\alpha_*A_{*,w}\tau_*$ have a minimal $R$-state positive
realization with $\tau_*\ge s_*\mathbf1>0$. Fix distinct core prefixes
$u_i$ and suffixes $v_j$, $1\le i,j\le R$, for which
\[
 H_{ij}=p(u_iv_j),\quad (M_a)_{ij}=p(u_iav_j),\quad
 h_j=p(v_j),\quad t_i=p(u_i)
\]
give an invertible $H$. Put $B_a=M_aH^{-1}$, $\beta=hH^{-1}$,
$q_i=p(u_i\Sigma^*)$, and $D_q=\diag(q_i)$. Then $q_i>0$ and
\[
 q=(I-\textstyle\sum_aB_a)^{-1}t,\qquad
 \widetilde B_a=D_q^{-1}B_aD_q,\quad
 \widetilde\beta=\beta D_q,\quad \widetilde t=D_q^{-1}t
\]
satisfy $\widetilde\beta\mathbf1=1$ and
$\sum_a\widetilde B_a\mathbf1+\widetilde t=\mathbf1$.
For any invertible $V$ with $V\mathbf1=\mathbf1$, define
\[
 \alpha=\widetilde\beta V^{-1},\qquad
 A_a=V\widetilde B_aV^{-1},\qquad \tau=V\widetilde t.
\]
If $\alpha\ge0$, every $A_a\ge0$, and $\tau\ge s\mathbf1$ for some
$s>0$, these finitely many inequalities certify an exact proper
$R$-state realization of $p$.
\end{theorem}
\begin{proof}
Let $P_{i,:}=\alpha_*A_{*,u_i}$ and let column $j$ of $C$ be
$A_{*,v_j}\tau_*$. Since $H=PC$ is invertible, both factors are
invertible. Thus $B_a=PA_{*,a}P^{-1}$, $\beta=\alpha_*P^{-1}$,
and $t=P\tau_*$. Row normalization and the stopping margin give
$(I-\sum_aA_{*,a})^{-1}\tau_*=\mathbf1$, hence $q=P\mathbf1$.
Every row of $P$ is nonzero and nonnegative, proving $q_i>0$.
The normalized matrix $\widetilde P=D_q^{-1}P$ is row-stochastic;
the displayed identities follow by similarity. A second similarity
by $V$ preserves every coefficient. Moreover $V^{-1}\mathbf1=\mathbf1$,
so the resulting initial vector and transition/stop rows sum to one.
The nonnegative aggregate has row sums at most $1-s$; hence it is
transient, every state terminates, and the word law has mass one.
\end{proof}
The certificate eliminates an infinite family of word-positivity tests
\emph{after} $V$ is obtained. It does not find $V$. Observable coordinates
can be signed; they are not hidden-state posterior probabilities.

\subsection{Contracting cones and finite trailing-word dictionaries}
\label{sec:trailing-word}
\begin{theorem}[Conditional finite-prefix construction under projective contraction]
\label{thm:trailing-word}
Suppose a supplied finite-dimensional proper cone $\mathcal C$ carries
positive linear letter maps and a linear map $x\mapsto f_x$ to
nonnegative summable word functions. Let $\ell(x)=\sum_w f_x(w)$ be
strictly positive on $\mathcal C\setminus\{0\}$; normalized states have
$\ell(x)=1$. Assume letter shifts and normalized updates agree with this
representation, the initial law is $p=f_{x_0}$, and every normalized
state has stopping probability at least $s>0$. Each nonzero letter map
sends $\mathcal C\setminus\{0\}$ into its interior, has Hilbert diameter
at most $\Delta<\infty$, and has projective contraction coefficient at
most $\gamma\in(0,1)$. Then for every $L\ge1$ there is a proper automaton
$q_L$, using at most $N_d(L)=\sum_{j=0}^Ld^j$ actual-prefix residuals, with
\begin{equation}
 \TV(p,q_L)\le\frac{1-s}{s}
       \bigl(e^{\Delta\gamma^{L-1}}-1\bigr).
 \label{eq:trailing-tv}
\end{equation}
The same bound holds after keeping one representative of each vertex of
the convex hull of the residuals of length at most $L$.
\end{theorem}
\begin{proof}
Use a state for every positive-mass prefix of length at most $L$, initialized at
$\varepsilon$. At state $u$, retain the source's stopping and letter
probabilities. For $|u|<L$, route letter $a$ to $ua$. For
$u=u_1\cdots u_L$, route it to $u_2\cdots u_La$. The latter source
state and the true successor are obtained by the \emph{same} sequence
of $L$ normalized letter maps, applied to $x_0$ and $x_{u_1}$.
The first map bounds their distance by $\Delta$, and the remaining
$L-1$ maps contract it to $D=\Delta\gamma^{L-1}$.

To pass to probability laws, take cone-order bounds
$m y\le_{\mathcal C}x\le_{\mathcal C}M y$. Positivity implies
$m f_y\le f_x\le M f_y$. Normalization gives $m\le1\le M$;
from $M/m\le e^D$ one gets $M\le e^D$ and hence
$\|f_x-f_y\|_1\le2(e^D-1)$, by summing the positive part.
Limits give the same bound when order extrema are not attained.
No coordinate $\ell^1$ bound on an arbitrary cone is used.

Each row defect is at most $2(1-s)(e^D-1)$. The constructed automaton
stops with probability at least $s$ at every visit, so its occupation
sum is at most $1/s$. The occupation-flow certificate gives
\eqref{eq:trailing-tv}. For pruning, first merge duplicate residuals,
retain the distinct vertices, and express the initial residual and each
routed successor as convex mixtures of retained vertices. The same row
defect and stopping estimates apply. Keeping the empty-prefix state
instead is always permitted within the original $N_d(L)$ budget.
\end{proof}
The Hilbert contraction estimate is the classical Birkhoff theorem
\citep{gaubert2015}. For fixed $s,\Delta,\gamma$ and $d\ge2$, the state bound is
$N(\epsilon)=O(\epsilon^{-\log d/\log(1/\gamma)})$:
use $e^z-1\le ze^\Delta$ for $0\le z\le\Delta$ and choose
$L=O(\log(1/\epsilon))$. This is exponential enumeration in depth,
not an efficient rank-$R$ discovery theorem. Exact transition queries,
cone data, and contraction bounds are inputs to this construction.

\subsection{Relative Doeblin anchors and a corrected simplex inclusion}
\label{sec:doeblin-sandwich}
Here $\Delta(Z)$ denotes the convex hull of the rows of a matrix $Z$.
\begin{lemma}[Inverse-entry bound for a nearly diagonal stochastic matrix]
\label{lem:inverse-entry}
If $\Pi\ge0$, $\Pi\mathbf1=\mathbf1$, and $\Pi_{ii}\ge1-\eta$
with $0\le\eta<1/2$, then $\Pi$ is invertible,
$\Pi^{-1}\mathbf1=\mathbf1$, and
\[
 (\Pi^{-1})_{ij}\ge-\frac{\eta}{1-2\eta}.
\]
\end{lemma}
\begin{proof}
For $\eta>0$, write $\Pi=(1-\eta)I+\eta K$ with $K$ row-stochastic.
Expand the inverse in powers of $-\eta K/(1-\eta)$. Each entry of
$K^n$ belongs to $[0,1]$. Dropping the even nonnegative terms bounds
every entry below by the sum of the odd negative coefficients,
$-\eta/(1-2\eta)$. The case $\eta=0$ is immediate.
\end{proof}

\begin{theorem}[Conditional Doeblin simplex inclusion and strict margins]
\label{thm:doeblin-sandwich}
Let a proper $R$-state law have $\tau_*\ge s_*\mathbf1>0$,
$\alpha_{*,j}\ge\mu$, and
\[
 A_{*,a,ij}\ge\mu w_{a,i},\qquad
 w_a=A_{*,a}\mathbf1,\qquad 0<\mu<1/R.
\]
Set $a_\mu=1-R\mu$, $J_\mu=a_\mu I+\mu\mathbf1\mathbf1^\top$,
and $V_*=\widetilde P^{-1}$ in Theorem~\ref{thm:observable-certificate}.
Suppose a supplied true anchor matrix has the form
$U=\Pi J_\mu V_*$, where $\Pi$ satisfies
Lemma~\ref{lem:inverse-entry}. Write
$c_{\max}=\max_j\sum_i\Pi_{ij}$ and choose a known upper bound
$C\ge c_{\max}$. Define
\[
 m_- =\frac{\eta}{1-2\eta},\qquad
 m_+=\frac{\mu}{a_\mu C}.
\]
If $m_-\le m_+$, every $m\in[m_-,m_+]$ gives an exact proper positive
realization through
\[
 D_m=(1+Rm)I-m\mathbf1\mathbf1^\top,\quad
 V=D_mU=WV_*,\quad W=D_m\Pi J_\mu.
\]
In particular,
$\Delta(J_\mu)\subseteq\Delta(W)\subseteq\Delta(I)$ and
$W\tau_*\ge s_*\mathbf1$. If $m_-<m<m_+$, put
$\delta_m=\min(m-m_-,m_+-m)$. If also
$w_{\min}=\min_{a,i}w_{a,i}>0$, all initial and transition entries are
at least
\begin{equation}
 \gamma_* =\frac{\delta_m}{1+Rm}\min(1,w_{\min})>0.
 \label{eq:sandwich-margin}
\end{equation}
One may use $C=1+(R-1)\eta$. A sufficient nonempty-interval condition is
\begin{equation}
 a_\mu(R-1)\eta^2+[1-(R-2)\mu]\eta-\mu\le0.
 \label{eq:sandwich-feasibility}
\end{equation}
\end{theorem}
\begin{proof}
Let $\widetilde W=D_m\Pi$. The inverse formula and row sums give
\[
 \widetilde W^{-1}=\frac{\Pi^{-1}+m\mathbf1\mathbf1^\top}{1+Rm}\ge0.
\]
It is row-stochastic. Since $J_\mu=\widetilde W^{-1}W$, this proves
the inner inclusion. For the outer inclusion one must retain column
sums $c_j$:
\[
 \widetilde W_{ij}=(1+Rm)\Pi_{ij}-mc_j,\qquad
 W_{ij}=a_\mu\widetilde W_{ij}+\mu
       \ge\mu-a_\mu m C\ge0.
\]
Also $W\mathbf1=\mathbf1$. For every row $z$ satisfying
$z_j\ge\mu(z\mathbf1)$,
\[
 (zJ_\mu^{-1})_j=\frac{z_j-\mu(z\mathbf1)}{a_\mu}\ge0.
\]
Consequently
\[
 \alpha=(\alpha_*J_\mu^{-1})\widetilde W^{-1},\qquad
 A_a=W(A_{*,a}J_\mu^{-1})\widetilde W^{-1},\qquad
 \tau=W\tau_*
\]
are nonnegative and normalized, with the asserted stopping margin.
If $m$ is strictly inside the interval, entries of $\widetilde W^{-1}$
are at least $\delta_m/(1+Rm)$ and $W_{ij}\ge a_\mu C\delta_m>0$.
The vector $\alpha_*J_\mu^{-1}$ has mass one. The row sums of
$W A_{*,a}J_\mu^{-1}$ equal $(Ww_a)_i\ge w_{\min}$, which proves
\eqref{eq:sandwich-margin}. Finally $c_j\le1+(R-1)\eta$; substituting
this bound into $m_-\le m_+$ gives \eqref{eq:sandwich-feasibility}.
\end{proof}
The margin in \eqref{eq:sandwich-margin} uses the stochasticity of $W$;
it does not require a cross-row letter-probability ratio. The hypotheses
concern relative anchors inside $\Delta(J_\mu)$, not nearly pure hidden
states in $\Delta(I)$. They are additional structural assumptions.

\subsection{Selection from a supplied relatively separable prefix pool}
\label{sec:relative-spa}
\begin{theorem}[Conditional relative SPA selection, with corrected noise budget]
\label{thm:relative-spa}
Let $G=J_\mu V_*$ be invertible, and let a supplied pool of $M$ true
prefix-coordinate rows form $F=HG$, where $H\ge0$ and
$H\mathbf1=\mathbf1$. Suppose there are distinct witness rows $i_r$
with $H_{i_r,r}\ge1-\eta_0$. Let
$\max_i\|\widehat F_{i,:}-F_{i,:}\|_2\le e_F$ and define
\begin{align*}
 \sigma&=\sigma_R(G),\qquad \kappa=\|G\|_2/\sigma,
 & C_{\rm SPA}&=1+80\kappa^2,\\
 e&=e_F+\sqrt2\eta_0\|G\|_2,
 & \theta_R&=\min\{(2\sqrt{R-1})^{-1},1/4\}.
\end{align*}
For $R\ge2$, if $e<\sigma\theta_R/C_{\rm SPA}$, ordinary SPA on
$\widehat F^\top$ returns $R$ indices in $O(MR^2)$ exact-real arithmetic
operations. After a permutation, their \emph{true} relative coordinates
satisfy
\begin{equation}
 H_{\widehat i_r,r}\ge1-\eta,\qquad
 \eta=\frac{e_F+C_{\rm SPA}e}{\sigma}.
 \label{eq:spa-relative-budget}
\end{equation}
If $\eta<1/2$ and \eqref{eq:sandwich-feasibility} is strict, these true
selected anchors satisfy Theorem~\ref{thm:doeblin-sandwich}.
\end{theorem}
\begin{proof}
Replace the witness rows of $F$ by the corresponding rows of $G$,
leaving other rows unchanged, to obtain an exactly separable matrix
$F^0$. The replacement changes each row by at most
$\sqrt2\eta_0\|G\|_2$, since a simplex vector $h$ with
$h_r\ge1-\eta_0$ obeys $\|h-e_r\|_2\le\sqrt2\eta_0$.
Apply the squared-Euclidean specialization of
\citet[Theorem 3 and Remark 2]{gillis2014} to $(F^0)^\top$ and
$\widehat F^\top$, using the larger factor $C_{\rm SPA}$.
Its hypotheses allow signed generator coordinates. The selected
\emph{observed} rows are within $C_{\rm SPA}e$ of distinct rows of $G$.
Their true rows are within $e_F+C_{\rm SPA}e$. Multiplication by
$G^{-1}$ bounds each relative-coordinate error by
$(e_F+C_{\rm SPA}e)/\sigma$, proving \eqref{eq:spa-relative-budget}.
\end{proof}
This is an application of a classical selection theorem, not a new SPA
algorithm. Its cost starts with the pool matrix already available.
The pool size, generation cost, prefix masses, and uniform estimation
error are separate inputs. The estimated rows $\widehat F_{\widehat i_r,:}$
are not the unknown true rows $F_{\widehat i_r,:}$; replacing the latter
in $V=D_mU$ requires an additional numerical error certificate.

\begin{proposition}[Uniform observable coordinates for a fixed prefix pool]
\label{prop:pool-estimation}
Use a fixed core with $|u_i|,|v_j|\le L_0$ and
$\sigma_R(H)\ge\sigma_0>0$, and prefixes in a fixed
pool $\mathcal P\subseteq\Sigma^{\le L_0}$ with masses at least $p_*>0$.
From $n$ independent words, observe the first $2L_0$ letters together
with whether and when the word terminates by that depth. Estimate
$H$, the core prefix masses $q$, the pool prefix masses $q_w$, and
$h_{w,j}=p(wv_j)$ by counts. Put
\[
 K=R^2+R+MR+M,\quad
 t=\sqrt{\log(2K/\delta)/(2n)},\quad h_0=\sqrt R/\sigma_0.
\]
If $Rh_0t\le1/2$ and $t\le p_*/2$, then with probability at least
$1-\delta$, all estimates
\[
 \widehat f_w=(\widehat h_w/\widehat q_w)\widehat H^{-1}D_{\widehat q}
\]
exist and have row $\ell^2$ error at most
\begin{equation}
 e_F\le\frac{8Rh_0t}{p_*}+(2Rh_0^2+h_0)t.
 \label{eq:pool-estimation}
\end{equation}
\end{proposition}
\begin{proof}
A union bound and Hoeffding give error at most $t$ for all $K$ entries.
The Neumann bound gives $\|\widehat H^{-1}\|_\infty\le2h_0$ and
$\|\widehat H^{-1}-H^{-1}\|_\infty\le2Rh_0^2t$.
Both conditional rows $h_w/q_w$ and $\widehat h_w/\widehat q_w$ are
nonnegative with mass at most one, because the $v_j$ are distinct.
Each coordinate differs by at most $4t/p_*$. Telescope
$(h_w/q_w)H^{-1}D_q$; the conditional-row error contributes at most
$8Rh_0t/p_*$ in $\ell^1$, and the remaining product error contributes
at most $(2Rh_0^2+h_0)t$. The $\ell^2$ norm is no larger.
\end{proof}
The uniform event permits choosing indices by SPA on the same pool.
It does not permit choosing the pool or core after observing these
data without an additional uniform argument or sample split.
