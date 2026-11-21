\subsection{Pathwise stability and a complete fixed-basis sampling construction}
\label{sec:fixed-basis-sampling}
The following construction uses exact rational arithmetic for intermediate
linear algebra and a specified dyadic output. It is a finite-bit algorithm,
but it is not a theorem about unrestricted $b$-bit floating-point LU.

\begin{proposition}[Pathwise perturbation and normalized initial estimates]
\label{prop:pathwise-stability}
Let $p,\widehat q$ be proper row-normalized positive automata on the same
$R$ states, with $\tau\ge s\mathbf1>0$ for the first automaton. Then
\begin{equation}
 \TV(p,\widehat q)\le\frac12\|\widehat\alpha-\alpha\|_1
 +\frac1{2s}\left(\|\widehat\tau-\tau\|_\infty+
                   \sum_a\|\widehat A_a-A_a\|_\infty\right).
 \label{eq:pathwise-general}
\end{equation}
If $\alpha_i,A_{a,ij}\ge\gamma>0$, $R\ge2$, and an unnormalized
estimate $a^0$ satisfies $\|a^0-\alpha\|_\infty\le\delta_\alpha
\le\gamma/[2(R+1)]$, its normalization $\bar\alpha=a^0/(a^0\mathbf1)$
is positive and
\begin{equation}
 \|\bar\alpha-\alpha\|_\infty
 \le\frac{(R-1)\delta_\alpha}{1-R\delta_\alpha}\le\gamma/2.
 \label{eq:initial-normalization}
\end{equation}
For transition estimates with per-letter maximum-entry errors $e_a$,
$\max_a e_a\le\gamma/2$ and $R\sum_a e_a\le s/2$ suffice for
$\bar A_{a,ij}\ge\gamma/2$ and
$\bar\tau=\mathbf1-\sum_a\bar A_a\mathbf1\ge(s/2)\mathbf1$.
\end{proposition}
\begin{proof}
For each $w=a_1\cdots a_m$, telescope
\begin{align*}
 \widehat\alpha\widehat A_w\widehat\tau-\alpha A_w\tau
 ={}&(\widehat\alpha-\alpha)\widehat A_w\widehat\tau
       +\alpha A_w(\widehat\tau-\tau)\\
 &+\sum_{k=1}^m\alpha A_{a_1\cdots a_{k-1}}
       (\widehat A_{a_k}-A_{a_k})
       \widehat A_{a_{k+1}\cdots a_m}\widehat\tau.
\end{align*}
The suffix sum $\widehat g=\sum_v\widehat A_v\widehat\tau$ is at most
$\mathbf1$ componentwise: its finite partial sums are
$\mathbf1-\widehat T^{n+1}\mathbf1$. Equality requires termination
from every state, and is not needed here. The original occupation row
$\nu=\alpha(I-T)^{-1}$ has $\nu\mathbf1\le1/s$.
Take absolute values, apply nonnegative matrix majorants, and sum over
all prefixes, the differing letter, and suffixes. Tonelli and these
two bounds give \eqref{eq:pathwise-general}. This also follows from
Lemma~\ref{lem:stopped-coupling}.

Write $a^0=\alpha+e$. In coordinate $i$ the normalized error is
$[e_i(1-\alpha_i)-\alpha_i\sum_{j\ne i}e_j]/(1+\sum_j e_j)$.
Its numerator is at most $(R-1)\delta_\alpha$ in absolute value;
$R\gamma\le1$ proves the last inequality in
\eqref{eq:initial-normalization}. Positivity of $a^0$ and of its sum
is immediate. Finally each row's total transition error is at most
$R\sum_a e_a$; recomputing STOP as the remainder proves the last claim.
\end{proof}

\begin{theorem}[Fixed rational basis: iid estimation, exact computation, dyadic output]
\label{thm:fixed-basis-sampling}
Fix the core in Theorem~\ref{thm:observable-certificate}, with
$|u_i|,|v_j|\le L_0$, $\sigma_R(H)\ge\sigma_0>0$, and
$q_i\ge q_*>0$. A \emph{supplied rational} invertible $V$ satisfies
$V\mathbf1=\mathbf1$ and certifies the true law with
$\alpha_i,A_{a,ij}\ge\gamma>0$, $\tau_i\ge s>0$, where the
supplied lower bounds $\gamma,s$ are rational.
The core and $V$ are fixed independently of the sample (or obtained
from an independent training sample with these true margins certified).
For $0<\delta<1$, observe $n$ iid words censored at $D=2L_0+1$, retaining whether and
when termination occurs by $D$, including at $D$ itself.

Estimate $H,M_a,h$ by exact completed-word counts and $q$ by prefix
counts, each divided by $n$. All subsequent intermediate operations
in this theorem are exact rational operations. Define
\begin{align*}
 t_n&=c_0/\sqrt n,&
 c_0&=\sqrt{\tfrac12\log(4R(R+2)/\delta)},\\
 b_n&=c_1/\sqrt n,&
 c_1&=\sqrt d+\sqrt{2\log(2R^2/\delta)},\\
 h_0&=\sqrt R/\sigma_0,&v_-&=\|V^{-1}\|_\infty,
 \qquad v_+=\|V\|_\infty,\\
 K_\alpha&=v_-[(2R+1)h_0+2Rh_0^2],\\
 K_0&=\frac{2Rh_0}{q_*^2}
          +\frac{4R^2h_0^2+4Rh_0}{q_*},\\
 E_\alpha&=K_\alpha c_0,&
 E_A&=v_+v_-\left(\frac{2h_0R^2}{q_*}c_1+K_0c_0\right).
\end{align*}
It suffices to take
\begin{equation}
 n\ge\max\left\{(2Rh_0c_0)^2,(2c_0/q_*)^2,
 (4E_\alpha/\gamma)^2,
 \left(2E_A/\min(\gamma,s)\right)^2\right\}.
 \label{eq:fixed-basis-n}
\end{equation}
Check that $\widehat H$ is invertible and all $\widehat q_i>0$;
otherwise return \texttt{UNRESOLVED}. Form
\begin{align*}
 a^0&=\widehat h\widehat H^{-1}D_{\widehat q}V^{-1},
 &\bar\alpha&=a^0/(a^0\mathbf1),\\
 \bar A_a&=V D_{\widehat q}^{-1}\widehat M_a
                  \widehat H^{-1}D_{\widehat q}V^{-1},
 &\bar\tau&=\mathbf1-\sum_a\bar A_a\mathbf1.
\end{align*}
Before dividing by $a^0\mathbf1$, check that this sum is positive.
Check normalization and the margins
$\bar\alpha_i,\bar A_{a,ij}\ge\gamma/2$ and
$\bar\tau_i\ge s/2$ exactly; on failure return
\texttt{UNRESOLVED}. The supplied checking thresholds can be taken
rational. Otherwise, with $h_b=2^{-b}\le\gamma/4$, round
every transition entry down to a multiple of $h_b$ and assign the
remaining row mass exactly to STOP. Round the first $R-1$ entries of
$\bar\alpha$ down and assign the remaining mass to its last entry.
The resulting dyadic law $q_b$ has exactly $R$ positive states and,
with probability at least $1-\delta$, is returned and satisfies
\begin{equation}
 \TV(p,q_b)\le
 \frac{e_\alpha}{1-e_\alpha}+(R-1)h_b
          +\frac{e_A+dRh_b}{s},\qquad
 e_\alpha=E_\alpha/\sqrt n,\quad e_A=E_A/\sqrt n.
 \label{eq:fixed-basis-tv}
\end{equation}
In particular the rate is
$O_{R,\sigma_0,q_*,V,s,\gamma}
(\sqrt{(d+\log(R/\delta))/n}+d2^{-b})$.
The intermediate rational algorithm has polynomial bit complexity
in $n,D,d,R,b$ and the rational input bit length of $V$;
the output uses $O(dR^2b)$ bits for its probabilities.
\end{theorem}
\begin{proof}
For fixed $(i,j)$, the events $W=u_iav_j$ over $a$ are disjoint.
Cauchy--Schwarz bounds the expected sum of absolute frequency errors
by $\sqrt{d/n}$. Changing one observation changes that sum by at most
$2/n$. McDiarmid and a union bound over $R^2$ pairs, together with
Hoeffding for the $R^2+2R$ base entries, give a common event of
probability at least $1-\delta$ on which
\begin{equation}
 \|\widehat H-H\|_{\max},\ \|\widehat h-h\|_\infty,
 \|\widehat q-q\|_\infty\le t_n,
 \qquad
 \sum_a\|\widehat M_a-M_a\|_\infty\le R^2b_n.
 \label{eq:joint-count-event}
\end{equation}
The last inequality sums over both row and column indices. A
per-pair union bound cannot replace this step by exchanging a sum
over letters with a maximum over rows.

Since $\|H^{-1}\|_\infty\le h_0$, the first two sample thresholds give
\begin{align*}
 \|\widehat H^{-1}\|_\infty&\le2h_0,&
 \|\widehat H^{-1}-H^{-1}\|_\infty&\le2Rh_0^2t_n,\\
 \|D_{\widehat q}^{-1}\|_\infty&\le2/q_*,&
 \|D_{\widehat q}^{-1}-D_q^{-1}\|_\infty&\le2t_n/q_*^2.
\end{align*}
Distinct suffixes imply $\|h\|_1,\|\widehat h\|_1\le1$.
For each fixed core prefix, the words $u_iav_j$ over $(a,j)$ are
distinct, so
$\sum_a\|M_a\|_\infty,\sum_a\|\widehat M_a\|_\infty\le R$.
Also $\|D_q\|_\infty,\|D_{\widehat q}\|_\infty\le1$.
Telescoping the three factors in $hH^{-1}D_q$ gives
$\|a^0-\alpha\|_1\le K_\alpha t_n=e_\alpha$.
Telescoping $D_q^{-1}M_aH^{-1}D_q$ and summing over $a$ gives
\[
 \sum_a\|\bar A_a-A_a\|_\infty
 \le v_+v_-\left[
 \frac{2h_0}{q_*}\sum_a\|\widehat M_a-M_a\|_\infty+K_0t_n\right]
 \le e_A.
\]
The remaining thresholds imply $e_\alpha\le\gamma/4$ and
$e_A\le\min(\gamma,s)/2$. Thus $a^0$ is positive, and its sum is
within $e_\alpha$ of one. Its normalized entries are at least
$(\gamma-e_\alpha)/(1+e_\alpha)\ge\gamma/2$;
$\|\bar\alpha-\alpha\|_1\le2e_\alpha/(1-e_\alpha)$.
The transition entries are at least $\gamma/2$ and the recomputed
stopping probabilities at least $s/2$. All exact checks therefore pass.

Downward rounding changes at most $dR$ entries in each transition row,
with total change at most $dRh_b$; STOP absorbs this mass. Initial TV
changes by at most $(R-1)h_b$. The original source has occupation
sum at most $1/s$. Apply Lemma~\ref{lem:stopped-coupling}, or
\eqref{eq:pathwise-general} with the row-remainder identity, to obtain
\eqref{eq:fixed-basis-tv}. The rounding restriction leaves initial
and transition entries at least $\gamma/4$ and never reduces STOP.

For complexity, a trie of the $O(dR^2)$ queried words and prefixes
permits counting in $O(nD+dR^2D)$ symbol operations. Counts have
$O(\log n)$ bits. Fraction-free elimination computes the rational
inverses; determinants of the integer matrices have polynomial bit
length in $R,\log n$ and the input bit length of $V$. The remaining
$O(dR^3)$ dense rational arithmetic operations, normalization, and
integer floor divisions also have polynomial bit cost. This proves
the stated polynomial bound; it does not assert a small uniform
degree or a fixed-word-size implementation. Dyadic probabilities
are stored as their integer numerators, with remainders determined
exactly, giving the output bound.
\end{proof}

\paragraph{Scope of the numerical theorem.}
The source itself need not have rational coefficients: empirical
counts and the supplied basis are rational. Exact rational
intermediates avoid an unsupported universal backward-error constant
for partially pivoted LU. The result concerns the mathematical law
encoded by the dyadic probabilities, not bitwise equivalence to a
floating-point executable. A rational certified basis with strict
margins is a substantive input. Neither its existence for every
positive realization nor its recovery from data is asserted here.
Combining noisy SPA-selected rows with basis construction needs
control of the \emph{estimated basis}, not just selection of indices.
All these distinctions persist even though the fixed-basis sampling
and output-rounding result is now complete.
