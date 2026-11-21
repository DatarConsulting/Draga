\subsection{Packed startup states remove the fixed-point initialization restriction}
\label{sec:packed-startup}
\begin{theorem}[Homogeneous addresses from an arbitrary effective initial state]
\label{thm:packed-startup}
Retain the separated homogeneous similarities and future-observability
assumptions of Theorem~\ref{thm:address-memory}, but replace its
fixed-point assumption by an arbitrary $x_0\in\mathcal X$ admitting
certified $q$-bit evaluation in polynomial time and $O_H(q)$ scratch
and immutable parameter bits. Assume also that $h$ has an effectively
evaluable $C_{\rm obs}$-Lipschitz extension to a fixed neighborhood
of $\mathcal X$, with the same evaluation resource bounds as before.
Then the same optimal persistent memory
\[
 B(p)=\frac{\log d}{\log(1/\lambda)}p+O_H(1)
       =B_f^*(2^{-p})+O_H(1)
\]
is achieved without a transition table or a supplied address of $x_0$.
Updates use $O_H(p)$ bit operations. Readout uses
$O_H(p\mathsf M(p))$ bit operations plus evaluation of $x_0$ and $h$,
and $O_H(p)$ scratch and immutable parameter bits.
\end{theorem}
\begin{proof}
Choose $m$ as in Theorem~\ref{thm:address-memory}, and retain the
entire input word until length $m$, and thereafter its last $m$
letters. The state set is $Q_m=\bigcup_{\ell=0}^m\Sigma^\ell$.
Its cardinality is $(d^{m+1}-1)/(d-1)$, so it uses
$m\log_2d+O_d(1)$ bits. Storing a separate length counter is
unnecessary: for a length-$\ell$ word of base-$d$ value $v$, put
\[
 S_\ell=\frac{d^\ell-1}{d-1},\qquad I=S_\ell+v.
\]
The length classes occupy the disjoint consecutive intervals
$[S_\ell,S_{\ell+1})$. Initialize $I=0$. On input
$a\in\{0,\ldots,d-1\}$, the exact update is
\[
 I\leftarrow
 \begin{cases}
 dI+a+1,&I<S_m,\\
 S_m+\big(d(I-S_m)+a\big)\bmod d^m,&I\ge S_m.
 \end{cases}
\]
This needs only the immutable $O_H(p)$-bit constants $S_m,d^m$
and fixed-$d$ integer arithmetic; hence the bit cost is $O_H(p)$.
No time-dependent external counter is used.

The representative of a stored word $u$ is $\tau_u(x_0)$, in the
sequential convention of Theorem~\ref{thm:address-memory}.
Before the window fills, this is the true state. Afterwards the
representative and the true state share their last $m$ contractions,
and so are within $\lambda^m\operatorname{diam}\mathcal X$.
The truncation output error is at most $2^{-p-1}$.

For readout, a scratch copy of $I$ decodes the letters backwards by
$a=(I-1)\bmod d$ and $I\leftarrow\lfloor(I-1)/d\rfloor$
until zero; a scratch word stores the reversed letters for replay.
There are at most $m=O_H(p)$ such operations, each on $O_H(p)$ bits.
Evaluate $x_0$ and replay the affine maps with initial error
$e_0\le 2^{-p}/(8C_{\rm obs})$ and per-step error
$e_r\le(1-\lambda)2^{-p}/(8C_{\rm obs})$.
The total state error is at most
$e_0+e_r/(1-\lambda)\le2^{-p}/(4C_{\rm obs})$.
A final readout error at most $2^{-p-2}$ completes the budget.
Fixed model constants ensure the computation stays in the stipulated
neighborhood; finitely many small precisions can use finer arithmetic.
The stated storage and readout bounds follow.
For every $x_0\in\mathcal X$, its finite images meet every cylinder,
so reachability is dense. The same packing lower bound as in
Theorem~\ref{thm:address-memory} applies.
\end{proof}
This removes an initialization restriction within the homogeneous
class. It does not resolve unequal contraction ratios or discover
effective coordinates for a general regular attractor.
