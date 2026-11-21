\subsection{Normalized transfer convergence on a strict log-Lipschitz cone}
\label{sec:normalized-transfer}
The endpoint value $L=c/(1-\gamma)$ suggested in the new patch gives
only cone invariance. Strict inclusion requires a strictly larger
$L$. Moreover the relevant order is that of the regularity cone,
not only pointwise positivity in $C([-1,1])$.

\begin{theorem}[Quantitative normalized positive eigenfunction]
\label{thm:normalized-transfer}
On $I=[-1,1]$, let
$(\mathcal L f)(x)=\sum_{a=1}^d w_a(x)f(G_a(x))$, where
$G_a:I\to I$ have Lipschitz constants at most $\gamma<1$,
$0<w_-\le w_a(x)\le w_+$, and
$|\log w_a(x)-\log w_a(y)|\le c|x-y|$.
Choose $L>c/(1-\gamma)$ and put $L'=c+\gamma L<L$.
Define
\[
 t_* =\frac{L-L'}{L+L'}e^{-2(L+2L')},\qquad
 \Delta=2\log(1/t_*),\qquad \vartheta=\tanh(\Delta/4)<1.
\]
There is a unique eigenfunction $h>0$ in the log-Lipschitz cone of
constant $L$, normalized by $h(0)=1$. It satisfies
\[
 \mathcal L h=\lambda h,\quad
 e^{-L'}\le h(x)\le e^{L'},\quad
 d w_-e^{-L'}\le\lambda\le d w_+e^{L'}.
\]
Writing $V_n=\mathcal L^n1$, for $n\ge1$,
\begin{equation}
 \sup_{x\in I}\left|\frac{V_n(x)/V_n(0)}{h(x)}-1\right|
       \le e^{\Delta\vartheta^{n-1}}-1.
 \label{eq:normalized-transfer}
\end{equation}
\end{theorem}
\begin{proof}
Let $K_L$ be the closed convex cone of nonnegative Lipschitz functions
with $f(x)\le e^{L|x-y|}f(y)$ for all $x,y$; its only function
vanishing at a point is zero. Each positive summand of $\mathcal L f$
has log-Lipschitz constant at most $c+\gamma L$, and their sum has
the same bound. Thus $\mathcal L(K_L\setminus\{0\})\subset K_{L'}$.

For $f,g\in K_{L'}$ normalized at zero, $g(y)/f(y)\ge e^{-2L'}$,
and hence also $g(y)/f(y)\ge e^{-4L'}$. For $r=|x-y|\in(0,2]$,
\[
 \frac{e^{Lr}-e^{L'r}}{e^{Lr}-e^{-L'r}}
 \ge\frac{L-L'}{L+L'}e^{-2L}.
\]
This follows by expressing numerator and denominator as integrals
of $r e^{ur}$. Therefore
$e^{Lr}g(y)-g(x)\ge
t_*[e^{Lr}f(y)-f(x)]$. Pointwise positivity also holds, so
$g-t_*f\in K_L$, and symmetrically $f-t_*g\in K_L$.
The image diameter in the Hilbert metric of $K_L$ is at most
$2\log(1/t_*)=\Delta$.

For projective contraction, use the Lipschitz Banach space: the strict
constant $L'<L$ places these positive image functions in the interior
of $K_L$. The normalized slice $f(0)=1$ of $K_L$ is convex and compact in the
uniform norm by uniform bounds and equicontinuity. The continuous
map $f\mapsto\mathcal L f/(\mathcal L f)(0)$ has a fixed point
by Schauder. Birkhoff's projective contraction theorem
\citep{gaubert2015} gives contraction factor at most $\vartheta$ on
comparable rays. Every eigenfunction in $K_L$ lies in the image,
whose diameter is finite, so this also gives uniqueness.
Both $\mathcal L1$ and $\mathcal Lh$ are image points; iteration
gives projective distance at most $\Delta\vartheta^{n-1}$.
Cone order implies pointwise order. Since both representatives have
value one at zero, their pointwise ratios are between $e^{-D}$ and
$e^D$ when their projective distance is at most $D$. This proves
\eqref{eq:normalized-transfer}. The fixed point lies in $K_{L'}$;
its normalization gives the bounds on $h$, and evaluation at zero
gives those on $\lambda$.
\end{proof}

\begin{corollary}[A certified polynomial eigenfunction approximator]
\label{cor:normalized-polynomial}
Assume the two-radius analytic hypotheses of Theorem~\ref{thm:oracle}
and the hypotheses of Theorem~\ref{thm:normalized-transfer}, with
fixed effective margins. Suppose the section coefficients can be
enclosed to $B$ bits in time polynomial in $m,B$ and their input
description. Then, for every $p\ge1$, a terminating algorithm returns
a dyadic polynomial $\widehat h_p$, positive on $I$, together with
certified error $\|\widehat h_p-h\|_\infty\le2^{-p}$.
Its horizon, degree and coefficient precision are $O_H(p)$.
Compilation has polynomial bit cost. Horner evaluation takes
$O_H(p\mathsf M(p))$ bit operations and $O_H(p)$ scratch bits,
with $O_H(p^2)$ bits of stored polynomial coefficients; here
$\mathsf M(p)$ is the cost of multiplying $O(p)$-bit integers.
These costs concern the eigenfunction approximator, not the full R6 compiler.
\end{corollary}
\begin{proof}
Choose $n=O_H(p)$ so that the right side of
\eqref{eq:normalized-transfer}, multiplied by $e^{L'}$, is at most
$2^{-p-2}$. Let $a_*=2g_-^{s_+}$ be the lower growth bound of
Theorem~\ref{thm:oracle}. For a rounded polynomial iterate with
$E_n=\|\widehat P_n-V_n\|_\infty\le a_*^n/2$, normalization gives
\[
 \left\|\frac{\widehat P_n}{\widehat P_n(0)}-
                  \frac{V_n}{V_n(0)}\right\|_\infty
 \le 2(1+e^{L'})E_n/a_*^n.
\]
Indeed $V_n/V_n(0)\in K_{L'}$ normalized at zero and
$\widehat P_n(0)\ge a_*^n-E_n$.
Take $E_n\le a_*^n2^{-p-2}/[2(1+e^{L'})]$.
Equation~\eqref{eq:transfer-rounded-error}, with
$\delta_m=O_H(R^{-m})$, then permits $m=O_H(p)$ and per-step
certified arithmetic error $2^{-O_H(p)}$. The polynomial coefficient
norms grow at most exponentially in $n$ and normalization divides
by a quantity at least $a_*^n/2$, also exponential in $n$.
Thus $O_H(p)$ coefficient bits suffice after allowing for the
$O(p)$ operations per row and final dyadic rounding of the normalized
coefficients, whose error is included in the budget. This can be
implemented by interval arithmetic with directed dyadic rounding
and the stipulated coefficient enclosures; double precision until
the explicit error tests pass. The preceding bounds prove termination.

If necessary request the smaller final error
$\min(2^{-p},e^{-L'}/2)$; this changes precision only by a model
constant and proves positivity from $h\ge e^{-L'}$.
There are $O(nm^2)=O_H(p^3)$ dense scalar operations after section
construction, on $O_H(p)$-bit numbers; this and the coefficient oracle
give polynomial compilation cost. There are $O(p)$ coefficients of
$O(p)$ bits. Horner evaluation at a rational argument known to the
required $O_H(p)$ bits uses $O(p)$ multiplications; the coefficient
norm bound controls error accumulation. Its working accumulator uses
$O_H(p)$ bits. Immutable coefficients are charged separately.
\end{proof}
For the homographic weights of Section~\ref{sec:oracle}, one can use
$c=\max_a 2s_+|\xi_a|/(1-|\xi_a|)$ and a verified branch derivative
bound $\gamma<1$. Rational fixed exponents and rational positive
parameters give effective power enclosures by standard interval
root or logarithm/exponential algorithms. Alternatively their
certified coefficient oracle is an explicit input to the corollary.

\paragraph{What this closes for R6.}
The normalized target, a positive lower bound, geometric convergence,
and its polynomial approximation are now supplied under explicit
effective hypotheses. A state compiler has additional acceptance
conditions, inverse-coordinate operations and inductive update
relations. The above theorem does not verify them. In particular its
fast Horner cost cannot be relabeled the cost of every R6 state update.
