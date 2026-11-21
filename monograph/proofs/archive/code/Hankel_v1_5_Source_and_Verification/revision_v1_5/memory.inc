\subsection{An optimal table-free executable on homogeneous separated attractors}
\label{sec:address-memory}
\begin{theorem}[Homogeneous addresses with all storage and update costs charged]
\label{thm:address-memory}
Suppose $\mathcal X\subset\mathbb R^k$ is the attractor of $d\ge2$
similarities $\tau_a(x)=\lambda O_ax+t_a$ of a common ratio
$0<\lambda<1$, with pairwise separated images. Let the dimension
$k$ and the rational affine maps be fixed. Assume $\mathcal X$ has
positive diameter, $x_0$ is the fixed point of $\tau_0$, and the
future-observability inequalities of Theorem~\ref{thm:ahlfors-memory}
hold for an output $h$ admitting certified $b$-bit evaluation in
polynomial time and $O_H(b)$ scratch bits. Then there is a table-free
finite-state executable with error at most $2^{-p}$ and
\[
 B(p)=\frac{\log d}{\log(1/\lambda)}p+O_H(1)
       =B_f^*(2^{-p})+O_H(1)
\]
persistent state bits. Its symbol update costs $O_H(p)$ bit operations.
For fixed-dimensional rational affine maps, readout takes
$O_H(p\mathsf M(p))$ bit operations plus the cost of evaluating $h$,
with $O_H(p)$ scratch and immutable parameter bits. No exponential
transition table is stored.
\end{theorem}
\begin{proof}
Put $D=\operatorname{diam}\mathcal X$ and choose
\[
 m=\max\left\{1,
 \left\lceil\frac{p+\log_2(2C_{\rm obs}D)}{\log_2(1/\lambda)}
                                      \right\rceil\right\}.
\]
The persistent state is a word $u$ of length $m$, initialized to $0^m$;
on a letter append it and discard the oldest symbol. Its mathematical
representative is $x(u)=\tau_{u_m}\circ\cdots\circ\tau_{u_1}(x_0)$.
Before $m$ inputs, the initial zeros fix $x_0$, so this representative
equals the true state. After $m$ inputs both states share the same
last $m$ maps, giving distance at most $D\lambda^m$ and output
error at most $2^{-p-1}$.

Pack the word as a base-$d$ integer $I\in\{0,\ldots,d^m-1\}$.
The exact update is $I\leftarrow(dI+a)\bmod d^m$.
It uses $\lceil m\log_2d\rceil$ persistent bits, without a startup
counter; $m,d^m$ and the program are immutable. Fixed-$d$ multiplication,
comparison and subtraction on this integer cost $O_H(p)$ bit
operations, not constant time on a machine with a fixed word width.
To read out, copy the index into scratch, recover its $m$ symbols,
and evaluate the $m$ affine maps. Per-step certified rounding error
$\epsilon_r$ contributes at most $\epsilon_r/(1-\lambda)$.
Choosing $O_H(p)$ fractional bits and a final output error at most
$2^{-p-2}$ leaves total numerical error at most $2^{-p-1}$.
The stated readout and storage costs follow for fixed-dimensional
affine maps and the stipulated output evaluator.

For completeness the separated homogeneous attractor is Ahlfors
regular of dimension $s=\log d/\log(1/\lambda)$: assign mass $d^{-m}$
to each depth-$m$ cylinder. A cylinder containing $x$ and of diameter
less than $r$ supplies the lower ball-mass bound. At a depth with
$\lambda^m$ comparable to $r$, distinct cylinders have separation
bounded below by a constant times $\lambda^m$; Euclidean packing
bounds the number intersecting a ball of radius $r$ by a model
constant, giving the upper bound. Reachable states are dense, since
the images of $x_0$ lie in every cylinder at every depth.
Theorem~\ref{thm:ahlfors-memory} supplies the matching lower bound.
\end{proof}
\begin{theorem}[An effective contractive coordinate chart]
\label{thm:chart-memory}
Suppose predictive states admit a supplied global chart
$\phi:[0,1]^k\to\mathcal X$ such that the conjugate updates
$F_a=\phi^{-1}\tau_a\phi$ map the cube to itself and are
$\lambda$-Lipschitz in maximum norm, $\lambda<1$.
Assume reachable states are dense and
\[
 c\|z-z'\|_\infty\le d_{\rm fut}(\phi z,\phi z')
                         \le C\|z-z'\|_\infty.
\]
Suppose the initial coordinate, the $F_a$ and $h\circ\phi$ have
certified $q$-bit evaluators using polynomial time, $O_H(q)$ scratch
bits and $O_H(q)$ immutable parameter bits, with fixed $k,d$.
Then an executable without a transition table attains error $2^{-p}$
with $kp+O_H(1)=B_f^*(2^{-p})+O_H(1)$ persistent bits and polynomial
compilation, update and readout bit cost.
\end{theorem}
\begin{proof}
Set $\epsilon=2^{-p}$ and choose a mesh $\delta=2^{-q}$ with
$\delta\le(1-\lambda)\epsilon/(4C)$, so $q=p+O_H(1)$.
Store the $k$ integer grid coordinates in $\{0,\ldots,2^q\}$.
Evaluate each update to maximum-norm error at most $\delta/2$,
clip to the cube and round to the nearest grid point. Clipping is
nonexpansive; the total one-step numerical error is at most $\delta$.
Initialize with error at most $\delta$. Thus
$e_{t+1}\le\lambda e_t+\delta$ and
$e_t\le\delta/(1-\lambda)$ for every history.
Output error from state tracking is at most $\epsilon/4$; evaluating
the readout to error at most $\epsilon/2$ meets the target.
Each coordinate needs $q+1$ bits, giving the asserted upper bound
and the stated computational costs by the supplied evaluators.
A grid packing of the cube at spacing greater than $2\epsilon/c$,
and reachable density, forces $\Omega_H(\epsilon^{-k})$ distinct
executable states by future observability. This gives the lower bound.
\end{proof}
The chart theorem assumes contraction in the executable coordinates;
bi-Lipschitz equivalence alone can enlarge a contraction constant
beyond one. Discovering such a chart is not part of the theorem.

The fixed-point initialization in Theorem~\ref{thm:address-memory}
can be replaced by an effectively
supplied length-$m$ address of $x_0$ whose nested cylinder contains it,
with its acquisition cost charged. Arbitrary initial coordinates do
not automatically supply such an address.

\paragraph{Why unequal ratios are a different target.}
For $\tau_0(x)=x/4$ and $\tau_1(x)=x/2+1/2$, the separated attractor
has dimension $s$ determined by $4^{-s}+2^{-s}=1$, hence
$s=\log_2((1+\sqrt5)/2)<1$.
A fixed-length register based on the worst contraction $1/2$ stores
$p+O(1)$ bits, rather than $sp+O(1)$. Thus the patch's optimal
coefficient does not follow for unequal ratios. Uniform bi-Lipschitz
bounds alone also do not prove Ahlfors regularity for arbitrary
nonconformal systems. Adaptive address coding or effective charts
need separate theorems with their own update and storage costs.
