\subsection{Reversible integer lifting with an actual coding contract}
\label{sec:integer-lifting}
Integer lifting is classical \citep{calderbank1998}. Its relevance here
is that it repairs injectivity of the proposed residual transform;
it does not by itself establish a Hankel-specific compression gain.

\begin{lemma}[Constructing lifting factors from a supplied orthogonal basis]
\label{lem:rotation-lifts}
A supplied computable $U\in SO(D)$ can be approximated arbitrarily
closely by a product of at most $3D(D-1)$ rational shears. A planar
rotation of angle $\theta$ is
\[
 \begin{pmatrix}1&-\tan(\theta/2)\\0&1\end{pmatrix}
 \begin{pmatrix}1&0\\\sin\theta&1\end{pmatrix}
 \begin{pmatrix}1&-\tan(\theta/2)\\0&1\end{pmatrix}.
\]
Split angles if necessary so that $|\theta|\le\pi/2$; each shear
coefficient then has magnitude at most one. If all $\ell$ coefficients
are rounded within $\delta$ while staying in $[-1,1]$, the product
error in maximum-row-sum norm is at most $\ell2^{\ell-1}\delta$.
\end{lemma}
\begin{proof}
Givens elimination expresses a special orthogonal matrix as at most
$D(D-1)/2$ planar rotations, allowing angles up to $\pi$.
Splitting each into at most two proves the shear count by the displayed
identity. Telescope the product; all factors before and after rounding
have norm at most two, and each factor error has norm at most $\delta$.
\end{proof}

\begin{theorem}[Bijective lifting and a covariance bound including correlated errors]
\label{thm:integer-lifting}
Let $S_r=I+s_re_{i_r}e_{j_r}^\top$, $i_r\ne j_r$, be rational
shears, $1\le r\le\ell$, and use a fixed deterministic nearest-integer
rule. Apply $z_{i_r}\leftarrow z_{i_r}+\operatorname{round}(s_rz_{j_r})$
in the stated order. The resulting map $T:\mathbb Z^D\to\mathbb Z^D$
is a bijection, whose inverse applies the same rounded values with
subtraction in reverse order. Put $U_*=S_\ell\cdots S_1$ and
\[
 b_j=\frac12\sum_{r=1}^\ell
       |(S_\ell\cdots S_{r+1})_{j,i_r}|.
\]
For a random integer vector $Z$ of finite second moment and covariance
$M$, the output obeys
\begin{equation}
 \operatorname{Var}(T(Z)_j)
 \le\left(\sqrt{(U_*MU_*^\top)_{jj}}+b_j\right)^2.
 \label{eq:lifting-variance}
\end{equation}
If a desired orthogonal $U$ diagonalizes $M$ with diagonal entries
$\lambda_j$, the right side is at most
\[
 \left(\sqrt{\lambda_j}+
    \|(U_*-U)_{j,:}\|_2\sqrt{\tr M}+b_j\right)^2.
\]
These assertions apply separately to conditional laws given a core.
\end{theorem}
\begin{proof}
The donor coordinate $j_r$ is unchanged by a shear, so subtraction
recovers its integer target exactly. Reverse composition proves
bijectivity, for any deterministic rounding convention.
Writing the rounding error at step $r$ as $\epsilon_r e_{i_r}$,
$|\epsilon_r|\le1/2$, gives
\[
 T(z)=U_*z+\sum_{r=1}^\ell
            S_\ell\cdots S_{r+1}e_{i_r}\epsilon_r(z).
\]
Thus the $j$th error lies in $[-b_j,b_j]$ and its standard deviation
is at most $b_j$. The triangle inequality in centered $L^2$ proves
\eqref{eq:lifting-variance}, without independence. For the last
bound, split $U_*=U+(U_*-U)$ and apply the same inequality, using
$\operatorname{Var}(vZ)\le\|v\|_2^2\tr M$.
\end{proof}
The claim that each lift adds at most $1/4$ to the \emph{output}
variance is false. Let $Z_1=Z_2$ be Bernoulli$(1/2)$ and lift with
$s=51/100$. The exact linear target has variance $(151/100)^2/4$,
whereas the integer target is $2Z_1$, with variance one. The increase
is $17199/40000>1/4$, because the rounding error is correlated with
the signal. The propagated-error terms above are indispensable.

\begin{theorem}[Complete expected length against a specified baseline]
\label{thm:actual-lifting-code}
Let $S$ be a finite-bit integer vector and $C=c(S)$ a core. Given $C$,
an integer predictor and a specified lifting transform produce
$Z'=T_C(S-D(C))$. The decoder receives the core and every model,
predictor, transform and framing description needed for inversion;
their total nonpayload length is at most $B_{\rm meta}$ bits.
Fix the coding models before coding the evaluated source; any dependence
on the transmitted core is explicitly conditional. Let prefix codes
for $C$ and for each $Z'_j$ conditional on $C$ have
Shannon lengths for specified full-support models $Q_C,Q_{j|C}$.
Use base-two KL and put
\[
 D_C=D_2(P_C\|Q_C),\qquad
 D_j(c)=D_2(P_{Z'_j|c}\|Q_{j|c}).
\]
If $v_j(c)$ is any valid bound on the conditional variance, the
expected complete length is at most
\begin{align}
 \mathbb E L\le B_{\rm meta}+H(C)+D_C+1+D
 +\mathbb E_C\sum_{j=1}^D\left[
 \tfrac12\log_2(2\pi e(v_j(C)+1/12))+D_j(C)\right].
 \label{eq:actual-lifting-code}
\end{align}
If this explicit upper bound is smaller than a proved lower bound
on the expected length of a \emph{specified} baseline coder, the new
coder improves that baseline in expectation at exact reconstruction.
\end{theorem}
\begin{proof}
The inverse lifting steps and the integer predictor reconstruct $S$
exactly from $C,Z'$. Expected Shannon lengths are bounded by cross
entropy plus one for each code: cross entropy is entropy plus the
displayed model mismatch. Proposition~\ref{prop:discrete-entropy-bound}
bounds each marginal conditional entropy by its variance term.
Sum these inequalities, average over the core, and add every
nonpayload description. The final comparison is then immediate.
\end{proof}
This is a coding theorem for the declared models, not an assertion
that their probabilities or description lengths are free to obtain.
An implementation must also specify finite integer widths. For
bounded source support, the lifting recurrence bounds every intermediate
magnitude by $\|S_r\|_\infty$ times the previous bound plus $1/2$;
these explicit bounds determine exact numerator, product and accumulator
widths. In-place lifting stores $D$ integers plus multiplication scratch;
holding additional input/output buffers must be charged. No reduction
of update latency or peak memory is implied by \eqref{eq:actual-lifting-code}.

\begin{proposition}[The universal saving threshold is impossible]
\label{prop:no-universal-code-gain}
For the deterministic exact representation above,
$H(C,Z')=H(S)$. Consequently a valid complete-code upper bound cannot
lie below a valid lower bound on $H(S)$. In particular the patch's
criterion for beating \emph{every} lossless coder by placing its
complete length below an entropy lower bound has no realizable instance.
\end{proposition}
\begin{proof}
The encoder is deterministic and injective, and has a deterministic
inverse on its image; entropy is invariant. Any uniquely decodable
code has expected length at least $H(S)$. If $U$ bounds its complete
length above and $L$ bounds $H(S)$ below, then
$U\ge\mathbb E L_{\rm code}\ge H(S)\ge L$, so $U<L$ is impossible.
\end{proof}
For comparison, let $S=(Z,Z)$ with $Z$ uniform on $b$-bit integers.
The lift $(z_1,z_2)\mapsto(z_1,z_2-z_1)$ yields $(Z,0)$.
With the distribution and lift fixed, $b$ payload bits suffice instead
of a raw $2b$-bit layout; including a header of length $H_0$ still
saves bits if $H_0<b$. The joint entropy was already $b$, so this is
a gain over the specified raw layout, not over an optimal lossless code.

For a practical fixed frame, a one-bit selector between a new exact
codec and an exact baseline gives length at most
$1+\min(L_{\rm new},L_{\rm base})$, with all headers included.
This provides a deterministic fallback guarantee; it does not supply
held-out evidence for learned Hankel predictors.
