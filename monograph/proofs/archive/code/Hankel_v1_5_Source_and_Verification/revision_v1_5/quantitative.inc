\subsection{Quantitative truncation with two verified Stein witnesses}
\label{sec:quantitative-truncation}
The following result supplies one explicit answer to the quantitative
truncation target. Its constants are allowed to deteriorate when the
source approaches the obstruction of Theorem~\ref{thm:separation}.
It uses full matrix inequalities, not a spectral-radius margin.

\begin{theorem}[Conditional square-root bound for balanced truncation]
\label{thm:quantitative-truncation}
Let a minimal stable realization $(A_a,b,c)$ be balanced with
$P=Q=\Sigma=\diag(\Sigma_1,\Sigma_2)\succ0$. Retain the first $k$
coordinates and put $B_a=A_{a,11}$, $g(w)=c_1^*B_wb_1$.
Assume the tuple $B$ is independently certified mean-square stable.
Suppose positive definite matrices $X,Y$ and $0\le r_X,r_Y<1$ satisfy
\[
 \sum_a A_aXA_a^*\preceq r_X^2X,\qquad
 \sum_a A_a^*YA_a\preceq r_Y^2Y.
\]
Write $\kappa_X=\|X\|/\lambda_{\min}(X)$ and similarly for $Y$.
Then
\begin{equation}
 \|H_f-H_g\|\le\sqrt{\sigma_1\tr\Sigma_2}
 \left(\frac{\sqrt{\kappa_X}}{1-r_X}
             +\frac{\sqrt{\kappa_Y}}{1-r_Y}\right).
 \label{eq:quantitative-truncation}
\end{equation}
In particular $\tr\Sigma_2\le(n-k)\sigma_{k+1}$ gives an explicit
square-root rate with all other displayed parameters fixed.
\end{theorem}
\begin{proof}
Let $E:\mathbb F^k\to\mathbb F^n$ be the coordinate inclusion and
let $C,O$ and $C_B,O_B$ denote reachability and observation maps.
The leading Stein blocks show $P_B,Q_B\preceq\Sigma_1$: for example,
iterate $\Sigma_1\succeq b_1b_1^*+\sum_aB_a\Sigma_1B_a^*$ and
use stability of $B$. Hence $\|O\|=\sqrt{\sigma_1}$ and
$\|C_B\|\le\sqrt{\sigma_1}$.

For $z_w=A_wb-EB_wb_1$, the tree recursion is
\[
 z_\varepsilon=(0,b_2)^\top,\qquad
 z_{aw}=A_az_w+(0,A_{a,21}B_wb_1)^\top.
\]
On the Hilbert sum of state vectors indexed by words, with each state
norm $\|z\|_Y$, the map $(z_w)\mapsto(A_az_w)_{aw}$ has norm at
most $r_Y$. Root and nonroot forcing coordinates are disjoint. Their
squared Euclidean sum is
\[
 \|b_2\|^2+\sum_a\tr(A_{a,21}P_BA_{a,21}^*)
 \le\tr\Sigma_2,
\]
by the lower reachability Stein block. The Neumann series on the tree
and equivalence of the state norms therefore give
\[
 \|C-EC_B\|_{\rm HS}
 \le\frac{\sqrt{\kappa_Y\tr\Sigma_2}}{1-r_Y}.
\]
Apply the identical argument to $(A_a^*,c)$ with witness $X$, reversing
word order, to obtain
$\|O-O_BE^*\|_{\rm HS}\le
\sqrt{\kappa_X\tr\Sigma_2}/(1-r_X)$. Finally,
\[
 H_f-H_g=O(C-EC_B)+(O-O_BE^*)EC_B.
\]
The operator norm is at most the displayed Hilbert--Schmidt bounds
times $\sqrt{\sigma_1}$, proving the assertion.
\end{proof}
Both witnesses exist for each stable finite tuple, but their quantitative
constants need not be uniform. A witness for the full tuple does not
automatically give a witness for its leading block. For example,
\[
 A=\begin{pmatrix}2&1\\-4&-2\end{pmatrix},\qquad A^2=0,
 \qquad X=I+AA^*,\qquad X-AXA^*=I,
\]
whereas $A_{11}=2$ is unstable. This refutes the proposed inheritance
step. Balanced truncation and singular-value automata have established
antecedents \citep{sva2019}; the theorem above records the precise
conditional operator-norm estimate used here, without a novelty claim.

\begin{theorem}[A sufficient dyadic precision chosen before recertification]
\label{thm:explicit-dyadic}
Let an $n$-state tuple $(A_a,b,c)$ have witnesses as above. Choose
$r\in(0,1)$ and $\kappa\ge1$ with
$r_X,r_Y\le r$, $\kappa_X,\kappa_Y\le\kappa$, and put $a=1-r$.
Round every scalar entry to a dyadic number with absolute error at most
$h=2^{-b_0}$. For complex entries, replace these scalar-error bounds
by bounds on their complex moduli. Define
\[
 D_b=\frac{2\sqrt{\kappa n}}a+
                  \frac{2\kappa\sqrt d\,n\|b\|}{a^2},\qquad
 D_c=\frac{2\sqrt{\kappa n}}a+
                  \frac{2\kappa\sqrt d\,n\|c\|}{a^2},
\]
\[
 K= D_c\frac{\sqrt\kappa\|b\|}a
       +\frac{2\sqrt\kappa(\|c\|+1)}a D_b.
\]
The choices
\begin{equation}
 h\le\min\left\{n^{-1/2},
       \frac{a}{2\sqrt\kappa\sqrt d\,n},\frac\epsilon K\right\}
 \label{eq:explicit-dyadic}
\end{equation}
guarantee mean-square stability of the rounded tuple and
$\|H_f-H_{\widehat f}\|\le\epsilon$.
If $K=0$, the last condition is omitted. Rational inputs and rational
upper bounds for these constants give polynomial bit cost in the input
length and $b_0$ for rounding and exact joint-Gramian recertification.
\end{theorem}
\begin{proof}
In either witness norm the tree-operator perturbation has norm at most
$\sqrt\kappa\sqrt d\,nh\le a/2$. Thus both rounded tree operators
have norm at most $1-a/2$, proving stability. The resolvent identity,
$\|\delta b\|,\|\delta c\|\le\sqrt n h$, and the original and
rounded resolvent bounds $a^{-1}$ and $2a^{-1}$ give
\[
 \|\widehat C-C\|_{\rm HS}\le D_bh,\qquad
 \|\widehat O-O\|_{\rm HS}\le D_ch.
\]
Indeed, the two terms for reachability are the changed root vector
and the changed tree operator applied to the original reachability
sequence; their bounds are exactly the two terms of $D_bh$.
The dual argument gives $D_ch$. Also
$\|C\|\le\sqrt\kappa\|b\|/a$ and
$\|\widehat O\|\le2\sqrt\kappa(\|c\|+1)/a$.
Use $\widehat H-H=(\widehat O-O)C+\widehat O(\widehat C-C)$.
This proves the bound $Kh$. Rational Stein systems and their joint
error systems have polynomial dimension; exact linear algebra and
positive-semidefinite tests have polynomial bit cost. This is a
bound for specified witnesses, not for finding a well-conditioned
balancing transform or a uniformly well-conditioned truncation.
\end{proof}
