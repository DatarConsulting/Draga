\subsection{A second truncation certificate from positive Stein resolvents}
\label{sec:resolvent-truncation}
The extended proposal supplies a useful alternative to the tree estimate
above. Its positive-resolvent argument is valid after correcting the
block inequalities. This subsection retains that argument and exploits
the commuting reference Gramians; it does not supply a uniform AAK factor.

\begin{lemma}[Hermitian norm of a positive Stein resolvent]
\label{lem:positive-resolvent}
Let $\Phi(Z)=\sum_a T_aZT_a^*$ have spectral radius less than one,
and put $\mathcal R=(I-\Phi)^{-1}$ and $X=\mathcal R(I)$.
For every Hermitian $M$,
\[
 -\|M\|X\preceq\mathcal R(M)\preceq\|M\|X,
 \qquad \|\mathcal R(M)\|\le\|X\|\,\|M\|.
\]
Thus the induced spectral norm on the Hermitian subspace is exactly
$\|X\|$, with equality at $M=I$.
\end{lemma}
\begin{proof}
Apply the positive map $\sum_{j\ge0}\Phi^j$ to
$-\|M\|I\preceq M\preceq\|M\|I$. The order bounds imply the norm
bound, and $M=I$ attains it. This is the elementary Hermitian case of
the classical norm-at-identity property for positive maps
\citep{bourin2019}, not a new general operator-theory result.
\end{proof}

\begin{theorem}[Commuting-reference certificate for balanced truncation]
\label{thm:reference-truncation}
Let $(A_a,b,c)$ be minimal, mean-square stable, and balanced with
$P=Q=\Sigma=\diag(\Sigma_1,\Sigma_2)\succ0$, where $1\le k<n$
coordinates are retained. Set $B_a=A_{a,11}$ and
$g(w)=c_1^*B_wb_1$, and independently assume $B$ is mean-square
stable. Write $\delta=\|\Sigma_2\|=\sigma_{k+1}$ and
$\sigma_1=\|\Sigma\|$. Define the four Stein solutions
\[
 X_A-\sum_aA_aX_AA_a^*=I_n,\quad
 X_B-\sum_aB_aX_BB_a^*=I_k,
\]
\[
 Y_A-\sum_aA_a^*Y_AA_a=I_n,\quad
 Y_B-\sum_aB_a^*Y_BB_a=I_k,
\]
and $K_P=\max(\|X_A\|,\|X_B\|)$,
$K_Q=\max(\|Y_A\|,\|Y_B\|)$.
On the joint space with block sizes $(k,n-k,k)$ put
\[
 T_a=\diag(A_a,B_a),\quad \beta=(b_1,b_2,b_1)^\top,
 \quad\chi=(c_1,c_2,-c_1)^\top,
\]
\[
 P_0=\begin{pmatrix}\Sigma_1&0&\Sigma_1\\0&\Sigma_2&0\\
 \Sigma_1&0&\Sigma_1\end{pmatrix},\qquad
 Q_0=\begin{pmatrix}\Sigma_1&0&-\Sigma_1\\0&\Sigma_2&0\\
 -\Sigma_1&0&\Sigma_1\end{pmatrix}.
\]
Let
\[
 D_P=\beta\beta^*+\sum_aT_aP_0T_a^*-P_0,\qquad
 D_Q=\chi\chi^*+\sum_aT_a^*Q_0T_a-Q_0,
\]
and take any $e_P\ge K_P\|D_P\|$, $e_Q\ge K_Q\|D_Q\|$.
Then
\begin{equation}
 \|H_f-H_g\|^2\le
 \max\{(\delta+e_P)(\delta+e_Q),\,
                    2\sigma_1\max(e_P,e_Q)+e_Pe_Q\}.
 \label{eq:reference-truncation}
\end{equation}
In particular, define
\[
 M_P=\Big\|\sum_aA_{a,12}A_{a,12}^*\Big\|,\quad
 M_Q=\Big\|\sum_aA_{a,21}^*A_{a,21}\Big\|,
\]
\[
 c_P=K_P(2M_P+\sqrt{M_P}),\qquad
 c_Q=K_Q(2M_Q+\sqrt{M_Q}).
\]
One may take $e_P=c_P\delta$, $e_Q=c_Q\delta$, giving
\begin{equation}
 \|H_f-H_g\|^2\le
 \max\{\delta^2(1+c_P)(1+c_Q),\,
             2\sigma_1\delta\max(c_P,c_Q)+c_Pc_Q\delta^2\}.
 \label{eq:reference-square-root}
\end{equation}
This is an $O(\sqrt\delta)$ bound when the other displayed quantities
are bounded uniformly. It contains no explicit tail-dimension factor,
but its coupling and resolvent constants may depend on that dimension.
\end{theorem}
\begin{proof}
Let $P_E,Q_E$ be the exact joint error Gramians. Subtracting the
reference Stein equations gives
$P_E-P_0=(I-\Phi_T)^{-1}D_P$, and the dual identity for $Q_E$.
The identity forcings have solutions
$\diag(X_A,X_B)$ and $\diag(Y_A,Y_B)$, respectively.
Lemma~\ref{lem:positive-resolvent} therefore gives
\[
 0\preceq P_E\preceq P_0+e_PI=:\overline P,\qquad
 0\preceq Q_E\preceq Q_0+e_QI=:\overline Q.
\]
For positive semidefinite matrices the quantity
$\lambda_{\max}(P^{1/2}QP^{1/2})$ is separately monotone in $P,Q$:
use congruence for $Q$, and the equality of the nonzero spectra of
$P^{1/2}QP^{1/2}$ and $Q^{1/2}PQ^{1/2}$ for $P$.
Apply Theorem~\ref{thm:exact}. The reference matrices commute.
On the retained symmetric and antisymmetric subspaces their pairs
of eigenvalues are $(2\sigma_j,0)$ and $(0,2\sigma_j)$; on a tail
coordinate they are $(\sigma_j,\sigma_j)$.
The largest eigenvalue of $\overline P\overline Q$ is exactly the
right-hand side of~\eqref{eq:reference-truncation}.

For the explicit estimates put
$S=\sum_aA_{a,12}\Sigma_2A_{a,12}^*$ and
$U=\sum_aA_{a,22}\Sigma_2A_{a,12}^*$.
Direct multiplication of the balanced Stein equations gives
\[
 D_P=-\begin{pmatrix}0&0&S\\0&0&U\\S&U^*&S\end{pmatrix}.
\]
Here $\|S\|\le M_P\delta$, and operator Cauchy--Schwarz, together
with $\sum_aA_{a,22}\Sigma_2A_{a,22}^*\preceq\Sigma_2$, gives
$\|U\|\le\sqrt{M_P}\delta$.
Split $D_P$ into its $S$ and $U$ parts; their norms are at most
$2\|S\|$ and $\|U\|$. Hence
$\|D_P\|\le(2M_P+\sqrt{M_P})\delta$.
Apply the same calculation to $A_a^*,c$, followed by changing the
sign of the last joint block, to obtain the asserted $D_Q$ bound.
\end{proof}

\paragraph{Correct balanced block inequalities.}
The identities imply
\[
 b_2b_2^*+\sum_aA_{a,22}\Sigma_2A_{a,22}^*\preceq\Sigma_2,
 \qquad
 \Big\|\sum_aA_{a,12}\Sigma_2A_{a,12}^*\Big\|
                  \le\min\{\sigma_1,M_P\delta\}.
\]
They do not permit replacing the norm of a positive sum by the sum
of its spectral norms, nor multiplying the two independent bounds.
For a concrete balanced counterexample let
\[
 A_0=E_{12},\quad A_1=E_{23},\quad
 b=\sqrt t\,e_3,\quad c=\sqrt t\,e_1,\qquad 0<t<1,
\]
and retain $k=1$. All products of three letters vanish, and the
reachability and observation Gramians are both $tI_3$; in particular
the representation is stable and minimal. The retained tuple is zero.
Nevertheless,
\[
 \|b_2\|^2+\sum_a\|A_{a,22}\Sigma_2^{1/2}\|^2=2t>t,
 \qquad
 \Big\|\sum_aA_{a,12}\Sigma_2A_{a,12}^*\Big\|=t>t^2.
\]
The reference-Gramian theorem uses neither false strengthening.
