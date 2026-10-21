# Padé interpolation at singular reciprocal differences

`AnalyContPade` retains the original sample selection and multiplicative
Thiele fraction whenever its construction and evaluation are regular.
A zero reciprocal difference need not mean that the data have no rational
interpolant. The checked fallback changes coordinates for the same Padé
degree class, preserving the samples and their complex amplitudes.

## The failure and the degree class

For selected samples $(z_i,f_i)$ the original construction uses

\[
g_i^{(0)}=f_i,\qquad
g_i^{(p)}=
\frac{g_{p-1}^{(p-1)}-g_i^{(p-1)}}
     {(z_i-z_{p-1})g_i^{(p-1)}},\qquad a_p=g_p^{(p)}.
\]

Repeated values can give an intermediate zero while the remaining data
are nonconstant. For example, $f(z)=1/(1+z^2)$ at
$z=(-1,1,-2,2)$ gives $(1/2,1/2,1/5,1/5)$; no ordering avoids all repeated
values. A nonconstant function with a zero sample also defeats this finite
multiplicative representation. Neither event justifies replacing a small
matrix element by zero. Complex offdiagonal entries of a matrix-causal
self-energy can have both properties.

The forward convergents obey

\[
A_{-1}=0,\ A_0=a_0,\quad B_{-1}=B_0=1,\qquad
A_j=A_{j-1}+a_j(z-z_{j-1})A_{j-2},
\]

with the same recurrence for $B_j$. An $N$-term fraction therefore has
type at most $[\lfloor(N-1)/2\rfloor/\lfloor N/2\rfloor]$.
For sixteen samples this is $[7/8]$. Adding a constant to the data and
subtracting it after interpolation can raise the numerator degree to
eight; that operation does not preserve this class.

## A representation that permits zeros

Write $m=\lfloor N/2\rfloor$ and choose $m+1$ distinct support nodes
$s_j$ from the unchanged selected data. The barycentric form is

\[
r(z)=\frac{\sum_j w_j f(s_j)/(z-s_j)}
           {\sum_j w_j/(z-s_j)}.
\]

Each nonsupport sample supplies the homogeneous equation

\[
\sum_j w_j\frac{f_i-f(s_j)}{z_i-s_j}=0.
\]

For odd $N=2m+1$, these are $m$ equations for $m+1$ weights and the type
is $[m/m]$. For even $N=2m$, there are $m-1$ such equations. The additional
constraint

\[
\sum_j w_j f(s_j)=0
\]

cancels the leading numerator coefficient and gives $[m-1/m]$.
The implementation scales nodes and values, equilibrates the small matrix,
and obtains homogeneous null vectors by complete-pivot elimination in
`complex<long double>`, without forming normal equations. The pivot screen
is `64 * epsilon(long double)` after equilibration. It is a numerical rank
screen, not an absolute-amplitude cutoff. A finite collection of basis
vectors and their sum is checked; deterministic support shifts are tried
if the first layout does not pass. Exhausting these candidates still throws.

For even $N$, select the largest-magnitude support value as anchor $L$,
and explicitly set $w_L f_L=-\sum_{j\ne L}w_j f_j$. Evaluation uses

\[
\widetilde A(z)=\sum_{j\ne L}
  \frac{w_j f_j(s_j-s_L)}{z-s_j},\qquad
\widetilde B(z)=w_L+\sum_{j\ne L}
  \frac{w_j(z-s_L)}{z-s_j},\qquad r=\widetilde A/\widetilde B.
\]

This cancellation is structural, so a rounding residual in the weight
sum does not introduce a constant high-frequency tail. The derivative is
$(\widetilde A'\widetilde B-\widetilde A\widetilde B')/\widetilde B^2$;
differentiate the displayed terms directly. At a nonzero-weight support,
use the finite limiting value and derivative. A zero-weight support must
pass the actual remaining ratio, rather than being assigned its recorded
sample automatically.

Every candidate is evaluated at **all** original selected nodes and must
be finite and satisfy the existing gate

\[
|r(z_i)-f_i|\le64\,\epsilon_{\rm double}\max_j|f_j|.
\]

The public source records keep their original selection and order. No
frequency, parameter count, fitted physical model, staticization formula,
mixing rule or symmetry assumption changes. The ordinary successful path
retains its original floating-point operations.

## Removable singularities during evaluation

A nested denominator can vanish at a regular numerator root even when
construction succeeds. For example $(z-3)/(z+1)$ has a regular value zero
and derivative $1/4$ at $z=3$. If the original nested value or derivative
is nonfinite, evaluate the same forward convergents above, together with
their differentiated recurrences. Joint rescaling of both successive
numerator/denominator pairs and their derivatives controls intermediate
growth and leaves the quotient and its derivative unchanged. A genuine
pole remains nonfinite. Finite original evaluations are preserved.

## Numerical controls and the material boundary

`test_analycont` checks repeated values with no admissible finite Thiele
ordering, raw zero samples, complex matrix-causal offdiagonal examples,
analytic derivatives, a regular numerator root, and amplitude/phase
rescalings down to $10^{-20}$. A captured sixteen-point Bi2Se3 Gamma
offdiagonal provides an additional numerical regression. Its samples zero
and four are bit-identical, causing a second-stage reciprocal-difference
zero; the function is not constant. The reference at the Mode-B target is

\[
\Sigma_{12,13}(0)=
9.0136364054620232\times10^{-15}
-1.4319544233203666\times10^{-15}i\quad\mathrm{Ha}.
\]

The regression reference comes from an independent 90-decimal polynomial
$[7/8]$ interpolation of the exact binary64 inputs. Original successful
reverse and diagonal continuations are separately compared with the old
kernel. The selected-node fit and agreement of alternative representations
test implementation correctness. They do not establish causality of the
raw material self-energy, real-axis continuation accuracy for other
entries, QSGW self-consistency, an accepted spectrum or topology.
