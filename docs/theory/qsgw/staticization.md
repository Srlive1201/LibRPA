# QSGW static correlation potential without redundant continuation

The numerical input is a complete complex self-energy matrix at positive
imaginary frequencies, in the fixed reference-state basis. For each ordered
entry, let

$$
P_{ij}(z)=\operatorname{Pade}\{\Sigma^c_{ij}(i\omega_\ell)\}(z),\qquad
e_i=\epsilon_i-E_F.
$$

The same node selection, parameter counts and optional two-stage resampling
are used as in the original continuation. The imaginary-axis input is not
assumed Hermitian: $P_{ij}$ and $P_{ji}$ must be fitted independently.
At a real energy the Hermitian part is

$$
H_{ij}(E)=\frac{P_{ij}(E)+P_{ji}(E)^*}{2}.
$$

The existing two staticization modes are

$$
V^{c,A}_{ij}=\frac{H_{ij}(e_i)+H_{ij}(e_j)}2,
\qquad
V^{c,B}_{ij}=\begin{cases}H_{ii}(e_i),&i=j,\\H_{ij}(0),&i\ne j.\end{cases}
$$

This change preserves these equations. It introduces no real projection,
time-reversal, spatial-symmetry or collinearity assumption, and uses the
complete supplied state dimension for scalar, magnetic and SOC inputs.
Correctness of the upstream self-energy and symmetry conventions remains
an independent requirement.

## Pairwise evaluation

For a diagonal entry, fit $P_{ii}$ once and return
$(P_{ii}(e_i)+P_{ii}(e_i)^*)/2$. For each pair $i<j$, fit both ordered entries
once. Mode B evaluates them at zero; Mode A evaluates both at $e_i$ and
$e_j$, Hermitizes at each energy, and then averages. The lower entry is
$V^c_{ji}=V^{c*}_{ij}$, which follows from both original formulas.

Only the two fits for the current pair and the final matrix are retained.
The source matrices are borrowed as immutable pointers after their original
shape/frequency/finite-data validation. No input matrix is modified or copied
to implement the fit reuse. A constant amount of floating-point evaluation
order can change for lower-triangle Mode A entries; this is not a changed
staticization formula.

The optional resampling route still fits the original data, evaluates that
fit on the requested imaginary grid, and fits the resampled values using
the original two independent parameter counts. Each ordered entry performs
this chain once. Source, coefficient and required-target finite guards remain.
Mode B no longer evaluates unused off-diagonal state-energy values, so an
irrelevant pole at such a discarded target cannot abort that mode.

## Work and memory

Let $N$ be the number of states. Without resampling:

| Work per k point | Original Mode A | Original Mode B | Pairwise Mode A | Pairwise Mode B |
| --- | ---: | ---: | ---: | ---: |
| Padé fits | $N^3$ | $(N+1)N^2$ | $N^2$ | $N^2$ |
| Real-energy evaluations | $N^3$ | $(N+1)N^2$ | $2N^2-N$ | $N^2$ |

With resampling, there are two Padé fits per ordered entry and the first
stage's imaginary-grid evaluations. The fit-count ratios are unchanged.
The existing Thiele coefficient workspace, quadratic in the parameter count,
is unchanged.

The original implementation retained $N$ full $N\times N$ state-energy
matrices, an $O(N^3)$ auxiliary allocation. The pairwise route retains only
the $O(N^2)$ output and pairwise continuation work. The $O(n_\omega N^2)$
source data and driver-side maps over all k points are unchanged. MPI rank
ownership is also unchanged: this removes redundant work inside the existing
root stage rather than introducing a new parallel scheme.

For $N=220$ and 317 grid/path points, Mode B's fit count changes from
3,390,758,800 to 15,342,800, or 221 times fewer fits. Its removed state-energy
matrix cache alone occupies 170,368,000 bytes per k point. These are work and
allocation counts, not a measurement of whole-GW acceleration or peak memory.

## Numerical checks

The registered `test_qsgw_correlation_potential` evaluates an independent
analytic rational self-energy through spectral-projector formulas. It covers
asymmetric complex ordered entries, nonzero Fermi energy, different state
energies, both modes, direct/subsampled/two-stage fits, a full-spinor fixture,
and multiple collinear spin/k indices. Release assertions remain enabled.

Existing numerical end-to-end QSGW cases check the complete iteration traces
for molecular, periodic head-only and band/mixing/cut paths with multiple MPI
ranks. Agreement of these numerical controls establishes implementation
equivalence within their tolerances; it does not establish material convergence,
causal continuation or a physically accepted spectrum.
