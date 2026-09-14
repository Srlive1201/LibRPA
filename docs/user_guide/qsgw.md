# Experimental fixed-basis QSGW

The `qsgw` and `qsgw_band` driver tasks construct self-energies in the immutable
initial KS basis and update the eigenvalues and wavefunctions each iteration.
The effective Hamiltonian replaces the initial DFT exchange-correlation
potential with exchange and a Hermitian correlation potential. Hartree updates
are not implemented. Only same-grid analytic head correction is supported;
independent-grid head updates and wing updates are rejected.

QSGW reuses the G0W0 input readers and EXX/Sigma kernels. Task-specific input
and output live in `driver/qsgw`; numerical operations in `src/qsgw` accept
matrices and mean-field objects.

## Vxc input

QSGW needs the full DFT Vxc matrix, including off-diagonal elements; the diagonal
Vxc tables used by G0W0 are insufficient. The driver reads the full matrices
directly from `input_dir`, using spin and k-point filenames with these defaults:

| Producer | SCF example (spin 1, k point 1) | Band example (spin 1, k point 1) |
| --- | --- | --- |
| FHI-aims | `xc_matr_spin_1_kpt_000001.csc` | `band_vxc_mat_spin_1_k_00001.csc` |
| ABACUS, one spin channel | `vxck1_nao.txt` | `band_vxck1_nao.txt` |
| ABACUS, two spin channels | `vxck1s1_nao.txt` | `band_vxck1s1_nao.txt` |

No additional file list or manifest is required. Optional `prefix_vxc_scf` and
`prefix_vxc_band` override the prefixes shown above (`xc_matr`/`band_vxc_mat`
for FHI-aims and `vxc`/`band_vxc` for ABACUS). Prefixes are relative to
`input_dir` unless absolute. For a separate ABACUS band output directory,
`prefix_vxc_band = OUT.band/vxc` reads its native `vxck1_nao.txt`, etc., without
renaming; the default `band_vxc` prefix distinguishes files staged alongside
SCF input. Spin and k indices are one-based and must follow
the same ordering as the corresponding SCF or band reference. Vxc files do not
provide an independent k-coordinate check.

Use `constants_choice = aims` for FHI-aims ELSI matrices in Hartree, and
`constants_choice = internal` for ABACUS complex text matrices in Ry, which are
converted to Hartree. The default `qsgw_vxc_basis = state` reads matrices in the
initial KS-state basis. This includes standard ABACUS `out_mat_xc` output,
despite its `_nao` filename suffix. Set `qsgw_vxc_basis = nao` only for genuine
ABACUS AO matrices; the driver then projects them with the reference
wavefunctions. This setting applies to both SCF and band matrices. Matrix
dimensions, finiteness and Hermiticity are checked on input.

## Iteration and output

`qsgw_mixer = linear` applies `H_next = H_in + beta * (H_out - H_in)` with
`beta = qsgw_mixing_beta`; `none` accepts the new Hamiltonian directly.
The grid and band path use the same fraction. The band Hamiltonian cut is
applied before and after mixing, so excluded states obey the cut exactly.

`qsgw_iterations.dat` records the iteration number, largest eigenvalue change
(eV), Hamiltonian residual L2 and maximum norms (Ha), Fermi energy and gap (eV),
electron count, mixing fraction and convergence flag. Residuals are measured
before mixing on the SCF grid; the maximum is the largest complex-entry
magnitude. Stopping uses the largest eigenvalue change after `qsgw_min_iter`.
`qsgw_eigenvalues.dat` contains every grid/path eigenvalue at initialization and
after each update. Optional matrix diagnostics are enabled by
`qsgw_write_iteration_matrices`.

## Numerical checks

The small regression suite covers two updates of molecular H2O from FHI-aims
with linear mixing, Si k333 with analytic head correction, and
FHI-aims H2O through `qsgw_band` with mixing and Hamiltonian truncation.
An ABACUS H2O case checks one full update through the native complex Vxc input
and Ry conversion against results from before the input simplification.
The band case uses the default 16 minimax
frequencies and the same Gamma reference on grid and path. Tests compare actual
iteration energies, Hamiltonian residuals, electron count and eigenvalue
trajectories with the existing numeric table comparator used by G0W0.

A manual Si k666 band case is under
`regression_tests/manual/Si_k666_qsgw_band`. These cases check numerical workflows;
they do not establish convergence of material predictions.
