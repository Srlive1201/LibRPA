# Silicon atom-exchange symmetry regression

This case reuses the dataset of `g0w0_aims_Si_libri` unchanged and enables
`use_symmetry_exx`, `use_symmetry_gw`, and `use_symmetry_rpa`.
The two equivalent Si atoms are exchanged by 24 of the 48 spatial operations.
The 3x3x3 full-grid producer data have 18 states. Three numerical checks
cover all 18 states at all 27 k points:

1. All columns must match the verified symmetry-ON reference to 1e-4 eV.
2. EXX must match the independent full-grid reference to 2e-3 eV.
3. Re/Im Sigma_c and QP energies must match that reference to 1e-3 eV.

`refs/.../librpa/full_grid_reference.out` is copied byte-for-byte from the
pre-existing `g0w0_aims_Si_libri` symmetry-OFF reference. It was not generated
by the symmetry implementation under test. No existing reference or numerical
tolerance was changed. The additional ON snapshot was generated with
e0407be49c129ddda0008246dd9be5034a80419f; one and four MPI ranks agree in every
printed table entry. The runtime must explicitly report all three symmetry
switches as true. Missing data, NaN, and infinity fail.

The full-grid fixture is not exactly symmetry covariant even before GW:
its KS eigenvalues differ by up to 3.60756 meV within a spatial k star, and
the Gamma valence triplet is split by 2.44224 meV. The audit compares the
producer `band_out` eigenvalues using `bz_sampling_out` k points and all
48 column-convention operations in `stru_out`; see
`producer_symmetry_audit.json`. Thus an OFF/ON equality of 1e-4 eV is not
assumed for this archived numerical input. This observation does not identify
every source of the small GW residual.

Controlled results with unchanged input: the corrected ON/OFF difference is
1.16 meV in EXX, 0.75 meV in Re Sigma_c, 0.39 meV in Im Sigma_c, and 0.72 meV
in QP energy. A fresh OFF run reproduces the old OFF reference exactly at
printed precision. The uncorrected 7e40c5bb source differs by 1.39267 eV in
QP energy and fails both the independent comparison and the ON snapshot.

The test runs for both one and four MPI ranks. It complements the compact
ABACUS Si k444 occupied-density CTest, which explicitly exercises restoration
from irreducible k points and validates Bloch phase/atom routing independently
of the QP solver. The student ABACUS k444 example is a separate head/wing,
shrink, reader-v1 end-to-end validation.
