# Si k666 QSGW band validation

This manual case exercises the full QSGW band driver on a real, scalar PBE Si
dataset: 216 uniform k points, 201 band-path points, 44 bands, eight electrons.
It is intentionally outside the small default regression suite. The large RI
and Coulomb data must be supplied separately.

The input requests two fixed-basis updates, six minimax frequencies, analytic
same-grid head correction, no wing or Hartree update, no Hamiltonian mixing,
and no symmetry acceleration. Wigner–Seitz multiple-image BvK mapping
(`option_bvk_remap=1`) is selected explicitly in both versions. The Hamiltonian cut retains four occupied plus
ten unoccupied states (mode 2, zero shift). Shrunk auxiliary bases are used.
These settings test the original PR's functionality; they do not establish
frequency, basis, k-mesh, or self-consistency convergence.

Prepare a new directory with `python prepare.py DATASET NEW_DIRECTORY`, then
run `chi0_main.exe` with MPI from `NEW_DIRECTORY/librpa` on compute nodes.
The numerical inputs are read-only links to the supplied dataset; only the
three metadata manifests are copied and converted. Use `--legacy` to stage
the original PR's hashed manifests for a before/after comparison.

Compare both iterations with the original PR executable using the same
dataset, MPI layout and input parameters. Required results are the grid and
band eigenvalue trajectories, iteration summary, and both
`QSGW_band_spin_1_1.dat` and `QSGW_band_spin_1_2.dat` tables. Compare eigenvalues
to 2e-5 Ha and printed band energies to 5e-4 eV; record the measured maxima,
path and mesh gaps, electron count, and convergence flag separately.
Cross-check common grid/path k points in the initial and updated spectra and
plot KS plus each QSGW update. A completed run alone is not a numerical pass.

The repository contains no new k666 reference results until these calculations
have completed. Run-specific dataset locations, source/binary identities,
scheduler scripts and results belong in the workspace validation evidence.
