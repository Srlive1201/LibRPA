# H2O QSGW Gamma band regression

This case extends the existing FHI-aims H2O dataset with a single band-path Gamma
point equal to its SCF point. The RI/Coulomb data, reference eigenvalues and
full Vxc matrix are unchanged. `band_kpath_info` contains `24 24 1 1` and Gamma;
the band eigenvalue file adds a leading spin index to each `band_out` state row.
The band wavefunction file holds the same complex coefficients as the SCF text
file, reordered from AO/state to state/AO and written as binary complex doubles.
`band_vxc_mat_spin_1_k_00001.csc` is a byte-identical copy of the SCF Vxc
matrix, using the native band filename.

Two updates use linear mixing with beta=0.2 and cut mode 2, retaining five
occupied plus ten unoccupied states with zero shift. `nfreq` is omitted to
exercise the default of 16. This covers the grid/path fixed-basis projection,
synchronized linear update, exact cut and printed band tables.

References were calculated with the pre-cleanup QSGW implementation and the
16-frequency default, using Intel oneAPI, LibRI, four MPI ranks and one thread.
Both updated Gamma spectra agree between grid and path within 3.5e-13 eV.
The new iteration-summary reference omits obsolete metadata and mixing-history
fields; its physical numerical columns are taken unchanged from that run.
The eigenvalue tolerance is 2e-5 Ha (5.4422773e-4 eV), and printed band energies
use 5e-4 eV. Existing H2O and Si reference files are preserved verbatim and the
common table comparator selects their numerical rows.

This is a short workflow regression, not a molecular convergence result.
