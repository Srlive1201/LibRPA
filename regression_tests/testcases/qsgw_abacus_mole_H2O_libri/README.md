# H2O QSGW regression with ABACUS input

This molecular Gamma-point case uses the retained ABACUS H2O dataset: 23 states,
23 atomic orbitals and eight valence electrons. The native `vxck1_nao.txt`
contains the full Vxc in the initial KS-state basis and in Ry, despite the
`_nao` suffix. The driver reads it directly and converts to Hartree.

One QSGW update uses linear mixing with beta=0.2, 16 minimax frequencies and
LibRI, with head correction and symmetry acceleration disabled. References
were calculated with the implementation immediately before removal of the
Vxc manifest reader, using Intel oneAPI, four MPI ranks and one thread.
The numerical input files and reference rows for initialization and the first
update are copied unchanged from that baseline run.

The existing table comparator checks the updated eigenvalue spectrum,
Hamiltonian residuals, Fermi energies, gaps and electron count. The eigenvalue
tolerance is 2e-5 Ha (5.4422773e-4 eV). This is a regression of the full
calculation through native ABACUS input. The second update is excluded because
repeat runs of the pre-refactor executable differed by 0.0017 eV in that step,
exceeding the existing tolerance; the first-update comparison agreed within
2.4e-10 eV. This case does not establish self-consistency convergence or
reproducibility of later updates. Other QSGW cases cover two-update trajectories.
