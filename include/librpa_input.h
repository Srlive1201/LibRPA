#pragma once
/**
 * @file librpa_input.h
 * @brief Input data APIs for LibRPA.
 */

#include <stddef.h>
#include "librpa_enums.h"
#include "librpa_handler.h"

// C APIs
#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Set meanfield wavefunction dimension.
 * @param[in] h        Handler.
 * @param[in] nspins   Number of spin channels.
 * @param[in] nkpts    Number of k-points.
 * @param[in] nstates  Number of electronic states.
 * @param[in] nbasis   Number of atomic basis functions.
 * @param[in] nspinor  Number of spin components per wavefunction
 */
void librpa_set_scf_dimension(LibrpaHandler* h, int nspins, int nkpts, int nstates, int nbasis, int nspinor = 1);

/**
 * @brief Set occupation numbers, eigenvalues, and Fermi level of meanfield input.
 * @param[in] h        Handler.
 * @param[in] nspins   Number of spin channels.
 * @param[in] nkpts    Number of k-points.
 * @param[in] nstates  Number of electronic states.
 * @param[in] wg       Occupation numbers (size: nspins * nkpts * nstates).
 * @param[in] ekb      Eigenvalues (size: same as wg).
 * @param[in] efermi   Fermi level.
 */
void librpa_set_wg_ekb_efermi(LibrpaHandler* h, int nspins, int nkpts, int nstates,
                              const double* wg, const double* ekb, double efermi);

/**
 * @brief Set wavefunction coefficients (real/imag separate arrays).
 * @param[in] h              Handler.
 * @param[in] ispin          Spin index (0-based).
 * @param[in] ik             K-point index (0-based).
 * @param[in] nstates_local  Local number of states.
 * @param[in] nbasis_local   Local number of basis functions.
 * @param[in] wfc_real       Real part of wavefunction coefficients.
 * @param[in] wfc_imag       Imaginary part of wavefunction coefficients.
 */
void librpa_set_wfc(LibrpaHandler* h, int ispin, int ik, int nstates_local, int nbasis_local,
                    const double* wfc_real, const double* wfc_imag);

/**
 * @brief Set wavefunction coefficients of spinor format (real/imag separate arrays).
 * @param[in] h              Handler.
 * @param[in] ik             K-point index (0-based).
 * @param[in] nstates_local  Local number of states.
 * @param[in] nbasis_local   Local number of basis functions.
 * @param[in] wfc_up_real    Real part of wavefunction coefficients (spin-up component).
 * @param[in] wfc_up_imag    Imaginary part of wavefunction coefficients (spin-up component).
 * @param[in] wfc_dn_real    Real part of wavefunction coefficients (spin-down component).
 * @param[in] wfc_dn_imag    Imaginary part of wavefunction coefficients (spin-down component).
 */
void librpa_set_wfc_spinor(LibrpaHandler* h, int ik, int nstates_local, int nbasis_local,
                           const double* wfc_up_real, const double* wfc_up_imag,
                           const double* wfc_dn_real, const double* wfc_dn_imag);

/**
 * @brief Set wavefunction coefficients (packed complex array).
 * @param[in] h              Handler.
 * @param[in] ispin          Spin index (0-based).
 * @param[in] ik             K-point index (0-based).
 * @param[in] nstates_local  Local number of states.
 * @param[in] nbasis_local   Local number of basis functions.
 * @param[in] wfc_ri         Complex wavefunction coefficients (packed).
 */
void librpa_set_wfc_packed(LibrpaHandler* h, int ispin, int ik, int nstates_local, int nbasis_local,
                           const double* wfc_ri);

/**
 * @brief Set wavefunction coefficients (packed complex array).
 * @param[in] h              Handler.
 * @param[in] ik             K-point index (0-based).
 * @param[in] nstates_local  Local number of states.
 * @param[in] nbasis_local   Local number of basis functions.
 * @param[in] wfc_up_ri      Complex wavefunction coefficients of spin-up component (packed).
 * @param[in] wfc_dn_ri      Complex wavefunction coefficients of spin-down component (packed).
 */
void librpa_set_wfc_spinor_packed(LibrpaHandler* h, int ik, int nstates_local,
                                  int nbasis_local, const double* wfc_up_ri,
                                  const double* wfc_dn_ri);

/**
 * @brief Set atomic orbital basis size for wavefunctions.
 * @param[in] h         Handler.
 * @param[in] natoms    Number of atoms.
 * @param[in] nbs_wfc   Number of basis functions per atom.
 * @param[in] nshells   Optional number of angular shells per atom.
 * @param[in] l_shells  Optional concatenated angular momenta, grouped by atom.
 */
void librpa_set_ao_basis_wfc(LibrpaHandler* h, int natoms, const size_t* nbs_wfc,
                             const int* nshells = nullptr, const int* l_shells = nullptr);

/**
 * @brief Set auxiliary atomic orbital basis size.
 * @param[in] h         Handler.
 * @param[in] natoms    Number of atoms.
 * @param[in] nbs_aux   Number of auxiliary basis functions per atom.
 * @param[in] nshells   Optional number of angular shells per atom.
 * @param[in] l_shells  Optional concatenated angular momenta, grouped by atom.
 */
void librpa_set_ao_basis_aux(LibrpaHandler* h, int natoms, const size_t* nbs_aux,
                             const int* nshells = nullptr, const int* l_shells = nullptr);

/**
 * @brief Set shrink auxiliary atomic orbital basis size.
 * @param[in] h               Handler.
 * @param[in] natoms          Number of atoms.
 * @param[in] nbs_aux_shrink  Number of shrink auxiliary basis functions per atom.
 * @param[in] nshells         Optional number of angular shells per atom.
 * @param[in] l_shells        Optional concatenated angular momenta, grouped by atom.
 */
void librpa_set_ao_basis_aux_shrink(LibrpaHandler* h, int natoms,
                                    const size_t* nbs_aux_shrink,
                                    const int* nshells = nullptr,
                                    const int* l_shells = nullptr);

/**
 * @brief Set basis convention metadata used by symmetry-based reductions.
 *
 * This metadata is optional for calculations that do not use symmetry.
 *
 * @param[in] h       Handler.
 * @param[in] bloch_phase Bloch-sum phase sign; must be either +1 or -1.
 * @param[in] bloch_ratom Coefficient of atom position in the Bloch-sum phase; must be -1, 0, or +1.
 * @param[in] order   Ordering of real spherical harmonics inside each shell.
 * @param[in] nega_m  Coefficient convention for the m < 0 real-spherical-harmonic branch.
 * @param[in] posi_m  Coefficient convention for the m > 0 real-spherical-harmonic branch.
 */
void librpa_set_basis_convention(LibrpaHandler* h, int bloch_phase, int bloch_ratom,
                                 LibrpaAngularOrder order, LibrpaRshCoeff nega_m,
                                 LibrpaRshCoeff posi_m);

/**
 * @brief Set real-space symmetry operations.
 *
 * @param[in] h          Handler.
 * @param[in] n_symops   Number of symmetry operations.
 * @param[in] row_conv   Positive if rotations use the row-fractional convention.
 * @param[in] rotmats    Flattened rotation matrices, length 9*n_symops.
 * @param[in] trans      Optional flattened fractional translations, length 3*n_symops.
 */
void librpa_set_symmetry_operations(LibrpaHandler* h, int n_symops, int row_conv,
                                    const int* rotmats, const double* trans);

/**
 * @brief Set symmetry operations with explicit spin actions (magnetic or
 *        spin-space group input).
 *
 * The spatial part is normalized exactly like librpa_set_symmetry_operations.
 * Each operation additionally carries an antiunitary flag and an SU(2) spin
 * rotation, normalized internally to SymmetrySpinOperation entries whose
 * spatial_id indexes the operation in input order.
 *
 * @param[in] h          Handler.
 * @param[in] n_ops      Number of symmetry operations.
 * @param[in] row_conv   Positive if rotations use the row-fractional convention.
 * @param[in] rotmats    Flattened rotation matrices, length 9*n_ops.
 * @param[in] trans      Optional flattened fractional translations, length 3*n_ops.
 * @param[in] antiunitary Optional flags (0/1) per operation, length n_ops;
 *                        null means all unitary.
 * @param[in] spin_u     Optional flattened SU(2) matrices, 8 doubles per
 *                       operation as (re,im) x 4 in row-major 2x2 order.
 * @param[in] spin_action_source 0 = Identity (null spin_u means the identity),
 *                       1 = ExplicitSpinSpace (spin_u required, unitarity
 *                       checked), 2 = DerivedFromSpatialSOC (null spin_u means
 *                       LibRPA reconstructs U[det(Q)Q] from the spatial part,
 *                       requiring the lattice; a provided spin_u is validated
 *                       and adopted as the upstream authority).
 * @param[in] grey_group Nonzero to declare the table as the unitary block of a
 *                       grey group: LibRPA appends one antiunitary copy per
 *                       operation (unitary block first).
 */
void librpa_set_symmetry_spin_operations(
    LibrpaHandler* h, int n_ops, int row_conv,
    const int* rotmats, const double* trans,
    const int* antiunitary, const double* spin_u,
    int spin_action_source, int grey_group);

/**
 * @brief Set direct and reciprocal lattice vectors.
 * @param[in] h         Handler.
 * @param[in] lat_mat   Lattice vectors (3x3, column-major, in Bohr).
 * @param[in] G_mat     Reciprocal lattice vectors (3x3, column-major, in Bohr^-1).
 */
void librpa_set_latvec_and_G(LibrpaHandler* h, const double lat_mat[9], const double G_mat[9]);

/**
 * @brief Set atom types and Cartesian coordinates.
 * @param[in] h          Handler.
 * @param[in] natoms     Number of atoms.
 * @param[in] types      Atom types (species indices, 0-based).
 * @param[in] posi_cart  Cartesian coordinates (3*natoms, in Bohr).
 */
void librpa_set_atoms(LibrpaHandler* h, int natoms, const int* types, const double* posi_cart);

/**
 * @brief Set k-point grid vectors.
 *
 * @param[in] h      Handler.
 * @param[in] nk1    Number of k-points along direction 1.
 * @param[in] nk2    Number of k-points along direction 2.
 * @param[in] nk3    Number of k-points along direction 3.
 * @param[in] nkpts     Number of loaded SCF k-points.
 * @param[in] kvecs     K-point vectors (3*nkpts, Cartesian coordinates, in Bohr^-1).
 * @param[in] kweights  Optional k-point weights, normalized internally to sum to one.
 *                      Pass NULL if unavailable.
 */
void librpa_set_kgrids_kvec(LibrpaHandler* h, int nk1, int nk2, int nk3,
                            int nkpts, const double* kvecs, const double* kweights);

/**
 * @brief Set the mapping from loaded SCF k-points to Coulomb q-points.
 *
 * @param[in] h          Handler.
 * @param[in] nkpts      Number of loaded SCF k-points.
 * @param[in] map_q_ks   Mapping from each loaded SCF k-point to the representative q-point,
 *                       indexed in the loaded SCF k-list (0-based).
 *
 * Example: four loaded k-points where the first two and last are Coulomb q-points,
 *         and the third point maps to the second q-point, then map_q_ks is (0, 1, 1, 3).
 */
void librpa_set_kq_mapping(LibrpaHandler* h, int nkpts, const int* map_q_ks);

/**
 * @brief Set an explicit full-grid to loaded SCF representative mapping.
 *
 * Local operation without MPI collectives. Coordinates are Cartesian in Bohr^-1.
 * full_to_scf uses 0-based indices into the loaded SCF list. There is one Coulomb
 * representative per loaded SCF point. Does not change the SCF storage dimension.
 */
void librpa_set_irreducible_kgrids_kvec(LibrpaHandler* h, int nk1, int nk2, int nk3,
                                       int n_scf, const double* kvecs_scf,
                                       int n_full, const double* kvecs_full,
                                       const int* full_to_scf);

/**
 * @brief Set local RI coefficients.
 *
 * @param[in] h          Handler.
 * @param[in] routing    Parallel routing strategy.
 * @param[in] I          Atom I index (0-based).
 * @param[in] J          Atom J index (0-based).
 * @param[in] nbasis_i   Number of basis functions on atom I.
 * @param[in] nbasis_j   Number of basis functions on atom J.
 * @param[in] naux_mu    Number of auxiliary functions for mu index.
 * @param[in] R          Lattice vector [R1, R2, R3].
 * @param[in] Cs_in      RI triple coefficients.
 * @param[in] shrink_aux If positive, parse coefficients to the shrink auxiliary basis.
 */
void librpa_set_lri_coeff(LibrpaHandler* h, LibrpaParallelRouting routing, int I, int J,
                          int nbasis_i, int nbasis_j, int naux_mu, const int R[3],
                          const double* Cs_in, int shrink_aux = 0);

/**
 * @brief Set bare Coulomb matrix elements (atom-pair format).
 *
 * @param[in] h               Handler.
 * @param[in] ik              K-point index.
 * @param[in] I               Atom I index.
 * @param[in] J               Atom J index.
 * @param[in] naux_mu         Number of auxiliary functions for mu.
 * @param[in] naux_nu         Number of auxiliary functions for nu.
 * @param[in] Vq_real_in      Real part of Coulomb matrix.
 * @param[in] Vq_imag_in      Imaginary part of Coulomb matrix.
 * @param[in] vq_threshold    Threshold for screening.
 */
void librpa_set_aux_bare_coulomb_k_atom_pair(LibrpaHandler* h, int ik, int I, int J, int naux_mu,
                                             int naux_nu, const double* Vq_real_in,
                                             const double* Vq_imag_in, double vq_threshold);

/**
 * @brief Set bare Coulomb matrix elements (atom-pair format, packed complex).
 *
 * @param[in] h               Handler.
 * @param[in] ik              K-point index.
 * @param[in] I               Atom I index.
 * @param[in] J               Atom J index.
 * @param[in] naux_mu         Number of auxiliary functions for mu.
 * @param[in] naux_nu         Number of auxiliary functions for nu.
 * @param[in] Vq_ri_in        Coulomb matrix with interleaved real/imaginary parts.
 * @param[in] vq_threshold    Threshold for screening.
 */
void librpa_set_aux_bare_coulomb_k_atom_pair_packed(LibrpaHandler* h, int ik, int I, int J,
                                                    int naux_mu, int naux_nu,
                                                    const double* Vq_ri_in,
                                                    double vq_threshold);

/**
 * @brief Set truncated Coulomb matrix elements (atom-pair format).
 *
 * @param[in] h               Handler.
 * @param[in] ik              K-point index.
 * @param[in] I               Atom I index.
 * @param[in] J               Atom J index.
 * @param[in] naux_mu         Number of auxiliary functions for mu.
 * @param[in] naux_nu         Number of auxiliary functions for nu.
 * @param[in] Vq_real_in      Real part of Coulomb matrix.
 * @param[in] Vq_imag_in      Imaginary part of Coulomb matrix.
 * @param[in] vq_threshold    Threshold for screening.
 */
void librpa_set_aux_cut_coulomb_k_atom_pair(LibrpaHandler* h, int ik, int I, int J, int naux_mu,
                                            int naux_nu, const double* Vq_real_in,
                                            const double* Vq_imag_in, double vq_threshold);

/**
 * @brief Set truncated Coulomb matrix elements (atom-pair format, packed complex).
 *
 * @param[in] h               Handler.
 * @param[in] ik              K-point index.
 * @param[in] I               Atom I index.
 * @param[in] J               Atom J index.
 * @param[in] naux_mu         Number of auxiliary functions for mu.
 * @param[in] naux_nu         Number of auxiliary functions for nu.
 * @param[in] Vq_ri_in        Coulomb matrix with interleaved real/imaginary parts.
 * @param[in] vq_threshold    Threshold for screening.
 */
void librpa_set_aux_cut_coulomb_k_atom_pair_packed(LibrpaHandler* h, int ik, int I, int J,
                                                   int naux_mu, int naux_nu,
                                                   const double* Vq_ri_in,
                                                   double vq_threshold);

/**
 * @brief Set bare Coulomb matrix elements (2D block format).
 *
 * @param[in] h            Handler.
 * @param[in] ik           K-point index.
 * @param[in] mu_begin     Starting index for mu.
 * @param[in] mu_end       Ending index for mu.
 * @param[in] nu_begin     Starting index for nu.
 * @param[in] nu_end       Ending index for nu.
 * @param[in] Vq_real_in   Real part of Coulomb matrix.
 * @param[in] Vq_imag_in   Imaginary part of Coulomb matrix.
 */
void librpa_set_aux_bare_coulomb_k_2d_block(LibrpaHandler* h, int ik, int mu_begin, int mu_end,
                                            int nu_begin, int nu_end, const double* Vq_real_in,
                                            const double* Vq_imag_in);

/**
 * @brief Set bare Coulomb matrix elements (2D block format, packed complex).
 *
 * @param[in] h            Handler.
 * @param[in] ik           K-point index.
 * @param[in] mu_begin     Starting index for mu.
 * @param[in] mu_end       Ending index for mu.
 * @param[in] nu_begin     Starting index for nu.
 * @param[in] nu_end       Ending index for nu.
 * @param[in] Vq_ri_in     Coulomb matrix with interleaved real/imaginary parts.
 */
void librpa_set_aux_bare_coulomb_k_2d_block_packed(LibrpaHandler* h, int ik, int mu_begin,
                                                   int mu_end, int nu_begin, int nu_end,
                                                   const double* Vq_ri_in);

/**
 * @brief Set truncated Coulomb matrix elements (2D block format).
 *
 * @param[in] h            Handler.
 * @param[in] ik           K-point index.
 * @param[in] mu_begin     Starting index for mu.
 * @param[in] mu_end       Ending index for mu.
 * @param[in] nu_begin     Starting index for nu.
 * @param[in] nu_end       Ending index for nu.
 * @param[in] Vq_real_in   Real part of Coulomb matrix.
 * @param[in] Vq_imag_in   Imaginary part of Coulomb matrix.
 */
void librpa_set_aux_cut_coulomb_k_2d_block(LibrpaHandler* h, int ik, int mu_begin, int mu_end,
                                           int nu_begin, int nu_end, const double* Vq_real_in,
                                           const double* Vq_imag_in);

/**
 * @brief Set truncated Coulomb matrix elements (2D block format, packed complex).
 *
 * @param[in] h            Handler.
 * @param[in] ik           K-point index.
 * @param[in] mu_begin     Starting index for mu.
 * @param[in] mu_end       Ending index for mu.
 * @param[in] nu_begin     Starting index for nu.
 * @param[in] nu_end       Ending index for nu.
 * @param[in] Vq_ri_in     Coulomb matrix with interleaved real/imaginary parts.
 */
void librpa_set_aux_cut_coulomb_k_2d_block_packed(LibrpaHandler* h, int ik, int mu_begin,
                                                  int mu_end, int nu_begin, int nu_end,
                                                  const double* Vq_ri_in);

/**
 * @brief Set dielectric function on imaginary frequency axis.
 *
 * @param[in] h            Handler.
 * @param[in] nfreq        Number of frequency points.
 * @param[in] omegas_imag  Imaginary frequency values.
 * @param[in] dielect_func Dielectric function values.
 */
void librpa_set_dielect_func_imagfreq(LibrpaHandler* h, int nfreq, const double* omegas_imag,
                                      const double* dielect_func);

/**
 * @brief Set velocity/momentum matrix.
 *
 * The flattened matrix order is [spin][k][cartesian][row band][column band].
 * Values are currently used by LibRPA's internal analytic head/wing convention; code-specific
 * conversions from ABACUS or FHI-aims file formats are intentionally not applied here.
 *
 * @param[in] h              Handler.
 * @param[in] n_spins        Number of spin channels.
 * @param[in] n_kpts         Number of k-points used by analytic head/wing.
 * @param[in] n_states       Number of bands.
 * @param[in] velocity_real  Real parts of the flattened velocity matrices.
 * @param[in] velocity_imag  Imaginary parts of the flattened velocity matrices.
 */
void librpa_set_velocity_matrix(LibrpaHandler* h, int n_spins, int n_kpts, int n_states,
                                const double* velocity_real, const double* velocity_imag);

/**
 * @brief Set velocity/momentum matrix.
 *
 * The flattened matrix order is [spin][k][cartesian][row band][column band], with
 * interleaved real/imaginary values.
 *
 * @param[in] h           Handler.
 * @param[in] n_spins     Number of spin channels.
 * @param[in] n_kpts      Number of k-points used by analytic head/wing.
 * @param[in] n_states    Number of bands.
 * @param[in] velocity_ri Flattened matrix with interleaved real/imaginary parts.
 */
void librpa_set_velocity_matrix_packed(LibrpaHandler* h, int n_spins, int n_kpts, int n_states,
                                       const double* velocity_ri);

/**
 * @brief Set k-points for band structure calculations.
 *
 * @param[in] h              Handler.
 * @param[in] n_kpts_band    Number of band k-points.
 * @param[in] kfrac_band     Band k-point coordinates (fractional).
 */
void librpa_set_band_kvec(LibrpaHandler* h, int n_kpts_band, const double* kfrac_band);

/**
 * @brief Set occupation numbers and eigenvalues for band k-points.
 *
 * @param[in] h            Handler.
 * @param[in] n_spins      Number of spin channels.
 * @param[in] n_kpts_band  Number of band k-points.
 * @param[in] n_states     Number of states.
 * @param[in] occ          Occupation numbers.
 * @param[in] eig          Eigenvalues.
 */
void librpa_set_band_occ_eigval(LibrpaHandler* h, int n_spins, int n_kpts_band, int n_states,
                                const double* occ, const double* eig);

/**
 * @brief Set wavefunction for band k-point (separate real/imag).
 *
 * @param[in] h              Handler.
 * @param[in] ispin          Spin index.
 * @param[in] ik_band        Band k-point index.
 * @param[in] nstates_local  Local number of states.
 * @param[in] nbasis_local   Local number of basis functions.
 * @param[in] wfc_real       Real part of wavefunction.
 * @param[in] wfc_imag       Imaginary part of wavefunction.
 */
void librpa_set_wfc_band(LibrpaHandler* h, int ispin, int ik_band, int nstates_local,
                         int nbasis_local, const double* wfc_real, const double* wfc_imag);

/**
 * @brief Set wavefunction for band k-point (spinor format, separate real/imag).
 * @param[in] h              Handler.
 * @param[in] ik_band        Band k-point index.
 * @param[in] nstates_local  Local number of states.
 * @param[in] nbasis_local   Local number of atomic basis functions.
 * @param[in] wfc_up_real    Real part of wavefunction coefficients (spin-up component).
 * @param[in] wfc_up_imag    Imaginary part of wavefunction coefficients (spin-up component).
 * @param[in] wfc_dn_real    Real part of wavefunction coefficients (spin-down component).
 * @param[in] wfc_dn_imag    Imaginary part of wavefunction coefficients (spin-down component).
 */
void librpa_set_wfc_band_spinor(LibrpaHandler* h, int ik_band, int nstates_local,
                                int nbasis_local, const double* wfc_up_real, const double* wfc_up_imag,
                                const double* wfc_dn_real, const double* wfc_dn_imag);

/**
 * @brief Set wavefunction for band k-point (packed complex).
 *
 * @param[in] h              Handler.
 * @param[in] ispin          Spin index.
 * @param[in] ik_band        Band k-point index.
 * @param[in] nstates_local  Local number of states.
 * @param[in] nbasis_local   Local number of atomic basis functions.
 * @param[in] wfc_ri         Complex wavefunction (packed).
 */
void librpa_set_wfc_band_packed(LibrpaHandler* h, int ispin, int ik_band, int nstates_local,
                                int nbasis_local, const double* wfc_ri);

/**
 * @brief Set wavefunction for band k-point (packed complex).
 *
 * @param[in] h              Handler.
 * @param[in] ik_band        Band k-point index.
 * @param[in] nstates_local  Local number of states.
 * @param[in] nbasis_local   Local number of basis functions.
 * @param[in] wfc_up_ri      Complex wavefunction of spin-up component (packed).
 * @param[in] wfc_dn_ri      Complex wavefunction of spin-down component (packed).
 */
void librpa_set_wfc_band_spinor_packed(LibrpaHandler* h, int ik_band, int nstates_local,
                                       int nbasis_local, const double* wfc_up_ri, const double* wfc_dn_ri);

/**
 * @brief Reset band structure input data.
 *
 * @param[in] h Handler.
 */
void librpa_reset_band_data(LibrpaHandler* h);

/** Query the numerical dimensions needed to prepare a QSGW task's host input.
 * band_path selects band input when nonzero. This call is noncollective.
 */
void librpa_get_meanfield_dimensions(const LibrpaHandler* h, int band_path, int* nspins,
                                     int* nkpts, int* nstates, int* naos, int* nspinor);

#ifdef __cplusplus
}
#endif
