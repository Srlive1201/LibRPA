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
 * @brief Set SCF wavefunction dimension.
 * @param[in] h        Handler.
 * @param[in] nspins   Number of spin channels.
 * @param[in] nkpts    Number of k-points.
 * @param[in] nstates  Number of electronic states.
 * @param[in] nbasis   Number of basis functions.
 */
void librpa_set_scf_dimension(LibrpaHandler* h, int nspins, int nkpts, int nstates, int nbasis);

/**
 * @brief Set occupation numbers, eigenvalues, and Fermi level.
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
 * @brief Set atomic orbital basis size for wavefunctions.
 * @param[in] h        Handler.
 * @param[in] natoms   Number of atoms.
 * @param[in] nbs_wfc  Number of basis functions per atom.
 */
void librpa_set_ao_basis_wfc(LibrpaHandler* h, int natoms, const size_t* nbs_wfc);

/**
 * @brief Set auxiliary atomic orbital basis size.
 * @param[in] h       Handler.
 * @param[in] natoms  Number of atoms.
 * @param[in] nbs_aux Number of auxiliary basis functions per atom.
 */
void librpa_set_ao_basis_aux(LibrpaHandler* h, int natoms, const size_t* nbs_aux);

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
 * @param[in] pos_cart   Cartesian coordinates (3*natoms, in Bohr).
 */
void librpa_set_atoms(LibrpaHandler* h, int natoms, const int* types, const double* posi_cart);

/**
 * @brief Set k-point grid vectors.
 *
 * @param[in] h      Handler.
 * @param[in] nk1    Number of k-points along direction 1.
 * @param[in] nk2    Number of k-points along direction 2.
 * @param[in] nk3    Number of k-points along direction 3.
 * @param[in] kvecs  K-point vectors (nk1*nk2*nk3, Cartesian coordinates, in Bohr^-1).
 */
void librpa_set_kgrids_kvec(LibrpaHandler* h, int nk1, int nk2, int nk3, const double* kvecs);

/**
 * @brief Set the mapping from full k-point list to the irreducible sector.
 *
 * @param h         Handler.
 * @param nkpts     Number of k-points in the full Brillouin zone.
 * @param map_ibzk  Mapping to the k-point in the irreducible sector (0-based).
 *
 * Example: four-k-point case where the first two and last points are in the irreducible sector,
 *         and the third point is mapped to the second, then map_ibzk should be (0, 1, 1, 3).
 */
void librpa_set_ibz_mapping(LibrpaHandler* h, int nkpts, const int* map_ibzk);

/**
 * @brief Set local RI coefficients.
 *
 * @param h          Handler.
 * @param routing    Parallel routing strategy.
 * @param I          Atom I index (0-based).
 * @param J          Atom J index (0-based).
 * @param nbasis_i   Number of basis functions on atom I.
 * @param nbasis_j   Number of basis functions on atom J.
 * @param naux_mu    Number of auxiliary functions for mu index.
 * @param R          Lattice vector [R1, R2, R3].
 * @param Cs_in      RI triple coefficients.
 */
void librpa_set_lri_coeff(LibrpaHandler* h, LibrpaParallelRouting routing, int I, int J,
                          int nbasis_i, int nbasis_j, int naux_mu, const int R[3],
                          const double* Cs_in);

/**
 * @brief Set bare Coulomb matrix elements (atom-pair format).
 *
 * @param h               Handler.
 * @param ik              K-point index.
 * @param I               Atom I index.
 * @param J               Atom J index.
 * @param naux_mu         Number of auxiliary functions for mu.
 * @param naux_nu         Number of auxiliary functions for nu.
 * @param Vq_real_in      Real part of Coulomb matrix.
 * @param Vq_imag_in      Imaginary part of Coulomb matrix.
 * @param vq_threshold    Threshold for screening.
 */
void librpa_set_aux_bare_coulomb_k_atom_pair(LibrpaHandler* h, int ik, int I, int J, int naux_mu,
                                             int naux_nu, const double* Vq_real_in,
                                             const double* Vq_imag_in, double vq_threshold);

/**
 * @brief Set truncated Coulomb matrix elements (atom-pair format).
 *
 * @param h               Handler.
 * @param ik              K-point index.
 * @param I               Atom I index.
 * @param J               Atom J index.
 * @param naux_mu         Number of auxiliary functions for mu.
 * @param naux_nu         Number of auxiliary functions for nu.
 * @param Vq_real_in      Real part of Coulomb matrix.
 * @param Vq_imag_in      Imaginary part of Coulomb matrix.
 * @param vq_threshold    Threshold for screening.
 */
void librpa_set_aux_cut_coulomb_k_atom_pair(LibrpaHandler* h, int ik, int I, int J, int naux_mu,
                                            int naux_nu, const double* Vq_real_in,
                                            const double* Vq_imag_in, double vq_threshold);

/**
 * @brief Set bare Coulomb matrix elements (2D block format).
 *
 * @param h            Handler.
 * @param ik           K-point index.
 * @param mu_begin     Starting index for mu.
 * @param mu_end       Ending index for mu.
 * @param nu_begin     Starting index for nu.
 * @param nu_end       Ending index for nu.
 * @param Vq_real_in   Real part of Coulomb matrix.
 * @param Vq_imag_in   Imaginary part of Coulomb matrix.
 */
void librpa_set_aux_bare_coulomb_k_2d_block(LibrpaHandler* h, int ik, int mu_begin, int mu_end,
                                            int nu_begin, int nu_end, const double* Vq_real_in,
                                            const double* Vq_imag_in);

/**
 * @brief Set truncated Coulomb matrix elements (2D block format).
 *
 * @param h            Handler.
 * @param ik           K-point index.
 * @param mu_begin     Starting index for mu.
 * @param mu_end       Ending index for mu.
 * @param nu_begin     Starting index for nu.
 * @param nu_end       Ending index for nu.
 * @param Vq_real_in   Real part of Coulomb matrix.
 * @param Vq_imag_in   Imaginary part of Coulomb matrix.
 */
void librpa_set_aux_cut_coulomb_k_2d_block(LibrpaHandler* h, int ik, int mu_begin, int mu_end,
                                           int nu_begin, int nu_end, const double* Vq_real_in,
                                           const double* Vq_imag_in);

/**
 * @brief Set dielectric function on imaginary frequency axis.
 *
 * @param h               Handler.
 * @param nfreq           Number of frequency points.
 * @param omegas_imag     Imaginary frequency values.
 * @param dielectric_func Dielectric function values.
 */
void librpa_set_dielect_func_imagfreq(LibrpaHandler* h, int nfreq, const double* omegas_imag,
                                      const double* dielect_func);

/**
 * @brief Set k-points for band structure calculations.
 *
 * @param h              Handler.
 * @param n_kpts_band    Number of band k-points.
 * @param kfrac_band     Band k-point coordinates (fractional).
 */
void librpa_set_band_kvec(LibrpaHandler* h, int n_kpts_band, const double* kfrac_band);

/**
 * @brief Set occupation numbers and eigenvalues for band k-points.
 *
 * @param h            Handler.
 * @param n_spins      Number of spin channels.
 * @param n_kpts_band  Number of band k-points.
 * @param n_states     Number of states.
 * @param occ          Occupation numbers.
 * @param eig          Eigenvalues.
 */
void librpa_set_band_occ_eigval(LibrpaHandler* h, int n_spins, int n_kpts_band, int n_states,
                                const double* occ, const double* eig);

/**
 * @brief Set wavefunction for band k-point (separate real/imag).
 *
 * @param h              Handler.
 * @param ispin          Spin index.
 * @param ik_band        Band k-point index.
 * @param nstates_local  Local number of states.
 * @param nbasis_local   Local number of basis functions.
 * @param wfc_real       Real part of wavefunction.
 * @param wfc_imag       Imaginary part of wavefunction.
 */
void librpa_set_wfc_band(LibrpaHandler* h, int ispin, int ik_band, int nstates_local,
                         int nbasis_local, const double* wfc_real, const double* wfc_imag);

/**
 * @brief Set wavefunction for band k-point (packed complex).
 *
 * @param h              Handler.
 * @param ispin          Spin index.
 * @param ik_band        Band k-point index.
 * @param nstates_local  Local number of states.
 * @param nbasis_local   Local number of basis functions.
 * @param wfc_ri         Complex wavefunction (packed).
 */
void librpa_set_wfc_band_packed(LibrpaHandler* h, int ispin, int ik_band, int nstates_local,
                                int nbasis_local, const double* wfc_ri);

/**
 * @brief Reset band structure input data.
 *
 * @param h Handler.
 */
void librpa_reset_band_data(LibrpaHandler* h);

#ifdef __cplusplus
}
#endif
