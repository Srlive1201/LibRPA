#pragma once
#include <stddef.h>
#include "librpa_enums.h"
#include "librpa_handler.h"
#include "librpa_options.h"
#include "librpa_func.h"

// C APIs
#ifdef __cplusplus
extern "C" {
#endif

LIBRPA_C_H_FUNC_WRAP(void, librpa_set_scf_dimension,
                     int nspins, int nkpts, int nstates, int nbasis,
                     int st_istate, int nstates_local, int st_ibasis, int nbasis_local);

LIBRPA_C_H_FUNC_WRAP(void, librpa_set_wg_ekb_efermi,
                     int nspins, int nkpts, int nstates, const double* wg, const double* ekb, double efermi);

LIBRPA_C_H_FUNC_WRAP(void, librpa_set_wfc,
                     int ispin, int ik, int nstates_local, int nbasis_local,
                     const double* wfc_real, const double* wfc_imag);

LIBRPA_C_H_FUNC_WRAP(void, librpa_set_ao_basis_wfc, int natoms, const size_t *nbs_wfc);

LIBRPA_C_H_FUNC_WRAP(void, librpa_set_ao_basis_aux, int natoms, const size_t *nbs_aux);

LIBRPA_C_H_FUNC_WRAP(void, librpa_set_latvec_and_G,
                     const double lat_mat[9], const double G_mat[9]);

LIBRPA_C_H_FUNC_WRAP(void, librpa_set_atoms,
                     int natoms, const int *types, const double *posi_cart);

LIBRPA_C_H_FUNC_WRAP(void, librpa_set_kgrids_kvec,
                     int nk1, int nk2, int nk3, const double* kvecs);

/*! @brief Set the mapping from full k-point list to the irreducbile sector.
 *
 * @param nkpts     number of k-points in the full Brillouin zone
 * @param map_ibzk  mapping to the k-point in the irreducible sector

 * Example: four-k-point case where the first two and last points are in the irreducbile sector,
 *          and the third point is mapped to the second, then map_ibzk should be (0, 1, 1, 3)
 */
LIBRPA_C_H_FUNC_WRAP(void, librpa_set_ibz_mapping,
                     int nkpts, const int* map_ibzk);

LIBRPA_C_H_FUNC_WRAP(void, librpa_set_lri_coeff,
                     LibrpaParallelRouting routing,
                     int I, int J, int nbasis_i, int nbasis_j, int naux_mu, const int R[3], const double* Cs_in);

LIBRPA_C_H_FUNC_WRAP(void, librpa_set_aux_bare_coulomb_k_atom_pair,
                     int ik, int I, int J, int naux_mu, int naux_nu,
                     const double* Vq_real_in, const double* Vq_imag_in, double vq_threshold);

LIBRPA_C_H_FUNC_WRAP(void, librpa_set_aux_cut_coulomb_k_atom_pair,
                     int ik, int I, int J, int naux_mu, int naux_nu,
                     const double* Vq_real_in, const double* Vq_imag_in, double vq_threshold);

LIBRPA_C_H_FUNC_WRAP(void, librpa_set_aux_bare_coulomb_k_2d_block,
                     int ik, int mu_begin, int mu_end, int nu_begin, int nu_end,
                     const double* Vq_real_in, const double* Vq_imag_in);

LIBRPA_C_H_FUNC_WRAP(void, librpa_set_aux_cut_coulomb_k_2d_block,
                     int ik, int mu_begin, int mu_end, int nu_begin, int nu_end,
                     const double* Vq_real_in, const double* Vq_imag_in);

// LIBRPA_C_H_FUNC_WRAP(void, get_rpa_correlation_energy,
//                          double *rpa_corr, double *rpa_corr_irk_contrib);

#ifdef __cplusplus
}
#endif
