#pragma once
#include "librpa_handler.h"
#include "librpa_options.h"
#include "librpa_func.h"

// C APIs
#ifdef __cplusplus
extern "C" {
#endif

LIBRPA_C_H_FUNC_WRAP_WOPT(void, librpa_get_imaginary_frequency_grids,
                          double *omegas, double *weights);

LIBRPA_C_H_FUNC_WRAP_WOPT(double, librpa_get_rpa_correlation_energy,
                          int n_ibz_kpoints, double *rpa_corr_ibzk_contrib_re,  double *rpa_corr_ibzk_contrib_im);

//! Build exact-exchange matrix
/**
 * @param[in]  h                Pointer to LibRPA handler.
 * @param[in]  opts             Pointer to runtime options.
 */
LIBRPA_C_H_FUNC_WRAP_WOPT_NOPAR(void, librpa_build_exx);

//! Obtain exact-exchange potential for selected states.
/**
 * @param[in]  h                Pointer to LibRPA handler.
 * @param[in]  opts             Pointer to runtime options.
 * @param[in]  n_spins          Number of spin channels.
 * @param[in]  n_kpoints_local  Number of local k-points.
 * @param[in]  iks_local        Index of local k-points.
 *                              Each process can have different indices.
 *                              Must be a subset of k-points at which the eigenvetors are parsed.
 * @param[in]  i_state_low      Index of the first state to compute the potential (inclusive)
 * @param[in]  i_state_high     Index of the last state to compute the potential (exclusive)
 * @param[out] vexx             Exact-exchange potential for selected states.
 *                              It should be at least as long as n_spins * n_kpoints_local * (i_state_high - i_state_low).
 */
LIBRPA_C_H_FUNC_WRAP_WOPT(void, librpa_get_exx_pot_kgrid,
                          const int n_spins, const int n_kpoints_local,
                          const int *iks_local, int i_state_low, int i_state_high,
                          double *vexx);

//! Build self-energy matrix of G0W0, including the correlation and exchange contributions.
/**
 * @param[in]  h                Pointer to LibRPA handler.
 * @param[in]  opts             Pointer to runtime options.
 */
LIBRPA_C_H_FUNC_WRAP_WOPT_NOPAR(void, librpa_build_g0w0_sigma);

//! Obtain quasi-particle energies for selected states.
/**
 * @param[in]  h                Pointer to LibRPA handler.
 * @param[in]  opts             Pointer to runtime options.
 * @param[in]  n_spins          Number of spin channels.
 * @param[in]  n_kpoints_local  Number of local k-points.
 * @param[in]  iks_local        Index of local k-points.
 *                              Each process can have different indices.
 *                              Must be a subset of k-points at which the eigenvetors are parsed.
 * @param[in]  i_state_low      Index of the first state to compute the potential (inclusive)
 * @param[in]  i_state_high     Index of the last state to compute the potential (exclusive)
 * @param[in]  vxc              exchange-correlation potential of the selected states.
 * @param[in]  vexx             Exact-exchange potential for the selected states.
 *                              It should be at least as long as n_spins * n_kpoints_local * (i_state_high - i_state_low).
 *                              It can be obtained using librpa_get_exx_pot_kgrid.
 * @param[out] sigc_re          Real-part of the correlation self-energy for the selected states.
 *                              It should be at least as long as n_spins * n_kpoints_local * (i_state_high - i_state_low).
 * @param[out] sigc_im          Same as sigc_re, but for imaginary part.
 */
LIBRPA_C_H_FUNC_WRAP_WOPT(void, librpa_get_g0w0_qpe_kgrid,
                          const int n_spins, const int n_kpoints_local,
                          const int *iks_local, int i_state_low, int i_state_high,
                          const double *vxc, const double *vexx, double *sigc_re, double *sigc_im);

#ifdef __cplusplus
}
#endif
