#pragma once
#include "librpa_handler.h"
#include "librpa_options.h"
#include "librpa_func.h"

// C APIs
#ifdef __cplusplus
extern "C" {
#endif

LIBRPA_C_H_FUNC_WRAP_WOPT(double, librpa_get_rpa_correlation_energy,
                          int n_ibz_kpoints, double *rpa_corr_ibzk_contrib_re,  double *rpa_corr_ibzk_contrib_im);

//! Build exact-exchange matrix
LIBRPA_C_H_FUNC_WRAP_WOPT_NOPAR(void, librpa_build_exx);

// from i_state_low (included) to i_state_high (excluded)
// vexx should be at least as long as n_spins * n_kpoints_local * (i_state_high - i_state_low).
LIBRPA_C_H_FUNC_WRAP_WOPT(void, librpa_get_exx_pot_kgrid,
                          const int n_spins, const int n_kpoints_local,
                          const int *iks_local, int i_state_low, int i_state_high,
                          double *vexx);

//! Build self-energy matrix of G0W0
LIBRPA_C_H_FUNC_WRAP_WOPT_NOPAR(void, librpa_build_g0w0_sigma);

LIBRPA_C_H_FUNC_WRAP_WOPT(void, librpa_get_g0w0_qpe_kgrid,
                          const int n_kpoints_local, const int *iks_local,
                          int i_state_low, int i_state_high, const double *vxc,
                          double *qpe_re, double *qpe_im);

#ifdef __cplusplus
}
#endif
