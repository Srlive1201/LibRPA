#pragma once

#include <functional>
#include <vector>

#include "dataset.h"
#include "librpa_options.h"

namespace librpa_int
{
//! Optional task composition at the resident full-response boundary.
using ResponseObserver = std::function<void(const Chi0QMap &)>;

//! Build the live exchange, full response, screening and correlation self-energy.
void build_g0w0_sigma(Dataset &dataset, const LibrpaOptions &opts,
                      const ResponseObserver &on_response = {});

//! Continue and solve the diagonal QP equation in the supplied numerical state.
//! Sigma is [spin][collected k][frequency][selected state], energies in Hartree.
void solve_g0w0_qpe_diagonal(const MeanField &mf, const TFGrids &tfg, const LibrpaOptions &opts,
                             const std::vector<cplxdb> &sigc_diag,
                             const std::vector<int> &iks_collect, int n_kpts_this,
                             const int *iks_this, int i_state_low, int n_states_calc,
                             const double *vxc, const double *vexx, double *sigc_re,
                             double *sigc_im, double *eqp, const char *context,
                             int *solver_info = nullptr);

void get_g0w0_qpe_kgrid(Dataset &dataset, const LibrpaOptions &opts, const int n_spins,
                        const int n_kpts_this, const int *iks_this, int i_state_low,
                        int i_state_high, const double *vxc, const double *vexx, double *sigc_re,
                        double *sigc_im, double *eqp);

void get_g0w0_spectral_function_kgrid(Dataset &dataset, const LibrpaOptions &opts,
                                      const int n_spins, const int n_kpts_this, const int *iks_this,
                                      int i_state_low, int i_state_high, const int n_omegas,
                                      const double *omegas, const double *vxc, const double *vexx,
                                      double *spectral_function, double *sigc);

void get_g0w0_qpe_band_k(Dataset &dataset, const LibrpaOptions &opts, const int n_spins,
                         const int n_kpts_band_this, const int *iks_band_this, int i_state_low,
                         int i_state_high, const double *vxc_band, const double *vexx_band,
                         double *sigc_band_re, double *sigc_band_im, double *eqp_band);

void get_g0w0_spectral_function_band_k(Dataset &dataset, const LibrpaOptions &opts,
                                       const int n_spins, const int n_kpts_band_this,
                                       const int *iks_band_this, int i_state_low, int i_state_high,
                                       const int n_omegas, const double *omegas,
                                       const double *vxc_band, const double *vexx_band,
                                       double *spectral_function_band, double *sigc_band);

}  // namespace librpa_int
