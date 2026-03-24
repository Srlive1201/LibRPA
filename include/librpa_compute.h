#pragma once
/**
 * @file librpa_compute.h
 * @brief Computing data APIs for LibRPA.
 */

#include "librpa_handler.h"
#include "librpa_options.h"

// C APIs
#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Construct and return frequency grids for numerical integration.
 *
 * Generates the frequency points and quadrature weights based on the grid type
 * specified in the options (e.g., Gauss-Legendre, Minimax, etc.).
 *
 * @param[in]  h        Handler.
 * @param[in]  p_opts   Runtime options.
 * @param[out] omegas  Array of frequency points (size: opts->nfreq).
 * @param[out] weights Array of quadrature weights (size: opts->nfreq).
 */
void librpa_get_imaginary_frequency_grids(LibrpaHandler *h, const LibrpaOptions *p_opts,
                                          double *omegas, double *weights);

/**
 * @brief Compute RPA correlation energy.
 *
 * Calculates the RPA correlation energy using the input data (wavefunctions,
 * Coulomb matrices, etc.) that has been set via the input parsing functions.
 *
 * @param[in]  h                           Handler.
 * @param[in]  p_opts                      Runtime options.
 * @param[in]  n_ibz_kpoints                Number of irreducible k-points.
 * @param[out] rpa_corr_ibzk_contrib_re   Real part of correlation energy per IBZ k-point.
 * @param[out] rpa_corr_ibzk_contrib_im   Imaginary part of correlation energy per IBZ k-point.
 * @return Total RPA correlation energy (real part).
 */
double librpa_get_rpa_correlation_energy(LibrpaHandler *h, const LibrpaOptions *p_opts,
                                         int n_ibz_kpoints, double *rpa_corr_ibzk_contrib_re,
                                         double *rpa_corr_ibzk_contrib_im);

//! Build exact-exchange matrix
/**
 * @param[in]  h                Pointer to LibRPA handler.
 * @param[in]  opts             Pointer to runtime options.
 */
void librpa_build_exx(LibrpaHandler *h, const LibrpaOptions *p_opts);

//! Obtain exact-exchange potential for selected states.
/**
 * @param[in]  h                Pointer to LibRPA handler.
 * @param[in]  opts             Pointer to runtime options.
 * @param[in]  n_spins          Number of spin channels.
 * @param[in]  n_kpts_this      Number of k-points to compute on this process.
 * @param[in]  iks_this         (Global) index of k-points that this process compute.
 *                              Each process can have different indices.
 *                              Must be a subset of k-points at which the eigenvetors are parsed.
 * @param[in]  i_state_low      Index of the first state to compute the potential (inclusive)
 * @param[in]  i_state_high     Index of the last state to compute the potential (exclusive)
 * @param[out] vexx             Exact-exchange potential for selected states.
 *                              It should be at least as long as n_spins * n_kpts_local * (i_state_high - i_state_low).
 */
void librpa_get_exx_pot_kgrid(LibrpaHandler *h, const LibrpaOptions *p_opts,
                              const int n_spins, const int n_kpts_this, const int *iks_this,
                              int i_state_low, int i_state_high, double *vexx);

//! Obtain exact-exchange potential for selected states at band k-points.
/**
 * @param[in]  h                Pointer to LibRPA handler.
 * @param[in]  opts             Pointer to runtime options.
 * @param[in]  n_spins          Number of spin channels.
 * @param[in]  n_kpts_band_this Number of k-points to compute on this process.
 * @param[in]  iks_band_this    (Global) index of k-points that this process compute.
 *                              Each process can have different indices.
 *                              Must be a subset of band k-points at which the eigenvetors are parsed.
 * @param[in]  i_state_low      Index of the first state to compute the potential (inclusive)
 * @param[in]  i_state_high     Index of the last state to compute the potential (exclusive)
 * @param[out] vexx_band        Exact-exchange potential for selected states at band k-points.
 *                              It should be at least as long as n_spins * n_kpts_band_this * (i_state_high - i_state_low).
 */
void librpa_get_exx_pot_band_k(LibrpaHandler *h, const LibrpaOptions *p_opts,
                               const int n_spins, const int n_kpts_band_this, const int *iks_band_this,
                               int i_state_low, int i_state_high, double *vexx_band);

//! Build self-energy matrix of G0W0, including the correlation and exchange contributions.
/**
 * @param[in]  h                Pointer to LibRPA handler.
 * @param[in]  opts             Pointer to runtime options.
 */
void librpa_build_g0w0_sigma(LibrpaHandler *h, const LibrpaOptions *p_opts);

//! Obtain correlation self-energies for selected states.
/**
 * @param[in]  h                Pointer to LibRPA handler.
 * @param[in]  opts             Pointer to runtime options.
 * @param[in]  n_spins          Number of spin channels.
 * @param[in]  n_kpts_this      Number of k-points to compute on this process.
 * @param[in]  iks_this         (Global) index of k-points that this process compute.
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
 * @param[out] sigc_im          Same as sigc_re, but for the imaginary part.
 */
void librpa_get_g0w0_sigc_kgrid(LibrpaHandler *h, const LibrpaOptions *p_opts,
                                const int n_spins, const int n_kpts_this,
                                const int *iks_this, int i_state_low, int i_state_high,
                                const double *vxc, const double *vexx, double *sigc_re, double *sigc_im);

//! Obtain correlation self-energies for selected states at band k-points.
/**
 * @param[in]  h                Pointer to LibRPA handler.
 * @param[in]  opts             Pointer to runtime options.
 * @param[in]  n_spins          Number of spin channels.
 * @param[in]  n_kpts_band_this Number of k-points to compute on this process.
 * @param[in]  iks_band_this    (Global) index of k-points that this process compute.
 *                              Each process can have different indices.
 *                              Must be a subset of k-points at which the eigenvetors are parsed.
 * @param[in]  i_state_low      Index of the first state to compute the potential (inclusive)
 * @param[in]  i_state_high     Index of the last state to compute the potential (exclusive)
 * @param[in]  vxc_band         exchange-correlation potential of the selected states at band k-points.
 * @param[in]  vexx_band        Exact-exchange potential for the selected states at band k-points.
 *                              It should be at least as long as n_spins * n_kpts_band_this * (i_state_high - i_state_low).
 *                              It can be obtained using librpa_get_exx_pot_kgrid.
 * @param[out] sigc_band_re     Real-part of the correlation self-energy for the selected states.
 *                              It should be at least as long as n_spins * n_kpts_band_this * (i_state_high - i_state_low).
 * @param[out] sigc_band_im     Same as sigc_band_re, but for the imaginary part.
 */
void librpa_get_g0w0_sigc_band_k(LibrpaHandler *h, const LibrpaOptions *p_opts,
                                 const int n_spins, const int n_kpts_band_this,
                                 const int *iks_band_this, int i_state_low, int i_state_high,
                                 const double *vxc_band, const double *vexx_band, double *sigc_band_re, double *sigc_band_im);

#ifdef __cplusplus
}
#endif
