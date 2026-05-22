#pragma once
/**
 * @file librpa_options.h
 * @brief Runtime options for LibRPA calculations.
 *
 * This file defines the LibrpaOptions structure that controls runtime behavior
 * for RPA, exact exchange, and G0W0 calculations.
 */

#include "librpa_enums.h"

// C APIs
#ifdef __cplusplus
extern "C" {
#endif

/** Maximum length for string parameters (e.g., output directory path). */
#define LIBRPA_MAX_STRLEN 200

// NOTE: in case the data layout of LibrpaOptions is changed,
// its Fortran binding should be adapted accordingly.

/**
 * @brief Runtime options structure for LibRPA calculations.
 *
 * This structure contains all configurable parameters for controlling
 * RPA, exact exchange (EXX), and G0W0 calculations. Initialize with
 * librpa_init_options() and modify as needed before passing to
 * computation functions.
 *
 * @note The data layout must match the Fortran binding.
 */
typedef struct
{
    /* ============================================================================= */
    /* Common runtime control */

    //! Output directory
    char output_dir[LIBRPA_MAX_STRLEN];

    //! Scheme of parallelization.
    LibrpaParallelRouting parallel_routing;

    //! Verbose level for output
    LibrpaVerbose output_level;

    //! Threshold for real-space Coulomb matrices
    double vq_threshold;

    //! Flag of using spinor wavefunction to construct objects, like Green's function and self-energy.
    LibrpaSwitch use_spinor_wfc;

    //! Flag to specify parallel distribution of eigenvectors of SCF starting point
    LibrpaSwitch use_kpara_scf_eigvec;

    //! Type of time/frequency grids.
    LibrpaTimeFreqGrid tfgrids_type;

    //! Number of time/frequency points
    int nfreq;

    /* ============================================================================= */
    /* Parameters for time and freqeuncy grids.
     * Different grid types will use none/part/all of the parameters. */
    double tfgrids_freq_min;
    double tfgrids_freq_interval;
    double tfgrids_freq_max;
    double tfgrids_time_min;
    double tfgrids_time_interval;

    //! Minimal transition energy when generating minimax grids.
    /*!
     * When set negative (the default), the minimal energy is decided by the mean-field input
     * and set to the energy gap.
     * For now, it is suggested to manually set this option for gapless systems.
     */
    double minimax_emin;

    //! Maximal transition energy when generating minimax grids.
    /*!
     * This option is introduced for test of supercell to fix the minimax frequency points
     */
    double minimax_emax;

    /* ============================================================================= */
    /* Coulomb matrices usage */

    //! Switch of using full Coulomb interaction in exact-exchange operator.
    LibrpaSwitch use_fullcoul_exx;

    //! Switch of using full Coulomb interaction in $\varepsilon = 1 - v \chi^0$
    LibrpaSwitch use_fullcoul_eps;

    //! Switch of using full Coulomb interaction in $W^c = (\varepsilon^{-1} - 1) v$.
    LibrpaSwitch use_fullcoul_wc;

    /* ============================================================================= */
    /* Sum of states */

    //! Maximal number of bands for computing response function
    int n_bands_chi0;

    //! Maximal number of bands for computing correlation self-energy
    int n_bands_sigc;

    /* ============================================================================= */
    /* RPA specific */

    //! Threshold for real-space Green's function in response function calculation
    double gf_threshold;

    //! Flag to control whether to use ScaLAPACK to compute RPA correlation energy
    LibrpaSwitch use_scalapack_ecrpa;

    /* ============================================================================= */
    /* ABF compression */

    LibrpaSwitch use_shrink_abfs;
    LibrpaSwitch use_shrink_chi;

    /* ============================================================================= */
    /* GW specific */

    //! Number of parameters for analytic continuation
    int n_params_anacon;

    //! Flag of using ScaLAPACK for computing Wc from chi0
    LibrpaSwitch use_scalapack_gw_wc;

    //! Flag of using cholesky factorization for computing Wc from chi0
    LibrpaSwitch use_cholesky_gw_wc;

    //! Flag of replacing head of screened interaction by macroscopic dielectric function
    LibrpaSwitch replace_w_head;

    //! Option of computing dielectric function on imaginary axis
    /*!
     * Available values:
     * - 0: use data directly read from input
     * - 1: dielectric model fitting from input data
     * - 2: cubic-spline interpolation from input data
     */
    int option_dielect_func;

    //! Switch of using 2D dielectric function
    LibrpaSwitch use_2d_dielectric;

    //! Flag of loading correlation self-energy matrix (real-space, imaginary frequency) from file
    LibrpaSwitch load_sigc_from_file;

    //! Threshold of eigenvalues to perform square root of Coulomb matrices
    double sqrt_coulomb_threshold;

    /* ============================================================================= */
    /* LibRI related */

    //! Threshold of real-space LRI triple coefficients to compute response function using LibRI
    double libri_chi0_threshold_C;

    //! Threshold of real-space Green's function to compute response function using LibRI
    double libri_chi0_threshold_G;

    //! Threshold of real-space LRI triple coefficients to compute exact exchange using LibRI
    double libri_exx_threshold_C;

    //! Threshold of real-space density matrices to compute exact exchange using LibRI
    double libri_exx_threshold_D;

    //! Threshold of real-space Coulomb matrices to compute exact exchange using LibRI
    double libri_exx_threshold_V;

    //! Threshold of real-space LRI triple coefficients to compute GW correlation self-energy using LibRI
    double libri_g0w0_threshold_C;

    //! Threshold of real-space Green's function to compute GW correlation self-energy using LibRI
    double libri_g0w0_threshold_G;

    //! Threshold of correlation screened Coulomb matrix to compute GW correlation self-energy using LibRI
    double libri_g0w0_threshold_Wc;

    /* Output controls */
    //! Output correlation self-energy matrix (reciprocal space, imaginary frequency domain)
    LibrpaSwitch output_gw_sigc_mat;

    //! Output correlation self-energy matrix in NAO (real space, imaginary time domain)
    LibrpaSwitch output_gw_sigc_mat_rt;

    //! Output correlation self-energy matrix in NAO (real space, imaginary frequency domain)
    LibrpaSwitch output_gw_sigc_mat_rf;

    //! Switch of outputting Wc matrix in Abs (real space, imaginary frequency domain)
    /*!
     * Available values:
     * - 0: do not output
     * - 1: output lowerest frequency
     * - 2: output all frequencies
     */
    int option_output_Wc_Rf_mat;

} LibrpaOptions;

/**
 * @brief Initialize runtime options to default values.
 *
 * Sets all options to their default settings. Must be called before
 * modifying options and passing to computation functions.
 *
 * @param[out] opts Pointer to the options structure to initialize.
 */
void librpa_init_options(LibrpaOptions *opts);

/**
 * @brief Set the output directory for LibRPA results.
 *
 * @param[in,out] opts      Pointer to the options structure.
 * @param[in]     output_dir Path to directory where output files will be written.
 */
void librpa_set_output_dir(LibrpaOptions *opts, const char *output_dir);

#ifdef __cplusplus
}
#endif
