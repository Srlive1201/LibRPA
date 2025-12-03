#pragma once
#include "librpa_enums.h"

// C APIs
#ifdef __cplusplus
extern "C" {
#endif

#define LIBRPA_MAX_STRLEN 200

// NOTE: in case the data layout of LibrpaOptions is changed,
// its Fortran binding should be adapted acoordingly.

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

    //! Threshold for real-space LRI triple coefficients.
    double cs_threshold;

    //! Threshold for real-space Coulomb matrices
    double vq_threshold;

    //! Flag of parsing input and performing calculation in spin-orbit coupling formalism
    LibrpaSwitch use_soc;

    //! Type of time/frequency grids.
    LibrpaTimeFreqGrid tfgrids_type;

    //! Number of time/frequency points
    int nfreq;

    //! Parameters for time and freqeuncy grids. Different grid types will use none/part/all of the parameters.
    double tfgrids_freq_min;
    double tfgrids_freq_interval;
    double tfgrids_freq_max;
    double tfgrids_time_min;
    double tfgrids_time_interval;

    /* ============================================================================= */
    /* RPA specific */

    //! Threshold for real-space Green's function in response function calculation
    double gf_threshold;

    //! Flag to control whether to use ScaLAPACK to compute RPA correlation energy
    LibrpaSwitch use_scalapack_ecrpa;

    /* ============================================================================= */
    /* GW specific */

    //! Number of parameters for analytic continuation
    int n_params_anacon;

    //! Flag of using ScaLAPACK for computing Wc from chi0
    LibrpaSwitch use_scalapack_gw_wc;

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

} LibrpaOptions;

//! Initialize runtime options to default values
void librpa_init_options(LibrpaOptions *opts);

//! Helper function to set output directory in the runtime options
void librpa_set_output_dir(LibrpaOptions *opts, const char *output_dir);

#ifdef __cplusplus
}
#endif
