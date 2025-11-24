#pragma once
#include "librpa_enums.h"

// C APIs
#ifdef __cplusplus
extern "C" {
#endif

#define LIBRPA_MAX_STRLEN 200

typedef struct
{
    /* ============================================================================= */
    /* Common runtime control */

    //! Scheme of parallelization.
    LibrpaParallelRouting parallel_routing;

    //! Verbose level for output
    LibrpaVerbose output_level;

    //! Type of time/frequency grids.
    LibrpaTimeFreqGrid tfgrids_type;

    //! Number of frequency points
    int nfreq;

    //! Threshold for real-space LRI triple coefficients.
    double cs_threshold;

    //! Threshold for real-space Coulomb matrices
    double vq_threshold;

    //! Flag of parsing input and performing calculation in spin-orbit coupling formalism
    LibrpaSwitch use_soc;

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
    double libri_gw_threshold_C;

    //! Threshold of real-space Green's function to compute GW correlation self-energy using LibRI
    double libri_gw_threshold_G;

    //! Threshold of correlation screened Coulomb matrix to compute GW correlation self-energy using LibRI
    double libri_gw_threshold_Wc;

    // /* Output controls */
    // //! Output correlation self-energy matrix (reciprocal space, imaginary frequency domain)
    // LibrpaSwitch output_gw_sigc_mat;

    // //! Output correlation self-energy matrix in NAO (real space, imaginary time domain)
    // LibrpaSwitch output_gw_sigc_mat_rt;

    // //! Output correlation self-energy matrix in NAO (real space, imaginary frequency domain)
    // LibrpaSwitch output_gw_sigc_mat_rf;

} LibrpaOptions;

//! Initialize options to default values
void librpa_init_options(LibrpaOptions *opts);

#ifdef __cplusplus
}
#endif

// C++ APIs
#ifdef __cplusplus

namespace librpa
{

// Straighforward inheritance
// IMPORTANT: DO NOT add extra member variables here, which will break the inheritated data layout
// New control options should be put under the LibrpaOptions C structure
class Options : public ::LibrpaOptions
{
public:
    Options() { ::librpa_init_options(this); }
};

}

#endif
