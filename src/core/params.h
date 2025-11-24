/*!
 * @file params.h
 * @brief parameters for controlling LibRPA calculation
 */
#ifndef PARAMS_H
#define PARAMS_H

#include <string>

//! a simple struct to collect all runtime parameters
struct Params
{
    //! the task to perform in LibRPA.
    static std::string task;

    static std::string DFT_software;
    //! the path of file to store librpa mainly output
    static std::string output_file;

    //! the path of directory to store librpa output files
    static std::string output_dir;

    //! the number of frequency grid points
    static int nfreq;

    //! the type of time-frequency grids
    static std::string tfgrids_type;

    //! the number of parameters for analytic continuation
    static int n_params_anacon;

    //! type of parallel routing
    static std::string parallel_routing;

    //! threshold of R-space Green's function when construcing.
    static double gf_R_threshold;

    //! threshold of RI coefficient when parsing. The atomic block with maximal element smaller than it will be filtered.
    static double cs_threshold;

    //! threshold of Coulomb matrix when parsing. The atom-pair block of Coulomb matrix with maximal element smaller than it will be filtered
    static double vq_threshold;

    //! threshold to filter when computing the square root of Coulomb matrix
    static double sqrt_coulomb_threshold;

    //! Cs threshold parsed to RPA object of LibRI.
    static double libri_chi0_threshold_C;

    //! Green's function threshold parsed to RPA object of LibRI.
    static double libri_chi0_threshold_G;

    //! switch of using ScaLAPACK for EcRPA calculation
    static bool use_scalapack_ecrpa;

    //! Cs threshold parsed to EXX object of LibRI.
    static double libri_exx_threshold_C;

    //! Density matrix threshold parsed to EXX object of LibRI.
    static double libri_exx_threshold_D;

    //! Coulomb matrix threshold parsed to EXX object of LibRI.
    static double libri_exx_threshold_V;

    //! Cs threshold parsed to G0W0 object of LibRI.
    static double libri_g0w0_threshold_C;

    //! Green's function threshold parsed to G0W0 object of LibRI.
    static double libri_g0w0_threshold_G;

    //! Screened Coulomb matrix threshold parsed to EXX object of LibRI.
    static double libri_g0w0_threshold_Wc;

    //! switch of using ScaLAPACK for computing Wc from chi0
    static bool use_scalapack_gw_wc;

    //! switch of run-time debug mode
    static bool debug;

    //! switch of replacing head of screened interaction by macroscopic dielectric function
    static bool replace_w_head;

    //! option of computing dielectric function on imaginary axis
    /*!
     * Available values:
     * - 0: direct read from input
     * - 1: dielectric model fitting
     * - 2: cubic-spline interpolation
     */
    static int option_dielect_func;

    /* ==========================================================
     * output options
     */
    //! output correlation self-energy matrix (reciprocal space, imaginary frequency domain)
    static bool output_gw_sigc_mat;

    //! output correlation self-energy matrix in NAO (real space, imaginary time domain)
    static bool output_gw_sigc_mat_rt;

    //! output correlation self-energy matrix in NAO (real space, imaginary frequency domain)
    static bool output_gw_sigc_mat_rf;

    static void check_consistency();
    static void print();
};

#endif
