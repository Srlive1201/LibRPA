/*!
 * @file params.h
 * @brief parameters for controlling LibRPA calculation
 */
#ifndef PARAMS_H
#define PARAMS_H
#include <iostream>
#include <cstdarg>
#include <fstream>

#pragma once
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

    static bool binary_input;

    //! CS-matrix threshold parsed to RPA object of LibRI.
    static double libri_chi0_threshold_CSM;

    //! Cs threshold parsed to RPA object of LibRI.
    static double libri_chi0_threshold_C;

    //! Green's function threshold parsed to RPA object of LibRI.
    static double libri_chi0_threshold_G;

    //! switch of using ScaLAPACK for EcRPA calculation
    static bool use_scalapack_ecrpa;

    //! CS-matrix threshold parsed to EXX object of LibRI.
    static double libri_exx_threshold_CSM;

    //! Cs threshold parsed to EXX object of LibRI.
    static double libri_exx_threshold_C;

    //! Density matrix threshold parsed to EXX object of LibRI.
    static double libri_exx_threshold_D;

    //! Coulomb matrix threshold parsed to EXX object of LibRI.
    static double libri_exx_threshold_V;

    //! switch of using ScaLAPACK for computing Wc from chi0
    static bool use_scalapack_gw_wc;

    //! switch of run-time debug mode
    static bool debug;

    //! switch of output correlation self-energy matrix
    static bool output_gw_sigc_mat;

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

    static void check_consistency();
    static void print();
};


// NOTE:(MYZ) Can we move customPrint to other file?

static void customPrint(const char* format, ...) {
    va_list args;
    va_start(args, format);
    const bool output_stdout = Params::output_file == "stdout";

    char buffer[1024];
    vsnprintf(buffer, sizeof(buffer), format, args);

    static std::ofstream outputFile;
    std::ostream &os = output_stdout ? std::cout : outputFile;
    if (!output_stdout) outputFile.open(Params::output_file, std::ios_base::app);
    os << buffer;
    os.flush();
    va_end(args);
}

// #define printf customPrint

#endif
