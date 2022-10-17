/*
 * @file params.h
 * @brief parameters for controlling LibRPA calculation
 */
#pragma once
#include <string>

//! a simple struct to collect all runtime parameters
struct Params
{
    //! the task to perform in LibRPA.
    std::string task = "rpa";

    //! the number of frequency grid points
    int nfreq = 0;

    //! the type of time-frequency grids
    std::string tfgrids_type = "minimax";

    //! threshold of R-space Green's function when construcing.
    double gf_R_threshold = 1e-4;

    //! threshold of RI coefficient when parsing. The atomic block with maximal element smaller than it will be filtered.
    double cs_threshold = 1e-6;

    //! threshold of Coulomb matrix when parsing. The atom-pair block of Coulomb matrix with maximal element smaller than it will be filtered
    double vq_threshold = 0;

    //! threshold to filter when computing the square root of Coulomb matrix
    double sqrt_coulomb_threshold = 0;

    //! switch of using LibRI for chi0 calculation
    bool use_libri_chi0 = false;

    //! CS-matrix threshold parsed to RPA object of LibRI. 
    double libri_chi0_csm_threshold = 0.0;

    //! Cs threshold parsed to RPA object of LibRI. 
    double libri_chi0_threshold_C = 0.0;

    //! Green's function threshold parsed to RPA object of LibRI. 
    double libri_chi0_threshold_G = 0.0;

    //! switch of using ScaLAPACK for EcRPA calculation
    bool use_scalapack_ecrpa = false;

    void check_consistency();
    void print();
};

extern Params params;
