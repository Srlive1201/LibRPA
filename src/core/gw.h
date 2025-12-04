/*
 * @file gw.h
 * @brief facilities to calculate self-energy operator.
 */
#pragma once
#include "../math/matrix_m.h"
#include "../mpi/base_blacs.h"
#include "atom.h"
#include "atomic_basis.h"
#include "meanfield.h"
#include "geometry.h"
#include "pbc.h"
#include "ri.h"
#include "timefreq.h"

namespace librpa_int
{

class G0W0
{
private:
    bool is_rspace_built_;
    bool is_kspace_built_;

public:
    const MeanField &mf;
    const AtomicBasis& atbasis_wfc;
    const PeriodicBoundaryData &pbc;
    const TFGrids &tfg;
    const MpiCommHandler &comm_h;

    std::string output_dir;

    double libri_threshold_C;
    double libri_threshold_Wc;
    double libri_threshold_G;

    bool output_sigc_mat;
    bool output_sigc_mat_rt;
    bool output_sigc_mat_rf;

    //! frequency-domain reciprocal-space correlation self-energy, indices [ispin][freq][k][I][J](n_I, n_J)
    // std::map<int, std::map<double, std::map<Vector3_Order<double>, atom_mapping<Matz>::pair_t_old>>> sigc_is_f_k_IJ;

    //! frequency-domain reciprocal-space correlation self-energy, indices [ispin][freq][R][I][J](n_I, n_J)
    std::map<int, std::map<double, std::map<Vector3_Order<int>, atom_mapping<Matz>::pair_t_old>>> sigc_is_f_R_IJ;

    //! correlation self-energy matrix in the basis of KS states, indices [ispin][ik][freq](n_bands, n_bands)
    std::map<int, std::map<int, std::map<double, Matz>>> sigc_is_ik_f_KS;

    void build_sigc_matrix_KS(const std::map<int, std::map<int, ComplexMatrix>> &wfc_target,
                              const std::vector<Vector3_Order<double>> &kfrac_target,
                              const Atoms &geometry,
                              const BlacsCtxtHandler &blacs_ctxt_h);
public:
    // Constructors
    G0W0(const MeanField &mf_in,
         const AtomicBasis& atbasis_wfc_in,
         const PeriodicBoundaryData &pbc_in,
         const TFGrids &tfg_in, const MpiCommHandler &comm_h_in);
    // delete copy/move constructors
    G0W0(const G0W0 &s_g0w0) = delete;
    G0W0(G0W0 &&s_g0w0) = delete;

    // delete assignment copy/move
    G0W0 operator=(const G0W0 &s_g0w0) const = delete;
    G0W0 operator=(G0W0 &&s_g0w0) = delete;

    //! Reset the real-space matrices
    void reset_rspace();

    //! Reset the k-space matrices
    void reset_kspace();

    //! Build the real-space correlation self-energy matrix on imaginary frequencies with space-time method using LibRI
    void build_spacetime(
        const LibrpaParallelRouting parallel_routing,
        const AtomicBasis &atbasis_abf,
        const Cs_LRI &LRI_Cs,
        std::map<double, std::map<Vector3_Order<double>, Matz>> &Wc_freq_q,
        const ArrayDesc &ad_Wc);

    //! build the correlation self-energy matrix in Kohn-Sham basis at the SCF k-points
    void build_sigc_matrix_KS_kgrid(const Atoms &geometry, const BlacsCtxtHandler &blacs_ctxt_h);

    //! build the correlation self-energy matrix in Kohn-Sham basis at the SCF k-points
    void build_sigc_matrix_KS_band(const std::map<int, std::map<int, ComplexMatrix>> &wfc_band,
                                   const std::vector<Vector3_Order<double>> &kfrac_band,
                                   const Atoms &geometry,
                                   const BlacsCtxtHandler &blacs_ctxt_h);
};

} /* end of namespace librpa_int */
