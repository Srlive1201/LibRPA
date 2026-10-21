/*
 * @file gw.h
 * @brief facilities to calculate self-energy operator.
 */
#pragma once
#include "../math/matrix_m.h"
#include "../mpi/base_blacs.h"
#include "../mpi/kpoint_blacs_parallel_context.h"
#include "atom.h"
#include "atomic_basis.h"
#include "symmetry_context.h"
#include "meanfield.h"
#include "geometry.h"
#include "pbc.h"
#include "qpoint_view.h"
#include "ri.h"
#include "timefreq.h"

namespace librpa_int
{

class G0W0
{
private:
    using SigcRspaceMap =
        std::map<int, std::map<int, std::map<int, std::map<double, ap_p_map<std::map<Vector3_Order<int>, Matz>>>>>>;

    bool is_eigvec_k_distributed_;
    bool is_rspace_built_;
    bool is_kspace_built_;
    int output_sigc_ks_kf_band_index_;
    std::string sigc_kspace_source_;

    //! frequency-domain reciprocal-space correlation self-energy, indices [ispin][ispinor_bra][ispinor_ket][freq][R][I][J](n_I, n_J)
    // Sparse storage from LibRI calculation
    SigcRspaceMap sigc_is_f_IJ_R;

    void collect_sigc_rf_output_shards();
    void write_sigc_rf_output_files() const;

    void build_sigc_matrix_KS(const std::map<int, std::map<int, std::map<int, ComplexMatrix>>> &wfc_target,
                              const std::vector<Vector3_Order<double>> &kfrac_target,
                              const AtomPairBvKRemap<atom_t> &bvk_remap);

    void build_sigc_matrix_KS_blacs(const std::map<int, std::map<int, std::map<int, ComplexMatrix>>> &wfc_target,
                                    const std::vector<Vector3_Order<double>> &kfrac_target,
                                    const AtomPairBvKRemap<atom_t> &bvk_remap,
                                    const BlacsCtxtHandler &blacs_ctxt_h,
                                    bool use_gpu_replace_scalapack,
                                    const std::string &source);

public:
    const MeanField &mf;
    //! Array descriptor of wave functions saved in the MeanField object
    const ArrayDesc &desc_wfc;
    //! Array descriptor of band-path wave functions.
    const ArrayDesc &desc_band_wfc;
    const AtomicBasis& atbasis_wfc;
    const PeriodicBoundaryData &pbc;
    const SymmetryContext &symmetry_context;
    const bool use_symmetry_context;
    //! Representation of the supplied screening, independent of Sigma reduction.
    const SymmetryQPointView qpoint_view;
    const TFGrids &tfg;
    const MpiCommHandler &comm_h;
    const KPointBlacsParallelContext &kblacs_ctxt;
    const KPointBlacsParallelContext &band_kblacs_ctxt;

    std::string output_dir;

    double libri_threshold_C;
    double libri_threshold_Wc;
    double libri_threshold_G;
    int nbands_G = -1;
    int n_bands_exclude = 0;

    bool output_sigc_ks_mat_kf;
    bool output_sigc_ks_kf;
    int istate_output_mat_start;
    int istate_output_mat_end;
    bool output_sigc_mat_kf;
    bool output_sigc_mat_rt;
    bool output_sigc_mat_rf;
    bool output_wc_rf;
    bool output_wc_rf_atom_pair;
    int ifreq_output_wc_start;
    int ifreq_output_wc_end;

    //! frequency-domain reciprocal-space correlation self-energy, indices [ispin][freq][k][I][J](n_I, n_J)
    // std::map<int, std::map<double, std::map<Vector3_Order<double>, atom_mapping<Matz>::pair_t_old>>> sigc_is_f_k_IJ;

    //! descriptor for the distributed KS self-energy matrix blocks
    ArrayDesc desc_sigc_is_ik_f_KS;

    //! correlation self-energy matrix in the basis of KS states, indices [ispin][ik][freq](local n_bands, local n_bands)
    std::map<int, std::map<int, std::map<double, Matz>>> sigc_is_ik_f_KS;

    //! correlation self-energy diagonal in the basis of KS states, indices [ispin][ik][freq](n_bands)
    std::map<int, std::map<int, std::map<double, std::vector<cplxdb>>>> sigc_diag_is_ik_f_KS;

public:
    // Constructors
    G0W0(const MeanField &mf_in, const AtomicBasis &atbasis_wfc_in,
         const PeriodicBoundaryData &pbc_in, const SymmetryContext &symmetry_context_in,
         const TFGrids &tfg_in, const KPointBlacsParallelContext &kblacs_ctxt_in,
         const KPointBlacsParallelContext &band_kblacs_ctxt_in, const ArrayDesc &desc_wfc_in,
         const ArrayDesc &desc_band_wfc_in, bool is_eigvec_k_distributed,
         bool use_symmetry_context_in = true,
         const SymmetryQPointView *screening_qpoint_view = nullptr);
    // // delete copy/move constructors
    // G0W0(const G0W0 &s_g0w0) = delete;
    // G0W0(G0W0 &&s_g0w0) = delete;

    // // delete assignment copy/move
    // G0W0 operator=(const G0W0 &s_g0w0) const = delete;
    // G0W0 operator=(G0W0 &&s_g0w0) = delete;

    //! Reset the real-space matrices
    void reset_rspace();

    //! Reset the k-space matrices
    void reset_kspace();

    //! Check if the real-space self-energy matrix is built
    bool is_rspace_built() const { return is_rspace_built_; }

    //! Source tag for the currently cached K-S self-energy matrices.
    const std::string &sigc_kspace_source() const { return sigc_kspace_source_; }

    //! Read real-space imaginary-frequency correlation self-energy matrices from disk
    void read_sigc(const std::string &input_dir);

    //! Build the real-space correlation self-energy matrix on imaginary frequencies with space-time method using LibRI
    //!
    //! Wc_freq_q is destroyed on exit.
    void build_spacetime(
        const LibrpaParallelRouting parallel_routing,
        const AtomicBasis &atbasis_abf,
        const Cs_LRI &LRI_Cs,
        std::map<double, std::map<Vector3_Order<double>, Matz>> &Wc_freq_q,
        const ArrayDesc &ad_Wc,
        std::map<double,
                 atom_mapping<std::map<Vector3_Order<double>, Matz>>::pair_t_old>
            *Wc_freq_q_atom_pair = nullptr,
        std::map<Vector3_Order<double>, ComplexMatrix> *sinvS = nullptr,
        const AtomicBasis *basis_aux_compressed = nullptr,
        const AtomicBasis *basis_aux_unfold = nullptr,
        const BlacsCtxtHandler *blacs_ctxt_h = nullptr,
        const ArrayDesc *desc_wfc_in = nullptr);

    void build_sigc_matrix_KS_kgrid(const Atoms &geometry = Atoms());
    void build_sigc_matrix_KS_band(const std::map<int, std::map<int, std::map<int, ComplexMatrix>>> &wfc_band,
                                   const std::vector<Vector3_Order<double>> &kfrac_band,
                                   const AtomPairBvKRemap<atom_t> &bvk_remap,
                                   const std::vector<int> *output_iks = nullptr);
    void build_sigc_matrix_KS_kgrid_blacs(const BlacsCtxtHandler &blacs_ctxt_h,
                                          bool use_gpu_replace_scalapack = false);
    void build_sigc_matrix_KS_band_blacs(const std::map<int, std::map<int, std::map<int, ComplexMatrix>>> &wfc_band,
                                         const std::vector<Vector3_Order<double>> &kfrac_band,
                                         const AtomPairBvKRemap<atom_t> &bvk_remap,
                                         const BlacsCtxtHandler &blacs_ctxt_h,
                                         bool use_gpu_replace_scalapack = false,
                                         const std::vector<int> *output_iks = nullptr);
    void write_sigc_matrices_KS_binary(const std::string &output_dir,
                                       const std::string &source) const;
};

} /* end of namespace librpa_int */
