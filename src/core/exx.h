/*
 * @file exx.h
 * @brief utilities for computing exact exchange energies, including orbital and total energies.
 */
#include "../math/matrix_m.h"
#include "../mpi/base_blacs.h"
#include "../mpi/kpoint_blacs_parallel_context.h"
#include "atomic_basis.h"
#include "symmetry_context.h"
#include "meanfield.h"
#include "pbc.h"
#include "geometry.h"
#include "ri.h"

namespace librpa_int
{

class Exx
{
    private:
        bool is_mf_eigvec_k_distributed_;
        bool is_rspace_built_;
        bool is_kspace_built_;
        bool is_rspace_redist_for_KS_;
        bool is_rspace_redist_blacs_;

        void build_KS(const std::map<int, std::map<int, std::map<int, ComplexMatrix>>> &wfc_target,
                      const std::vector<Vector3_Order<double>> &kfrac_target,
                      const AtomPairBvKRemap<atom_t> &bvk_remap);
        void build_KS_blacs(
            const std::map<int, std::map<int, std::map<int, ComplexMatrix>>> &wfc_target,
            const std::vector<Vector3_Order<double>> &kfrac_target,
            const AtomPairBvKRemap<atom_t> &bvk_remap,
            const BlacsCtxtHandler &blacs_ctxt_h,
            bool use_gpu_replace_scalapack = false,
            const std::string &source = "kgrid");

    public:
        //! refenrence to the MeanField object to compute density matrix
        const MeanField &mf;
        //! Array descriptor of wave functions saved in the MeanField object
        const ArrayDesc &desc_wfc;
        //! Array descriptor of band-path wave functions.
        const ArrayDesc &desc_band_wfc;
        const AtomicBasis &atbasis_wfc;
        const PeriodicBoundaryData &pbc;
        const SymmetryContext &symmetry_context;
        const bool use_symmetry_context;
        const MpiCommHandler &comm_h;
        const KPointBlacsParallelContext &kblacs_ctxt;
        const KPointBlacsParallelContext &band_kblacs_ctxt;

        double libri_threshold_C;
        double libri_threshold_V;
        double libri_threshold_D;

        std::string output_dir;
        bool output_exx_mat_k = false;

        //! exact-exchange Hamiltonian in real space, dimension (nspins, nspinors, nspinors, I, J, R, nao_I, nao_J)
        std::map<int, std::map<int, std::map<int, std::map<atom_t, std::map<atom_t, std::map<Vector3_Order<int>, Matd>>>>>> exx_IJR;
        std::map<int, std::map<int, std::map<int, std::map<atom_t, std::map<atom_t, std::map<Vector3_Order<int>, Matz>>>>>> exx_IJR_cplx;
        //! exact-exchange Hamiltonian in the basis of KS states, dimension (nspins, n_kpoints, n_bands, n_bands)
        std::map<int, std::map<int, Matz>> exx_KS;
        //! Descriptor for distributed KS exact-exchange matrix blocks.
        ArrayDesc desc_exx_KS;
        //! exact-exchange energy of each state, dimension (nspins, n_kpoints, n_bands). This is actually the diagonal elements of Heex_KS.
        std::map<int, std::map<int, std::map<int, double>>> Eexx;

        Exx(const MeanField& mf_in,
            const AtomicBasis &atbasis_wfc_in,
            const PeriodicBoundaryData &pbc_in,
            const SymmetryContext &symmetry_context_in,
            const KPointBlacsParallelContext &kblacs_ctxt_in,
            const KPointBlacsParallelContext &band_kblacs_ctxt_in,
            const ArrayDesc &desc_wfc_in,
            const ArrayDesc &desc_band_wfc_in,
            bool is_mf_eigvec_k_distributed,
            bool use_symmetry_context_in = true);

        //! Build and store the real-space exchange matrix
        void build(const LibrpaParallelRouting routing,
                   const AtomicBasis &atbasis_abf, const Cs_LRI &Cs,
                   const atpair_R_mat_t& coul_mat);

        void build_KS_kgrid();
        // void build_KS0_kgrid();
        void build_KS_band(const std::map<int, std::map<int, std::map<int, ComplexMatrix>>> &wfc_band,
                           const std::vector<Vector3_Order<double>> &kfrac_band,
                           const AtomPairBvKRemap<atom_t> &bvk_remap);
        void build_KS_kgrid_blacs(const BlacsCtxtHandler &blacs_ctxt_h,
                                  bool use_gpu_replace_scalapack = false);
        // void build_KS0_kgrid_blacs();
        void build_KS_band_blacs(const std::map<int, std::map<int, std::map<int, ComplexMatrix>>> &wfc_band,
                                 const std::vector<Vector3_Order<double>> &kfrac_band,
                                 const AtomPairBvKRemap<atom_t> &bvk_remap,
                                 const BlacsCtxtHandler &blacs_ctxt_h,
                                 bool use_gpu_replace_scalapack = false,
                                 int band_index = 0);
        void write_exx_matrices_KS_binary(const std::string &output_dir,
                                          const std::string &source,
                                          int istate_start,
                                          int istate_end) const;
        void reset_rspace();
        void reset_kspace();
};

} /* end of namespace librpa_int */
