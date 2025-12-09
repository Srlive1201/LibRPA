/*
 * @file exx.h
 * @brief utilities for computing exact exchange energies, including orbital and total energies.
 */
#include "../math/matrix_m.h"
#include "../mpi/base_blacs.h"
#include "atomic_basis.h"
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

        void build_dmat_R(const Vector3_Order<int>& R);
        void build_dmat_R(const atom_t& I, const atom_t& J, const Vector3_Order<int>& R);
        void build_LibRI(const Cs_LRI &Cs,
                         const vector<Vector3_Order<int>> &Rlist,
                         const atpair_R_mat_t& coul_mat);

        void build_KS(const std::map<int, std::map<int, ComplexMatrix>> &wfc_target,
                      const std::vector<Vector3_Order<double>> &kfrac_target,
                      const Atoms &geometry);
        void build_KS_blacs(const std::map<int, std::map<int, ComplexMatrix>> &wfc_target,
                            const std::vector<Vector3_Order<double>> &kfrac_target,
                            const Atoms &geometry,
                            const BlacsCtxtHandler &blacs_ctxt_h);
    public:
        //! refenrence to the MeanField object to compute density matrix
        const MeanField& mf;
        const AtomicBasis &atbasis_wfc;
        const PeriodicBoundaryData &pbc;
        const MpiCommHandler &comm_h;

        double libri_threshold_C;
        double libri_threshold_V;
        double libri_threshold_D;

        //! Density matrix in lattice vector space, dimension (nspins, I, J, R, nao_I, nao_J)
        map<int, atpair_R_mat_t> dmat;
        //! exact-exchange Hamiltonian in real space, dimension (nspins, R, I, J, nao_I, nao_J)
        map<int, map<Vector3_Order<int>, map<atom_t, map<atom_t, Matd>>>> exx;
        //! exact-exchange Hamiltonian in the basis of KS states, dimension (nspins, n_kpoints, n_bands, n_bands)
        map<int, map<int, Matz>> exx_is_ik_KS;
        //! exact-exchange energy of each state, dimension (nspins, n_kpoints, n_bands). This is actually the diagonal elements of Heex_KS.
        map<int, map<int, map<int, double>>> Eexx;

        Exx(const MeanField& mf_in,
            const AtomicBasis &atbasis_wfc_in,
            const PeriodicBoundaryData &pbc_in,
            const MpiCommHandler &comm_h_in,
            bool is_mf_eigvec_k_distributed);

        //! Build and store the real-space exchange matrix
        void build(const LibrpaParallelRouting routing,
                   const AtomicBasis &atbasis_abf, const Cs_LRI &Cs,
                   const atpair_R_mat_t& coul_mat);

        void build_KS_kgrid();
        void build_KS_band(const std::map<int, std::map<int, ComplexMatrix>> &wfc_band,
                           const std::vector<Vector3_Order<double>> &kfrac_band,
                           const Atoms &geometry);
        void build_KS_kgrid_blacs(const BlacsCtxtHandler &blacs_ctxt_h);
        void build_KS_band_blacs(const std::map<int, std::map<int, ComplexMatrix>> &wfc_band,
                                 const std::vector<Vector3_Order<double>> &kfrac_band,
                                 const Atoms &geometry,
                                 const BlacsCtxtHandler &blacs_ctxt_h);
        void reset_rspace();
        void reset_kspace();
};

} /* end of namespace librpa_int */
