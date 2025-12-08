/*!
 @file chi0.h
 @brief Utlities to compute the independent response function
 */
#pragma once
#include <set>
#include <vector>

#include "../math/vector3_order.h"
#include "atom.h"
#include "atomic_basis.h"
#include "meanfield.h"
#include "pbc.h"
#include "ri.h"
#include "timefreq.h"

namespace librpa_int {

using std::vector;

//! Object to handle calculation of independent repsonse function (\f$\chi_0\f$)
class Chi0
{
    private:
        bool is_mf_eigvec_k_distributed_;
        size_t gf_save;
        size_t gf_discard;
        //! space-time Green's function in occupied space, [ispin][I][J][R][tau]
        /*!
         * @note: tau (index) less than zero correspond to occupied GF,
         *        and larger than zero correspond to unoccpued GF.
         * @note: May need to use ComplexMatrix for GF.
         */
        map<int, atom_mapping<map<Vector3_Order<int>, map<double, matrix>>>::pair_t_old> gf_is_R_tau;

        //! R on which the space-time GF are created, used for atom-pair and rtau routings
        vector<Vector3_Order<int>> Rlist_gf;

        //! Indices of G_{IJ}(R) to build, local to process, used for LibRI routing
        std::vector<std::pair<atpair_t, Vector3_Order<int>>> IJRs_gf_local;

        //! chi0 data in frequency domain and reciprocal space, [omega][q]
        map<double, map<Vector3_Order<double>, atom_mapping<ComplexMatrix>::pair_t_old>> chi0_q;

        void build_gf_Rt(Vector3_Order<int> R, double tau);

        // Free the intermeidate Green's functions
        void free_gf_Rt();

        //! Internal procedure to compute chi0_q by space-time method
        /*
         * @todo add threshold parameter. Maybe in the class level?
         */
        void build_chi0_q_space_time(const LibrpaParallelRouting routing, const Cs_LRI &Cs,
                                     const vector<atpair_t> &atpairs_ABF);

        // NOTE: the following three methods could be converted to static functions in chi0.cpp
        void build_chi0_q_space_time_atom_pair_routing(const Cs_LRI &Cs,
                                                       const vector<atpair_t> &atpairs_ABF);
        void build_chi0_q_space_time_R_tau_routing(const Cs_LRI &Cs,
                                                   const vector<atpair_t> &atpairs_ABF);
        void build_chi0_q_space_time_LibRI_routing(const Cs_LRI &Cs,
                                                   const vector<atpair_t> &atpairs_ABF);

        //! Internal procedure to compute chi0_q in the conventional method, i.e. in frequency domain and reciprocal space
        // TODO: implement the conventional method
        void build_chi0_q_conventional(const Cs_LRI &Cs,
                                       const vector<atpair_t> &atpairs_ABF);
        /*!
         * s_alpha and s_beta are the spin component of unoccupied Green's function, G_{alpha, beta}(tau)
         * correspondingly, occupied GF G_{beta, alpha}(-tau) will be used. itau must be positive.
         */
        matrix compute_chi0_s_munu_tau_R(const atpair_R_mat_t &Cs_IJR,
                                         int spin_channel,
                                         atom_t mu, atom_t nu, double tau, Vector3_Order<int> R);
        // copy some reshape method inside chi0, test performance
        /* matrix reshape_Cs(const size_t n1, const size_t n2, const size_t n3, const std::shared_ptr<matrix> &Cs); */
        /* matrix reshape_dim_Cs(const size_t n1, const size_t n2, const size_t n3, const std::shared_ptr<matrix> &Cs);//(n1*n2,n3) -> (n1,n2*n3) */
        /* matrix reshape_mat(const size_t n1, const size_t n2, const size_t n3, const matrix &mat); */
        /* matrix reshape_mat_21(const size_t n1, const size_t n2, const size_t n3, const matrix &mat); //(n1,n2*n3) -> (n1*n2,n3) */
    public:
        const MeanField &mf;
        const AtomicBasis &atbasis_wfc;
        const AtomicBasis &atbasis_abf;
        const PeriodicBoundaryData &pbc;
        const TFGrids &tfg;
        const MpiCommHandler &comm_h;

        double libri_threshold_C;
        double libri_threshold_G;
        double gf_threshold;

        Chi0(const MeanField &mf_in, const AtomicBasis &atbasis_wfc_in,
             const AtomicBasis &atbasis_abf_in, const PeriodicBoundaryData &pbc_in,
             const TFGrids &tfg_in, const MpiCommHandler &comm_h_in, bool is_mf_eigvec_k_distributed);
        ~Chi0() {};
        //! Build the independent response function in q-omega domain for ABFs on the atom pairs atpair_ABF and q-vectors in qlist
        void build(LibrpaParallelRouting routing,
                   const Cs_LRI &Cs,
                   const vector<atpair_t> &atpair_ABF);
        const map<double, map<Vector3_Order<double>, atom_mapping<ComplexMatrix>::pair_t_old>> & get_chi0_q() const { return chi0_q; }
        void free_chi0_q(const double freq, const Vector3_Order<double> q);
};

}
