/*!
 @file chi0.h
 @brief Utlities to compute the independent response function
 @warning Not work now
 */
#pragma once
#include <vector>
#include <set>
#include "meanfield.h"
#include "timefreq.h"
#include "ri.h"
#include "vector3_order.h"

using std::vector;

//! Object to handle calculation of independent repsonse function (\f$\chi_0\f$)
class Chi0
{
    private:
        size_t gf_save;
        size_t gf_discard;
        //! space-time Green's function in occupied space, [ispin][I][J][R][tau]
        /*!
         * @note: tau (index) less than zero correspond to occupied GF,
         *        and larger than zero correspond to unoccpued GF.
         * @note: May need to use ComplexMatrix for GF.
         */
        map<int, atom_mapping<map<Vector3_Order<int>, map<double, matrix>>>::pair_t_old> gf_is_R_tau;
        //! R on which the space-time GF are created.
        vector<Vector3_Order<int>> Rlist_gf;
        //! chi0 data in frequency domain and reciprocal space, [omega][q]
        map<double, map<Vector3_Order<double>, atom_mapping<ComplexMatrix>::pair_t_old>> chi0_q;
        void build_gf_Rt(Vector3_Order<int> R, double tau);
        //! Internal procedure to compute chi0_q by space-time method
        /*
         @todo add threshold parameter. Maybe in the class level?
         */
        void build_chi0_q_space_time(const atpair_R_mat_t &LRI_Cs,
                                     const Vector3_Order<int> &R_period,
                                     const vector<atpair_t> &atpairs_ABF,
                                     const vector<Vector3_Order<double>> &qlist);
        void build_chi0_q_space_time_atom_pair_routing(const atpair_R_mat_t &LRI_Cs,
                                                       const Vector3_Order<int> &R_period,
                                                       const vector<atpair_t> &atpairs_ABF,
                                                       const vector<Vector3_Order<double>> &qlist);
        void build_chi0_q_space_time_R_tau_routing(const atpair_R_mat_t &LRI_Cs,
                                                   const Vector3_Order<int> &R_period,
                                                   const vector<atpair_t> &atpairs_ABF,
                                                   const vector<Vector3_Order<double>> &qlist);
#ifdef __USE_LIBRI
        void build_chi0_q_space_time_LibRI_routing(const atpair_R_mat_t &LRI_Cs,
                                                   const Vector3_Order<int> &R_period,
                                                   const vector<atpair_t> &atpairs_ABF,
                                                   const vector<Vector3_Order<double>> &qlist);
#endif

        //! Internal procedure to compute chi0_q in the conventional method, i.e. in frequency domain and reciprocal space
        // TODO: implement the conventional method
        void build_chi0_q_conventional(const atpair_R_mat_t &LRI_Cs,
                                       const Vector3_Order<int> &R_period,
                                       const vector<atpair_t> &atpairs_ABF,
                                       const vector<Vector3_Order<double>> &qlist);
        /*!
         * s_alpha and s_beta are the spin component of unoccupied Green's function, G_{alpha, beta}(tau)
         * correspondingly, occupied GF G_{beta, alpha}(-tau) will be used. itau must be positive.
         */
        matrix compute_chi0_s_munu_tau_R(const atpair_R_mat_t &LRI_Cs,
                                          const Vector3_Order<int> &R_period, 
                                          int spin_channel,
                                          atom_t mu, atom_t nu, double tau, Vector3_Order<int> R);
        // copy some reshape method inside chi0, test performance
        /* matrix reshape_Cs(const size_t n1, const size_t n2, const size_t n3, const std::shared_ptr<matrix> &Cs); */
        /* matrix reshape_dim_Cs(const size_t n1, const size_t n2, const size_t n3, const std::shared_ptr<matrix> &Cs);//(n1*n2,n3) -> (n1,n2*n3) */
        /* matrix reshape_mat(const size_t n1, const size_t n2, const size_t n3, const matrix &mat); */
        /* matrix reshape_mat_21(const size_t n1, const size_t n2, const size_t n3, const matrix &mat); //(n1,n2*n3) -> (n1*n2,n3) */
    public:
        double gf_R_threshold;
        MeanField mf;
        const vector<Vector3_Order<double>> &klist;
        TFGrids tfg;
        Chi0(const MeanField &mf_in,
             const vector<Vector3_Order<double>> &klist_in,
             unsigned n_tf_grids):
            mf(mf_in), klist(klist_in), tfg(n_tf_grids) { gf_R_threshold = 1e-9; }
        ~Chi0() {};
        //! Build the independent response function in q-omega domain for ABFs on the atom pairs atpair_ABF and q-vectors in qlist
        void build(const atpair_R_mat_t &LRI_Cs,
                   const vector<Vector3_Order<int>> &Rlist,
                   const Vector3_Order<int> &R_period,
                   const vector<atpair_t> &atpair_ABF,
                   const vector<Vector3_Order<double>> &qlist,
                   TFGrids::GRID_TYPES gt, bool use_space_time);
        const map<double, map<Vector3_Order<double>, atom_mapping<ComplexMatrix>::pair_t_old>> & get_chi0_q() const { return chi0_q; }
        
};

//! Compute the real-space independent reponse function in space-time method on a particular time
/*!
 @param[in] gf_occ_ab_t: occupied Green's function at time tau in a particular spin alpha-beta channel
 @param[in] gf_unocc_ab_t: same as above, but for unoccupied Green's function
 @param[in] LRI_Cs: LRI coefficients
 @param[in] Rlist: the integer list of unit cell coordinates
 @param[in] R_period: the periodicity of super cell
 @param[in] iRs: indices of R to compute
 @param[in] mu, nu: indices of atoms with ABF
 @retval mapping from unit cell index to matrix
 */
map<size_t, matrix> compute_chi0_munu_tau_LRI_saveN_noreshape(const map<size_t, atom_mapping<matrix>::pair_t_old> &gf_occ_ab_t,
                                                              const map<size_t, atom_mapping<matrix>::pair_t_old> &gf_unocc_ab_t,
                                                              const atpair_R_mat_t &LRI_Cs,
                                                              const vector<Vector3_Order<int>> &Rlist, const Vector3_Order<int> &R_period, 
                                                              const vector<int> iRs,
                                                              atom_t mu, atom_t nu);

