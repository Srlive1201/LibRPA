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
        unsigned gf_save;
        unsigned gf_discard;
        // char gf_R_timefreq;
        // char gf_k_timefreq;
        // char chi_R_timefreq;
        // char chi_k_timefreq;
        // [is][I][J][iR][itime](iwt1, iwt2)
        //! space-time Green's function in occupied space, [ispin][itau/iomega][iR][I][J]
        //! @note: should use complex matrix in general?
        map<size_t, map<size_t, map<size_t, atom_mapping<matrix>::pair_t_old>>> gf_occ_Rt;
        //! space-time Green's function in unoccupied (i.e. virtual) space
        map<size_t, map<size_t, map<size_t, atom_mapping<matrix>::pair_t_old>>> gf_unocc_Rt;
        /* //! chi0 data in real-space, [itau/iomega][iR] */
        /* map<int, map<int, atom_mapping<matrix>::pair_t_old>> chi0_R; */
        //! chi0 data in reci-space, [itau/iomega][ik]
        map<int, map<int, atom_mapping<ComplexMatrix>::pair_t_old>> chi0_q;
        void build_gf_Rt(size_t iR, Vector3_Order<int> R, size_t itau, char ov);
        //! Internal procedure to compute chi0_q by space-time method
        /*
         @todo add threshold parameter. Maybe in the class level?
         */
        void build_chi0_q_space_time(const atpair_R_mat_t &LRI_Cs,
                                     const vector<Vector3_Order<int>> &Rlist,
                                     const Vector3_Order<int> &R_period,
                                     const vector<atpair_t> &atpair_ABF,
                                     const vector<Vector3_Order<double>> &qlist);
        //! Internal procedure to compute chi0_q in the conventional method, i.e. in frequency domain and reciprocal space
        // TODO: implement the conventional method
        void build_chi0_q_conventional(const atpair_R_mat_t &LRI_Cs,
                                       const vector<Vector3_Order<int>> &Rlist,
                                       const Vector3_Order<int> &R_period,
                                       const vector<atpair_t> &atpair_ABF,
                                       const vector<Vector3_Order<double>> &qlist);
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
        //! Compute the independent response function in q-omega domain for ABFs on the atom pairs atpair_ABF and q-vectors in qlist
        void compute(const atpair_R_mat_t &LRI_Cs,
                     const vector<Vector3_Order<int>> &Rlist,
                     const Vector3_Order<int> &R_period,
                     const vector<atpair_t> &atpair_ABF,
                     const vector<Vector3_Order<double>> &qlist,
                     TFGrids::GRID_TYPES gt, bool use_space_time);
};

//! Compute the real-space independent reponse function in space-time method
/*!
 @param[in] gf_occ_ab_Rt: occupied Green's function at time tau in a particular spin alpha-beta channel
 @param[in] gf_unocc_ab_Rt: same as above, but for unoccupied Green's function
 @param[in] Cs_mu: RI coefficient with atom ABF mu
 @param[in] Cs_nu: RI coefficient onth atom ABF nu
 @param[in] Rlist: the integer list of unit cell coordinates
 @param[in] R_period: the periodicity of super cell
 @param[in] mu, nu: indices of atoms with ABF
 @param[in] iR: unit cell index of chi0 to compute
 @param[in] tau: time (imaginary) of chi0 to compute
 */
matrix compute_chi0_munu_iR_tau(const map<size_t, atom_mapping<matrix>::pair_t_old> &gf_occ_ab_Rt,
                                const map<size_t, atom_mapping<matrix>::pair_t_old> &gf_unocc_ab_Rt,
                                const atpair_R_mat_t &Cs_mu, const atpair_R_mat_t &Cs_nu,
                                const vector<Vector3_Order<int>> &Rlist, const Vector3_Order<int> &R_period, 
                                atom_t mu, atom_t nu, int iR, double tau);
