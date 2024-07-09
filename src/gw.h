/*
 * @file gw.h
 * @brief facilities to calculate self-energy operator.
 */
#pragma once
#include "ri.h"
#include "timefreq.h"
#include "meanfield.h"
#include "atoms.h"
#include "matrix_m.h"

namespace LIBRPA
{

class G0W0
{
private:
    bool d_sigc_built;
public:
    const MeanField &mf;
    const vector<Vector3_Order<double>>& kfrac_list;
    const TFGrids &tfg;
    //! frequency-domain reciprocal-space correlation self-energy, indices [ispin][freq][k][I][J](n_I, n_J)
    std::map<int, std::map<double, std::map<Vector3_Order<double>, atom_mapping<Matz>::pair_t_old>>> sigc_is_f_k_IJ;
    //! correlation self-energy matrix in the basis of KS states, indices [ispin][ik][freq](n_bands, n_bands)
    std::map<int, std::map<int, std::map<double, Matz>>> sigc_is_ik_f_KS;
public:
    G0W0(const MeanField &mf,
         const vector<Vector3_Order<double>>& kfrac_list,
         const TFGrids &tfg);
    // delete copy/move constructors
    G0W0(const G0W0 &s_g0w0) = delete;
    G0W0(G0W0 &&s_g0w0) = delete;

    // delete assignment copy/move
    G0W0 operator=(const G0W0 &s_g0w0) const = delete;
    G0W0 operator=(G0W0 &&s_g0w0) = delete;

    //! using native tensor contraction
    // void build_spacetime(
    //     const atpair_R_mat_t &LRI_Cs,
    //     const map<double, atom_mapping<std::map<Vector3_Order<double>,
    //                                             matrix_m<complex<double>>>>::pair_t_old> &Wc_freq_q,
    //     const vector<Vector3_Order<int>> &Rlist,
    //     const Vector3_Order<int> &R_period);

    //! using LibRI
    void build_spacetime_LibRI(
        const Cs_LRI &LRI_Cs,
        const map<double, atom_mapping<std::map<Vector3_Order<double>,
                                                matrix_m<complex<double>>>>::pair_t_old> &Wc_freq_q,
        const vector<Vector3_Order<int>> &Rlist,
        const Vector3_Order<int> &R_period);

    //! build the correlation self-energy matrix in Kohn-Sham basis
    void build_sigc_matrix_KS();
};

}
