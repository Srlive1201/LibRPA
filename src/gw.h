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
public:
    const MeanField &mf;
    const vector<Vector3_Order<double>>& kfrac_list;
    const TFGrids &tfg;
    //! frequency-domain reciprocal-space correlation self-energy, indices [ispin][freq][k][I][J](n_I, n_J)
    std::map<int, std::map<double, std::map<Vector3_Order<double>, atom_mapping<matrix_m<std::complex<double>>>::pair_t_old>>> sigc_is_freq_k_IJ;
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

    //!
    void build_spacetime(
        const atpair_R_mat_t &LRI_Cs,
        const map<double, atom_mapping<std::map<Vector3_Order<double>,
                                                matrix_m<complex<double>>>>::pair_t_old> &Wc_freq_q,
        const vector<Vector3_Order<int>> &Rlist,
        const Vector3_Order<int> &R_period);
};

}
