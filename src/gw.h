/*
 * @file gw.h
 * @brief facilities to calculate self-energy operator.
 */
#pragma once
#include "atoms.h"
#include "matrix_m.h"
#include "meanfield.h"
#include "ri.h"
#include "timefreq.h"

namespace LIBRPA
{

class G0W0
{
   private:
    bool is_rspace_built_;
    bool is_kspace_built_;

   public:
    const MeanField &mf;
    const vector<Vector3_Order<double>> &kfrac_list;
    const TFGrids &tfg;
    const Vector3_Order<int> &period_;

    //! frequency-domain reciprocal-space correlation self-energy, indices
    //! [ispin][freq][k][I][J](n_I, n_J)
    // std::map<int, std::map<double, std::map<Vector3_Order<double>,
    // atom_mapping<Matz>::pair_t_old>>> sigc_is_f_k_IJ;

    //! frequency-domain reciprocal-space correlation self-energy, indices
    //! [ispin][isoc1][isoc2][freq][R][I][J](n_I, n_J)
    std::map<
        int,
        std::map<int, std::map<int, std::map<double, std::map<Vector3_Order<int>,
                                                              atom_mapping<Matz>::pair_t_old>>>>>
        sigc_is_f_R_IJ;

    //! correlation self-energy matrix in the basis of KS states, indices [ispin][ik][freq](n_bands,
    //! n_bands)
    std::map<int, std::map<int, std::map<double, Matz>>> sigc_is_ik_f_KS;

    void build_sigc_matrix_KS(
        const std::vector<std::vector<std::vector<ComplexMatrix>>> &wfc_target,
        const std::vector<Vector3_Order<double>> &kfrac_target);

   public:
    // Constructors
    G0W0(const MeanField &mf, const vector<Vector3_Order<double>> &kfrac_list,
         const TFGrids &tfg_in, const Vector3_Order<int> &period);
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

    //! Build the real-space correlation self-energy matrix on imaginary frequencies with space-time
    //! method using LibRI
    template <typename Tdata>
    void build_spacetime(
        const Cs_LRI &LRI_Cs,
        map<double,
            atom_mapping<std::map<Vector3_Order<double>, matrix_m<complex<double>>>>::pair_t_old>
            &Wc_freq_q,
        const vector<Vector3_Order<int>> &Rlist);

    //! build the correlation self-energy matrix in Kohn-Sham basis at the SCF k-points
    void build_sigc_matrix_KS_kgrid();
    void build_sigc_matrix_KS_kgrid0();
    //! build the correlation self-energy matrix in Kohn-Sham basis at the SCF k-points
    void build_sigc_matrix_KS_band(
        const std::vector<std::vector<std::vector<ComplexMatrix>>> &wfc_target,
        const std::vector<Vector3_Order<double>> &kfrac_band);
};

}  // namespace LIBRPA
