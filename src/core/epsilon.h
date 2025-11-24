#pragma once
#include "chi0.h"
#include "ri.h"
#include "atoms.h"
#include "../math/matrix_m.h"
#include "../mpi/base_blacs.h"

struct CorrEnergy
{
    enum type { RPA, MP2 };
    type etype;
    cplxdb value;
    map<Vector3_Order<double>, cplxdb> qcontrib;
};

CorrEnergy compute_RPA_correlation(const Chi0 &chi0, const atpair_k_cplx_mat_t &coulmat);

CorrEnergy compute_RPA_correlation_blacs(const Chi0 &chi0, const atpair_k_cplx_mat_t &coulmat);
CorrEnergy compute_RPA_correlation_blacs_2d( Chi0 &chi0,  atpair_k_cplx_mat_t &coulmat);
CorrEnergy compute_RPA_correlation_blacs_2d_gamma_only( Chi0 &chi0,  atpair_k_cplx_mat_t &coulmat);
CorrEnergy compute_MP2_correlation(const Chi0 &chi0, const atpair_k_cplx_mat_t &coulmat);

map<double, map<Vector3_Order<double>, atom_mapping<ComplexMatrix>::pair_t_old>> compute_Pi_q(const Chi0 &chi0, const atpair_k_cplx_mat_t &coulmat);
map<double, map<Vector3_Order<double>, atom_mapping<ComplexMatrix>::pair_t_old>> compute_Pi_q_MPI(const Chi0 &chi0, const atpair_k_cplx_mat_t &coulmat);
atom_mapping<ComplexMatrix>::pair_t_old gather_vq_row_q(const int &I, const atpair_k_cplx_mat_t &coulmat, const Vector3_Order<double> &ik_vec);

map<double, std::map<Vector3_Order<double>, Matz>>
compute_Wc_freq_q(Chi0 &chi0, const atpair_k_cplx_mat_t &coulmat_eps,
                  atpair_k_cplx_mat_t &coulmat_wc,
                  const std::vector<cplxdb> &epsmac_LF_imagfreq);

map<double, std::map<Vector3_Order<double>, Matz>>
compute_Wc_freq_q_blacs(Chi0 &chi0, const atpair_k_cplx_mat_t &coulmat_eps,
                        atpair_k_cplx_mat_t &coulmat_wc,
                        const std::vector<cplxdb> &epsilon_mac_imagfreq,
                        const librpa_int::ArrayDesc &ad);

//! Fourier transform screened Coulomb in q-space to R-space, but still in frequency domain
std::map<double, std::map<Vector3_Order<int>, Matz>>
FT_Wc_freq_q(
    std::map<double, std::map<Vector3_Order<double>, Matz>> &Wc_freq_q,
    const int &n_kpoints, const std::vector<Vector3_Order<int>> &Rlist,
    bool remove_freq_q = true, MAJOR major_out = MAJOR::AUTO);

std::map<double, std::map<Vector3_Order<int>, Matz>>
CT_FT_Wc_freq_q(
    std::map<double, std::map<Vector3_Order<double>, Matz>> &Wc_freq_q,
    const TFGrids &tfg, const int &n_kpoints, const std::vector<Vector3_Order<int>> &Rlist,
    bool remove_freq_q = true, MAJOR major_out = MAJOR::AUTO);

ComplexMatrix compute_Pi_freq_q_row_ri(const Vector3_Order<double> &ik_vec, const atom_mapping<ComplexMatrix>::pair_t_old &chi0_freq_q, const atpair_k_cplx_mat_t &Vq_loc, const int &I, const Vector3_Order<double> &q);
ComplexMatrix compute_Pi_freq_q_row(const Vector3_Order<double> &ik_vec, const atom_mapping<ComplexMatrix>::pair_t_old &chi0_freq_q, const atom_mapping<ComplexMatrix>::pair_t_old &Vq_row, const int &I);

cplxdb compute_pi_det_blacs(ComplexMatrix &loc_piT, const librpa_int::ArrayDesc& arrdesc_pi, int *ipiv, int &info);

cplxdb compute_pi_det_blacs_2d(Matz &loc_piT, const librpa_int::ArrayDesc &arrdesc_pi, int *ipiv, int &info);
double compute_pi_det_blacs_2d_gamma_only(matrix_m<double> &loc_piT, const librpa_int::ArrayDesc &arrdesc_pi, int *ipiv, int &info);

void test_libcomm_for_system(const atpair_k_cplx_mat_t &coulmat);
