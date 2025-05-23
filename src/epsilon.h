#pragma once
#include "atoms.h"
#include "base_blacs.h"
#include "chi0.h"
#include "dielecmodel.h"
#include "matrix_m.h"
#include "parallel_mpi.h"
#include "ri.h"
struct CorrEnergy
{
    enum type
    {
        RPA,
        MP2
    };
    type etype;
    complex<double> value;
    map<Vector3_Order<double>, complex<double>> qcontrib;
};

CorrEnergy compute_RPA_correlation(const Chi0 &chi0, const atpair_k_cplx_mat_t &coulmat);

CorrEnergy compute_RPA_correlation_blacs(const Chi0 &chi0, const atpair_k_cplx_mat_t &coulmat);
CorrEnergy compute_RPA_correlation_blacs_2d(Chi0 &chi0, atpair_k_cplx_mat_t &coulmat);
CorrEnergy compute_RPA_correlation_blacs_2d_gamma_only(Chi0 &chi0, atpair_k_cplx_mat_t &coulmat);
CorrEnergy compute_MP2_correlation(const Chi0 &chi0, const atpair_k_cplx_mat_t &coulmat);

map<double, map<Vector3_Order<double>, atom_mapping<ComplexMatrix>::pair_t_old>> compute_Pi_q(
    const Chi0 &chi0, const atpair_k_cplx_mat_t &coulmat);
map<double, map<Vector3_Order<double>, atom_mapping<ComplexMatrix>::pair_t_old>> compute_Pi_q_MPI(
    const Chi0 &chi0, const atpair_k_cplx_mat_t &coulmat);
atom_mapping<ComplexMatrix>::pair_t_old gather_vq_row_q(const int &I,
                                                        const atpair_k_cplx_mat_t &coulmat,
                                                        const Vector3_Order<double> &ik_vec);

map<double, atom_mapping<std::map<Vector3_Order<double>, matrix_m<complex<double>>>>::pair_t_old>
compute_Wc_freq_q(Chi0 &chi0, const atpair_k_cplx_mat_t &coulmat_eps,
                  atpair_k_cplx_mat_t &coulmat_wc,
                  const vector<std::complex<double>> &epsmac_LF_imagfreq);

map<double, atom_mapping<std::map<Vector3_Order<double>, matrix_m<complex<double>>>>::pair_t_old>
compute_Wc_freq_q_blacs(Chi0 &chi0, const atpair_k_cplx_mat_t &coulmat_eps,
                        atpair_k_cplx_mat_t &coulmat_wc,
                        const vector<std::complex<double>> &epsilon_mac_imagfreq);

//! Fourier transform screened Coulomb in q-space to R-space, but still in frequency domain
map<double, atom_mapping<std::map<Vector3_Order<int>, matrix_m<complex<double>>>>::pair_t_old>
FT_Wc_freq_q(
    const map<double,
              atom_mapping<std::map<Vector3_Order<double>, matrix_m<complex<double>>>>::pair_t_old>
        &Wc_freq_q,
    const TFGrids &tfg, const int &n_kpoints, const vector<Vector3_Order<int>> &Rlist);

map<double, atom_mapping<std::map<Vector3_Order<int>, matrix_m<complex<double>>>>::pair_t_old>
CT_FT_Wc_freq_q(
    const map<double,
              atom_mapping<std::map<Vector3_Order<double>, matrix_m<complex<double>>>>::pair_t_old>
        &Wc_freq_q,
    const TFGrids &tfg, const int &n_kpoints, const vector<Vector3_Order<int>> &Rlist);

map<double, atom_mapping<std::map<Vector3_Order<double>, matrix_m<complex<double>>>>::pair_t_old>
CT_FT_Wc_freq2time_q(
    const map<double,
              atom_mapping<std::map<Vector3_Order<double>, matrix_m<complex<double>>>>::pair_t_old>
        &Wc_freq_q,
    const TFGrids &tfg, const int &n_kpoints, const vector<Vector3_Order<int>> &Rlist,
    const vector<Vector3_Order<double>> &qlist);

atom_mapping<std::map<Vector3_Order<int>, matrix_m<complex<double>>>>::pair_t_old CT_FT_Wc_tau_R2q(
    const atom_mapping<std::map<Vector3_Order<double>, matrix_m<complex<double>>>>::pair_t_old
        &Wc_tau_q,
    const TFGrids &tfg, const int &n_kpoints, const vector<Vector3_Order<int>> &Rlist,
    const int &itau);

ComplexMatrix compute_Pi_freq_q_row_ri(const Vector3_Order<double> &ik_vec,
                                       const atom_mapping<ComplexMatrix>::pair_t_old &chi0_freq_q,
                                       const atpair_k_cplx_mat_t &Vq_loc, const int &I,
                                       const Vector3_Order<double> &q);
ComplexMatrix compute_Pi_freq_q_row(const Vector3_Order<double> &ik_vec,
                                    const atom_mapping<ComplexMatrix>::pair_t_old &chi0_freq_q,
                                    const atom_mapping<ComplexMatrix>::pair_t_old &Vq_row,
                                    const int &I);

complex<double> compute_pi_det_blacs(ComplexMatrix &loc_piT, const LIBRPA::Array_Desc &arrdesc_pi,
                                     int *ipiv, int &info);

complex<double> compute_pi_det_blacs_2d(matrix_m<complex<double>> &loc_piT,
                                        const LIBRPA::Array_Desc &arrdesc_pi, int *ipiv, int &info);
double compute_pi_det_blacs_2d_gamma_only(matrix_m<double> &loc_piT,
                                          const LIBRPA::Array_Desc &arrdesc_pi, int *ipiv,
                                          int &info);

void test_libcomm_for_system(const atpair_k_cplx_mat_t &coulmat);
