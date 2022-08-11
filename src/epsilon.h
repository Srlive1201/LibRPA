#include "chi0.h"
#include "ri.h"
#include "atoms.h"

struct CorrEnergy
{
    enum type { RPA, MP2 };
    type etype;
    complex<double> value;
    map<Vector3_Order<double>, complex<double>> qcontrib;
};

CorrEnergy compute_RPA_correlation(const Chi0 &chi0, const atpair_k_cplx_mat_t &coulmat);

CorrEnergy compute_RPA_correlation_blacs(const Chi0 &chi0, const atpair_k_cplx_mat_t &coulmat);
CorrEnergy compute_MP2_correlation(const Chi0 &chi0, const atpair_k_cplx_mat_t &coulmat);

map<double, map<Vector3_Order<double>, atom_mapping<ComplexMatrix>::pair_t_old>> compute_Pi_q(const Chi0 &chi0, const atpair_k_cplx_mat_t &coulmat);

map<double, atpair_k_cplx_mat_t> compute_Wc_freq_q(const Chi0 &chi0, const atpair_k_cplx_mat_t &coulmat_eps, atpair_k_cplx_mat_t &coulmat_wc);

map<double, atpair_R_cplx_mat_t> CT_FT_Wc_freq_q(const map<double, atpair_k_cplx_mat_t> &Wc_freq_q,
                                                 const TFGrids &tfg, vector<Vector3_Order<int>> Rlist);

atpair_R_cplx_mat_t FT_Vq(const atpair_k_cplx_mat_t &coulmat, vector<Vector3_Order<int>> Rlist);

atom_mapping<ComplexMatrix>::pair_t_old compute_Pi_freq_q(const Vector3_Order<double> &ik_vec, const atom_mapping<ComplexMatrix>::pair_t_old &chi0_freq_q, const atpair_k_cplx_mat_t &coulmat);

complex<double> compute_pi_det(map<size_t, map<size_t, ComplexMatrix>> &pi_freq_q, bool out_pi=0);
