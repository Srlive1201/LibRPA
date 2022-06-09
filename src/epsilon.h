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

CorrEnergy compute_MP2_correlation(const Chi0 &chi0, const atpair_k_cplx_mat_t &coulmat);

map<double, map<Vector3_Order<double>, atom_mapping<ComplexMatrix>::pair_t_old>> compute_Pi_q(const Chi0 &chi0, const atpair_k_cplx_mat_t &coulmat);

map<double, map<Vector3_Order<double>, atom_mapping<ComplexMatrix>::pair_t_old>> compute_Wc_freq_q(const Chi0 &chi0, const atpair_k_cplx_mat_t &coulmat_eps, const atpair_k_cplx_mat_t &coulmat_wc);

map<double, map<Vector3_Order<int>, atom_mapping<ComplexMatrix>::pair_t_old>>
Wc_tau_R_by_ct_ft(const map<double, map<Vector3_Order<int>, atom_mapping<ComplexMatrix>::pair_t_old>> &Wc_freq_q,
                  TFGrids tfg, vector<Vector3_Order<int>> Rlist);
