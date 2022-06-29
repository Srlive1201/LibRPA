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

complex<double> compute_pi_det(map<size_t, map<size_t, ComplexMatrix>> &pi_freq_q);