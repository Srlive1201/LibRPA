#pragma once
#include <functional>
#include <vector>
#include "atomic_basis.h"
#include "../math/matrix3.h"

namespace librpa_int {

//! double-dispersion Havriliak-Negami model
struct DoubleHavriliakNegami
{
    static const int d_npar;
    static const std::function<double(double, const std::vector<double> &)> func_imfreq;
    static const std::function<void(std::vector<double> &, double, const std::vector<double> &)>
        grad_imfreq;
};

std::vector<double> interpolate_dielec_func(int option, const std::vector<double> &frequencies_in,
                                            const std::vector<double> &df_in,
                                            const std::vector<double> &frequencies_target);
}

#include "../meanfield.h"
#include "ri.h"

class diele_func
{
   private:
    // ( alpha, beta, omega )
    std::vector<std::vector<std::vector<std::complex<double>>>> head;
    // ( alpha, lambda:n_abfs-n_singular-1, omega )
    std::vector<std::vector<std::vector<std::complex<double>>>> wing;
    // ( lambda: n_abfs-n_singular-1, mu: n_abfs)
    std::vector<std::vector<std::complex<double>>> Coul_vector;
    // ( lambda: n_abfs-n_singular-1 )
    std::vector<std::complex<double>> Coul_value;
    // ( mu: n_abfs, m: n_bands, n: n_bands, k )
    std::vector<std::vector<std::vector<std::map<Vector3_Order<double>, std::complex<double>>>>>
        Ctri_mn;
    // ( mu: n_abfs, i: i atom basis, j: j atom basis, k, I atom, J atom, R cell  )
    // Ctri_ij.data_libri[I][{J, k_array}](mu, i, j)
    librpa_int::Cs_LRI_clx Ctri_ij;

    const MeanField &meanfield_df;
    const std::vector<double> &omega;
    const std::vector<Vector3_Order<double>> &kfrac_band;
    const int n_basis, n_states, n_spin, n_abf;
    const librpa_int::Matrix3 &latvec_;
    size_t n_nonsingular;

   public:
    diele_func(const MeanField &mf, const std::vector<Vector3_Order<double>> &kfrac,
               const std::vector<double> &frequencies_target, const int nbasis, const int nstates,
               const int nspin, const int nabf, const librpa_int::Matrix3 &latvec)
        : meanfield_df(mf),
          omega(frequencies_target),
          kfrac_band(kfrac),
          n_basis(nbasis),
          n_states(nstates),
          n_spin(nspin),
          n_abf(nabf),
          latvec_(latvec)
    {};
    ~diele_func() {};
    void init_headwing(const librpa_int::AtomicBasis &atomic_basis_abf, double vq_threshold,
                       const librpa_int::atpair_k_cplx_mat_t &Vq_cut);
    void init_Cs(const librpa_int::Cs_LRI &Cs_data);
    // All calculation in unit: Ang and eV.
    void cal_head();
    double cal_factor(std::string name);
    void test_head();

    void cal_wing(const librpa_int::Cs_LRI &Cs_data, const librpa_int::AtomicBasis &atomic_basis_wfc);  // atpair_k_cplx_mat_t &Vq_cut, Cs_LRI &Cs_data
    // tranform Cs_ij(R) to Cs_ij(k)
    void FT_R2k(const librpa_int::Cs_LRI &Cs_data);
    void Cs_ij2mn(const librpa_int::AtomicBasis &atomic_basis_wfc);
    void get_Xv(const librpa_int::AtomicBasis &atomic_basis_abf, double vq_threshold,
                const librpa_int::atpair_k_cplx_mat_t &Vq_cut);  // diagonalize Vq_cut(q=0)
    void test_wing();
};
