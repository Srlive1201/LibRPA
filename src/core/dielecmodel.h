#pragma once
#include <functional>
#include <vector>
#include "atomic_basis.h"
#include "../math/matrix3.h"
#include "../math/matrix_m.h"
#include "../mpi/base_blacs.h"
#include "../meanfield.h"
#include "pbc.h"
#include "ri.h"

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

// All calculation in unit: Bohr and Ha.
class diele_func
{
   private:
    // ( omega, alpha * beta  )
    std::vector<matrix_m<std::complex<double>>> head;
    // ( omega, mu:n_abfs * alpha ) for comparison with FHI-aims
    std::vector<matrix_m<std::complex<double>>> wing_mu;
    // ( omega, lambda:n_lambda * alpha)
    std::vector<matrix_m<std::complex<double>>> wing;
    // std::vector<std::vector<std::vector<std::complex<double>>>> wing;

    // ( i:n_lambda, j:n_lambda )
    matrix_m<std::complex<double>> body_inv;
    // ( i:3, j:3 )
    matrix_m<std::complex<double>> Lind;
    // ( i:n_lambda, j:3 )
    matrix_m<std::complex<double>> bw;
    // ( i:3, j:n_lambda )
    matrix_m<std::complex<double>> wb;
    // ( i:n_lambda, j:n_lambda )
    matrix_m<std::complex<double>> chi0;
    // ( lambda: n_nonsingular-1, mu: n_abfs)
    std::vector<std::vector<std::complex<double>>> Coul_vector;
    // ( lambda: n_nonsingular-1 )
    std::vector<std::complex<double>> Coul_value;
    // ( mu: n_abfs, m: n_bands, n: n_bands, k )
    std::vector<std::vector<std::vector<std::map<Vector3_Order<double>, std::complex<double>>>>>
        Ctri_mn;
    // ( mu: n_abfs@I, i: i atom basis, j: j atom basis, k, I atom, J atom, q cell  )
    // Ctri_ij.data_libri[I][{J, k_array}](mu, i, j)
    Cs_LRI_clx Ctri_ij;
    // ( mu: n_abfs@I, i: i atom basis, j: j atom basis, k, I atom, J atom, R cell  )
    // Ctri_ij.data_libri[I][{J, R}](mu, i, j)
    // used for reduce all mpi Cs_data to Cs_IJR
    Cs_LRI Cs_IJR;

    MeanField &meanfield_df;
    const std::vector<double> &omega;
    const std::vector<Vector3_Order<double>> &kfrac_band;
    int n_basis, n_states, n_spin, n_abf, nk;
    const librpa_int::PeriodicBoundaryData &pbc_;
    const librpa_int::AtomicBasis &atomic_basis_wfc_;
    const librpa_int::AtomicBasis &atomic_basis_abf_;
    size_t n_nonsingular;
    // lebedev-quadrature, qw has absorbed 4Pi.
    std::vector<double> qx_leb, qy_leb, qz_leb, qw_leb;
    // gamma reciprocal lattice vector, (27-1)*3
    std::vector<Vector3_Order<double>> g_enclosing_gamma;
    std::vector<double> q_gamma;
    double vol_gamma;

   public:
    diele_func(MeanField &mf, const std::vector<Vector3_Order<double>> &kfrac,
               const librpa_int::AtomicBasis &atomic_basis_wfc,
               const librpa_int::AtomicBasis &atomic_basis_abf,
               const std::vector<double> &frequencies_target, const int nbasis, const int nstates,
               const int nspin, const int nabf, const librpa_int::PeriodicBoundaryData &pbc)
        : meanfield_df(mf),
          omega(frequencies_target),
          kfrac_band(kfrac),
          n_basis(nbasis),
          n_states(nstates),
          n_spin(nspin),
          n_abf(nabf),
          pbc_(pbc),
          atomic_basis_wfc_(atomic_basis_wfc),
          atomic_basis_abf_(atomic_basis_abf)
    {};
    ~diele_func() {};
    void init(double vq_threshold, const librpa_int::atpair_k_cplx_mat_t &Vq);
    void init_Cs(const librpa_int::Cs_LRI &Cs_data);
    void init_wing();
    // All calculation in unit: Ang and eV.
    // void set(MeanField &mf, std::vector<Vector3_Order<double>> &kfrac,
    //          std::vector<double> frequencies_target, int nbasis, int nstates, int nspin);

    void cal_head();
    double cal_factor(std::string name);
    void test_head();
    std::vector<double> get_head_vec();

    void cal_wing(const librpa_int::Cs_LRI &Cs_data);  // atpair_k_cplx_mat_t &Vq, Cs_LRI &Cs_data
    // tranform Cs_ij(R) to Cs_ij(k)
    void FT_R2k(const librpa_int::Cs_LRI &Cs_data);
    void Cs_ij2mn();
    // diagonalize Vq(q=0)
    void get_Xv(double vq_threshold,
                const librpa_int::atpair_k_cplx_mat_t &Vq);  // diagonalize Vq(q=0)
    std::complex<double> compute_wing(int alpha, int iomega, int mu);
    std::complex<double> compute_Cs_ij2mn(int mu, int m, int n, int ik);
    std::complex<double> compute_Cijk(const librpa_int::Cs_LRI &Cs_data, int mu, int I, int i, int J, int j, int ik);
    // compute wing in ABF representation
    // transform wing from ABF to Coulomb representation
    void wing_mu_to_lambda(matrix_m<std::complex<double>> &sqrtveig_blacs,
                           ArrayDesc &desc_nabf_nabf_opt);
    // tranform Cs_ij(R) to Cs_ij(k)
    // diagonalize real Vq_cut(q=0)
    void get_Xv_real(double vq_threshold, const librpa_int::atpair_k_cplx_mat_t &Vq);
    // diagonalize complex Vq_cut(q=0)
    void get_Xv_cpl(double vq_threshold, const librpa_int::atpair_k_cplx_mat_t &Vq);
    void test_wing();
    // set wing=0 for debug
    void set_0_wing();

    ArrayDesc get_body_inv(matrix_m<std::complex<double>> &chi0_block,
                            ArrayDesc &desc_nabf_nabf_opt);
    void construct_L(const int ifreq, ArrayDesc &desc_body);

    // Lebedev-Laikov quadrature
    void get_Leb_points();
    void get_g_enclosing_gamma();
    void calculate_q_gamma();
    void cal_eps(const int ifreq, ArrayDesc &desc_nabf_nabf_opt, ArrayDesc &desc_body);
    // not used now due to performance optimization
    // std::complex<double> compute_chi0_inv_00(const int ifreq);
    // std::complex<double> compute_chi0_inv_ij(const int ifreq, int i, int j);
    void rewrite_eps(matrix_m<std::complex<double>> &chi0_block, const int ifreq,
                     ArrayDesc &desc_nabf_nabf_opt);
    void assign_chi0(matrix_m<std::complex<double>> &chi0_block, ArrayDesc &desc_nabf_nabf_opt);
};

}
