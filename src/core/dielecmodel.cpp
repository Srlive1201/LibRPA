#include "dielecmodel.h"

#include <cassert>
#include <cmath>

#include "../math/fitting.h"
#include "../math/interpolate.h"
#include "../math/vec.h"
#include "../math/matrix_m.h"
#include "../math/utils_matrix_m_mpi.h"
#include "../mpi/base_blacs.h"
#include "../utils/error.h"
#include "../utils/constants.h"
#include "../utils/profiler.h"
#include "../utils/utils_mem.h"
#include "atomic_basis.h"
#include "pbc.h"
#include "ri.h"
#include <lebedev_quadrature.hpp>
#ifdef LIBRPA_USE_LIBRI
#include <RI/comm/mix/Communicate_Tensors_Map_Judge.h>
#include <RI/global/Tensor.h>
#else
#include "../utils/libri_stub.h"
#endif

#include "libri_utils.h"
#include "matrix_m_parallel_utils.h"
using RI::Tensor;
using RI::Communicate_Tensors_Map_Judge::comm_map2_first;

namespace librpa_int {

const int DoubleHavriliakNegami::d_npar = 8;

const std::function<double(double, const std::vector<double> &)>
    DoubleHavriliakNegami::func_imfreq = [](double u, const std::vector<double> &pars)
{
    return 1.0 + (pars[0] - 1.0) / std::pow(1.0 + std::pow(u * pars[3], pars[1]), pars[2]) +
           (pars[4] - 1.0) / std::pow(1.0 + pow(u * pars[7], pars[5]), pars[6]);
};

const std::function<void(std::vector<double> &, double, const std::vector<double> &)>
    DoubleHavriliakNegami::grad_imfreq =
        [](std::vector<double> &grads, double u, const std::vector<double> &pars)
{
    using std::pow;
    using std::log;
    grads[0] = 1.0 / pow(1.0 + pow(u * pars[3], pars[1]), pars[2]);
    grads[1] = (pars[0] - 1.0) * (-pars[2]) / pow(1.0 + pow(u * pars[3], pars[1]), pars[2] + 1) *
               log(u * pars[3]) * pow(u * pars[3], pars[1]);
    grads[2] = (1.0 - pars[0]) * log(1.0 + pow(u * pars[3], pars[1])) /
               pow(1.0 + pow(u * pars[3], pars[1]), pars[2]);
    grads[3] = (pars[0] - 1.0) * (-pars[2]) / pow(1.0 + pow(u * pars[3], pars[1]), pars[2] + 1) *
               pars[1] / pars[3] * pow(u * pars[3], pars[1]);
    grads[4] = 1.0 / pow(1.0 + pow(u * pars[7], pars[5]), pars[6]);
    grads[5] = (pars[4] - 1.0) * (-pars[6]) / pow(1.0 + pow(u * pars[7], pars[5]), pars[6] + 1) *
               log(u * pars[7]) * pow(u * pars[7], pars[5]);
    grads[6] = (1.0 - pars[4]) * log(1.0 + pow(u * pars[7], pars[5])) /
               pow(1.0 + pow(u * pars[7], pars[5]), pars[6]);
    grads[7] = (pars[4] - 1.0) * (-pars[6]) / pow(1.0 + pow(u * pars[7], pars[5]), pars[6] + 1) *
               pars[5] / pars[7] * pow(u * pars[7], pars[5]);
};

std::vector<double> interpolate_dielec_func(int option, const std::vector<double> &frequencies_in,
                                            const std::vector<double> &df_in,
                                            const std::vector<double> &frequencies_target)
{
    using librpa_int::DoubleHavriliakNegami;
    std::vector<double> df_target;

    switch (option)
    {
        case 0: /* No extrapolation, copy the input data to target */
        {
            assert(frequencies_in.size() == frequencies_target.size());
            df_target = df_in;
            break;
        }
        case 1: /* Use spline interpolation */
        {
            df_target = librpa_int::interp_cubic_spline(
                    frequencies_in, df_in, frequencies_target);
            break;
        }
        case 2: /* Use dielectric model for fitting */
        {
            librpa_int::LevMarqFitting levmarq;
            // use double-dispersion Havriliak-Negami model
            // initialize the parameters as 1.0
            std::vector<double> pars(DoubleHavriliakNegami::d_npar, 1);
            pars[0] = pars[4] = df_in[0];
            df_target = levmarq.fit_eval(pars, frequencies_in, df_in,
                                         DoubleHavriliakNegami::func_imfreq,
                                         DoubleHavriliakNegami::grad_imfreq,
                                         frequencies_target);
            break;
        }
        default:
            throw LIBRPA_RUNTIME_ERROR("Unsupported value for option");
    }

    return df_target;
}

// void diele_func::set(MeanField &mf, std::vector<Vector3_Order<double>> &kfrac,
//                      std::vector<double> frequencies_target, int nbasis, int nstates, int nspin)
// {
//     meanfield_df = mf;
//     kfrac_band = kfrac;
//     omega = frequencies_target;
//     n_basis = nbasis;
//     n_states = nstates;
//     n_spin = nspin;
//     init();
// };

void diele_func::init(double vq_threshold, const librpa_int::atpair_k_cplx_mat_t &Vq)
{
    this->n_abf = atomic_basis_abf_.nb_total;
    this->nk = this->kfrac_band.size();
    int n_omega = this->omega.size();
    get_Xv_cpl(vq_threshold, Vq);

    this->head.clear();
    this->head.resize(n_omega);
    for (int iomega = 0; iomega != n_omega; iomega++)
    {
        head[iomega].resize(3, 3, MAJOR::COL);
    }
};

void diele_func::init_wing()
{
    int n_omega = this->omega.size();
    this->n_abf = atomic_basis_abf_.nb_total;
    get_Xv_cpl();
    this->wing_mu.clear();
    this->wing_mu.resize(n_omega);
    this->wing.clear();
    this->wing.resize(n_omega);
    this->Lind.resize(3, 3, MAJOR::COL);
    this->bw.resize(n_nonsingular - 1, 3, MAJOR::COL);
    this->wb.resize(3, n_nonsingular - 1, MAJOR::COL);
    for (int iomega = 0; iomega != n_omega; iomega++)
    {
        wing_mu[iomega].resize(n_abf, 3, MAJOR::COL);
        wing[iomega].resize(n_nonsingular - 1, 3, MAJOR::COL);
    }
    get_Leb_points();
    get_g_enclosing_gamma();
    calculate_q_gamma();
    if (mpi_comm_global_h.is_root())
        std::cout << "* Success: initalize and calculate lebdev points and q_gamma." << std::endl;
}

// intraband term is not considered
void diele_func::cal_head()
{
    using global::profiler;

    const bool use_soc = false;

    profiler.start("cal_head");

    std::complex<double> tmp;
    int nocc = 0;
    double dielectric_unit = cal_factor("head");

    for (int ispin = 0; ispin != n_spin; ispin++)
    {
        auto &wg = this->meanfield_df.get_weight()[ispin];
        auto &eigenvalues = this->meanfield_df.get_eigenvals()[ispin];
        auto &velocity = this->meanfield_df.get_velocity()[ispin];
        for (int ik = 0; ik != nk; ik++)
        {
            for (int iocc = 0; iocc != n_states; iocc++)
            {
                for (int iunocc = 0; iunocc != n_states; iunocc++)
                {
                    if (iocc < iunocc)
                    {
                        double egap =
                            (eigenvalues(ik, iocc) - eigenvalues(ik, iunocc));  // * HA2EV;
                        double factor;
                        if (use_soc)
                            factor = wg.c[iocc] - wg.c[iunocc];
                        else
                            factor = (wg.c[iocc] - wg.c[iunocc]) / 2 * n_spin;
                        if (factor > 1.e-8)
                        {
                            for (int alpha = 0; alpha != 3; alpha++)
                            {
                                for (int beta = 0; beta != 3; beta++)
                                {
                                    for (int iomega = 0; iomega != this->omega.size(); iomega++)
                                    {
                                        double omega_ev = this->omega[iomega];  // * HA2EV;
                                        tmp = 2.0 * factor * velocity[ik][alpha](iunocc, iocc) *
                                              velocity[ik][beta](iocc, iunocc) /
                                              (egap * egap + omega_ev * omega_ev) / egap;
                                        this->head.at(iomega)(alpha, beta) -= tmp;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    for (int alpha = 0; alpha != 3; alpha++)
    {
        for (int beta = 0; beta != 3; beta++)
        {
            for (int iomega = 0; iomega != this->omega.size(); iomega++)
            {
                if (use_soc)
                    this->head.at(iomega)(alpha, beta) *= dielectric_unit;
                else
                    this->head.at(iomega)(alpha, beta) *= dielectric_unit * 2.0 / n_spin;
                if (alpha == beta)
                {
                    this->head.at(iomega)(alpha, beta) += std::complex<double>(1.0, 0.0);
                }
            }
        }
    }
    global::ofs_myid << "* Success: calculate head term." << std::endl;
    profiler.start("cal_head");
};

double diele_func::cal_factor(std::string name)
{
    using librpa_int::TWO_PI;
    using librpa_int::BOHR2ANG;

    double dielectric_unit;
    //! Bohr to A
    const double primitive_cell_volume = std::abs(pbc_.latvec.Det());
    // latvec.print();
    if (name == "head")
    {
        // abacus
        /*dielectric_unit = TWO_PI * hbar / h_divide_e2 / primitive_cell_volume /
                          this->meanfield_df.get_n_kpoints() * 1.0e30 / epsilon0 / eV;*/
        // aims
        dielectric_unit = 2 * TWO_PI / primitive_cell_volume;
    }
    else if (name == "wing")
    {
        // abacus
        /* dielectric_unit = TWO_PI * hbar / h_divide_e2 * sqrt(2 * TWO_PI / primitive_cell_volume)
         * 1.0e15 / this->meanfield_df.get_n_kpoints() / epsilon0 / TWO_PI / eV; */
        // aims
        dielectric_unit = 2 * sqrt(2 * TWO_PI / primitive_cell_volume);  // bohr
    }
    else
        throw std::logic_error("Unsupported value for head/wing factor");
    return dielectric_unit;
};

void diele_func::set_0_wing()
{
    int n_lambda = this->n_nonsingular - 1;
    for (int alpha = 0; alpha != 3; alpha++)
    {
        for (int iomega = 0; iomega != this->omega.size(); iomega++)
        {
            for (int mu = 0; mu != n_abf; mu++)
            {
                this->wing_mu.at(iomega)(mu, alpha) = 0.0;
            }
        }
    }
};

void diele_func::cal_wing(const librpa_int::Cs_LRI &Cs_data)
{
    using global::profiler;

    const bool use_soc = false;

    profiler.start("cal_wing_mu");
    init_wing();
    int n_lambda = this->n_nonsingular - 1;
    std::vector<std::complex<double>> local_wing_mu;
    local_wing_mu.resize(this->omega.size() * 3 * n_abf, 0.0);

    // #pragma omp parallel for schedule(dynamic) collapse(2)
    for (int mu = 0; mu < n_abf; ++mu)
    {
        for (int ik = 0; ik != nk; ik++)
        {
            // mpi_comm_global_h.barrier();
            // profiler.start("transform_Cs2mnk");
            auto desc_C_mnk = transform_Cs2mnk(ik, mu);
            // profiler.stop("transform_Cs2mnk");
            auto &desc_nband_nband = desc_C_mnk.first;
            auto &C_mnk = desc_C_mnk.second;
            // if (mu == 0 && mpi_comm_global_h.is_root())
            // {
            //     auto &velocity = this->meanfield_df.get_velocity();
            //     auto i_3 = desc_nband_nband.indx_g2l_r(3);
            //     auto j_4 = desc_nband_nband.indx_g2l_c(4);

            //     std::cout << "C,p: " << C_mnk(i_3, j_4) << "," << velocity[isp][ik][0](4, 3)
            //               << std::endl;
            // }
            // profiler.start("compute_wing");
            for (int iomega = 0; iomega != this->omega.size(); iomega++)
            {
                for (int alpha = 0; alpha != 3; alpha++)
                {
                    for (int isp = 0; isp != n_spin; isp++)
                    {
                        std::complex<double> tmp =
                            compute_wing(alpha, iomega, mu, ik, isp, desc_nband_nband, C_mnk);
                        const int index = mu * this->omega.size() * 3 + iomega * 3 + alpha;
                        local_wing_mu[index] += tmp;
                    }
                }
            }
            // profiler.stop("compute_wing");
        }
    }
    profiler.start("Comm_wing");
    MPI_Allreduce(MPI_IN_PLACE, local_wing_mu.data(), static_cast<int>(local_wing_mu.size()),
                  MPI_CXX_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
    profiler.stop("Comm_wing");
    double dielectric_unit = cal_factor("wing");

    for (int alpha = 0; alpha != 3; alpha++)
    {
        for (int mu = 0; mu != n_abf; mu++)
        {
            for (int iomega = 0; iomega != this->omega.size(); iomega++)
            {
                const int index = mu * this->omega.size() * 3 + iomega * 3 + alpha;
                this->wing_mu.at(iomega)(mu, alpha) = local_wing_mu[index];
                if (use_soc)
                    this->wing_mu.at(iomega)(mu, alpha) *= -dielectric_unit;
                else
                    this->wing_mu.at(iomega)(mu, alpha) *= -dielectric_unit * 2.0 / n_spin;
            }
        }
    }
    // transform_mu_to_lambda();
    if (mpi_comm_global_h.is_root()) std::cout << "* Success: calculate wing term." << std::endl;

    // this->wing_mu.clear();
    this->Coul_vector.clear();
    this->Coul_value.clear();
    this->meanfield_df.get_velocity().clear();
    release_free_mem();
    profiler.stop("cal_wing_mu");
};

std::pair<ArrayDesc, matrix_m<complex<double>>> diele_func::transform_Cs2mnk(const int ik,
                                                                             const int mu)
{
    using global::profiler;

    const bool use_shrink_abfs = false;
    const int n_soc = meanfield_df.get_n_soc();
    const int Mu = atomic_basis_abf_.get_i_atom(mu);
    const int mu_local = atomic_basis_abf_.get_local_index(mu, Mu);
    const int n_ao_Mu = atomic_basis_wfc_.get_atom_nb(Mu);

    ArrayDesc desc_nband_nao(blacs_ctxt_global_h);
    desc_nband_nao.init_1b1p(n_states, n_basis, 0, 0);
    ArrayDesc desc_nband_Mu(blacs_ctxt_global_h);
    desc_nband_Mu.init_1b1p(n_states, n_ao_Mu, 0, 0);
    ArrayDesc desc_Mu_nao(blacs_ctxt_global_h);
    desc_Mu_nao.init_1b1p(n_ao_Mu, n_basis, 0, 0);
    ArrayDesc desc_nao_nao(blacs_ctxt_global_h);
    desc_nao_nao.init_1b1p(n_basis, n_basis, 0, 0);
    ArrayDesc desc_nband_nband(blacs_ctxt_global_h);
    desc_nband_nband.init_1b1p(n_states, n_states, 0, 0);

    auto C_nao_nao = init_local_mat<complex<double>>(desc_Mu_nao, MAJOR::COL);
    auto C_nband_nband = init_local_mat<complex<double>>(desc_nband_nband, MAJOR::COL);

    const auto kfrac = kfrac_band[ik];
    const std::function<complex<double>(const int &, const std::pair<int, std::array<int, 3>> &)>
        fourier = [kfrac](const int &I, const std::pair<int, std::array<int, 3>> &J_Ra)
    {
        const auto &Ra = J_Ra.second;
        Vector3<double> R_IJ(Ra[0], Ra[1], Ra[2]);
        const auto ang = (kfrac * R_IJ) * TWO_PI;
        return complex<double>{std::cos(ang), std::sin(ang)};
    };
    // profiler.start("fourier");
    const auto set_IJ_nao_nao = get_necessary_IJ_from_block_2D(
        atomic_basis_wfc_, atomic_basis_wfc_, desc_nao_nao);
    auto s0_s1 = get_s0_s1_for_comm_map2_first(set_IJ_nao_nao);
    s0_s1.first.clear();
    s0_s1.first.insert(0);

    std::map<libri_types<int, int>::TAC, Tensor<double>> data_libri;
    std::map<int, std::map<std::pair<int, std::array<int, 3>>, Tensor<double>>> Cs_I_JR_local;
    // profiler.start("fourier_assign");
    if (use_shrink_abfs)
    {
        data_libri = Cs_shrinked_data.data_libri.at(Mu);
    }
    else
    {
        data_libri = Cs_data.data_libri.at(Mu);
    }
    const auto &n_Mu = atomic_basis_wfc_.get_atom_nb(Mu);
    for (const auto &JR_C : data_libri)
    {
        const auto J = JR_C.first.first;
        const auto R = JR_C.first.second;
        const auto &n_J = atomic_basis_wfc_.get_atom_nb(J);

        const std::initializer_list<std::size_t> shape{static_cast<std::size_t>(n_Mu),
                                                       static_cast<std::size_t>(n_J)};

        Cs_I_JR_local[0][{J, R}] = Tensor<double>(shape);

#pragma omp parallel for schedule(dynamic) collapse(2)
        for (int i = 0; i < n_Mu; i++)
        {
            for (int j = 0; j < n_J; j++)
            {
                Cs_I_JR_local[0][{J, R}](i, j) = JR_C.second(mu_local, i, j);
            }
        }
    }
    // profiler.stop("fourier_assign");
    // profiler.start("fourier_comm");
    auto Cs_I_JR =
        comm_map2_first(mpi_comm_global_h.comm, Cs_I_JR_local, s0_s1.first, s0_s1.second);
    // profiler.stop("fourier_comm");
    Cs_I_JR_local.clear();
    data_libri.clear();
    C_nao_nao.zero_out();
    C_nband_nband.zero_out();
    std::map<atom_t, size_t> atom_Mu;
    AtomicBasis atomic_basis_wfc_row;
    atom_Mu[0] = atomic_basis_wfc_.get_atom_nb(Mu);
    atomic_basis_wfc_row.set(atom_Mu);
    // profiler.start("fourier_collect");
    collect_block_from_IJ_storage_tensor_transform_triple(
        C_nao_nao, desc_Mu_nao, atomic_basis_wfc_row, atomic_basis_wfc_, fourier, Cs_I_JR, Mu);
    // profiler.stop("fourier_collect");
    Cs_I_JR.clear();
    // profiler.stop("fourier");
    // if (ik == 0 && mu == 1)
    // {
    //     for (int i = 0; i < n_ao_Mu; i++)
    //     {
    //         int J = 0;
    //         int j = 4;
    //         int ii = atomic_basis_wfc.get_global_index(0, i);
    //         int jj = atomic_basis_wfc.get_global_index(J, j);
    //         auto ii_loc = desc_Mu_nao.indx_g2l_r(ii);
    //         auto jj_loc = desc_Mu_nao.indx_g2l_c(jj);
    //         if (ii_loc >= 0 && jj_loc >= 0)
    //             ofs_myid << "Cij: " << i << ", " << C_nao_nao(ii_loc, jj_loc) << std::endl;
    //     }
    // }
    // profiler.start("scalapack_multiply");
    // prepare wave function BLACS
    for (int ispin = 0; ispin != n_spin; ispin++)
    {
        for (int is1 = 0; is1 != n_soc; is1++)
        {
            for (int is2 = 0; is2 != n_soc; is2++)
            {
                const auto &wfc_isp1_k = meanfield_df.get_eigenvectors()[ispin][is1][ik];
                ComplexMatrix wfc_Mu = ComplexMatrix(n_states, n_ao_Mu);
                // #pragma omp parallel for schedule collapse(2)
                for (int n = 0; n < n_states; n++)
                {
                    for (int i = 0; i < n_ao_Mu; i++)
                    {
                        const auto i_Mu = atomic_basis_wfc_.get_global_index(Mu, i);
                        wfc_Mu(n, i) = wfc_isp1_k(n, i_Mu);
                    }
                }
                const auto &wfc_isp2_k = meanfield_df.get_eigenvectors()[ispin][is2][ik];
                // blacs_ctxt_global_h.barrier();
                auto wfc1_block =
                    get_local_mat(conj(wfc_Mu).c, MAJOR::ROW, desc_nband_Mu, MAJOR::COL);
                auto wfc2_block =
                    get_local_mat(wfc_isp2_k.c, MAJOR::ROW, desc_nband_nao, MAJOR::COL);
                // if (ik == 0 && mu == 1)
                // {
                //     for (int i = 0; i < n_ao_Mu; i++)
                //     {
                //         auto ii_loc = desc_nband_Mu.indx_g2l_r(3);
                //         auto jj_loc = desc_nband_Mu.indx_g2l_c(i);
                //         if (ii_loc >= 0 && jj_loc >= 0)
                //             ofs_myid << "wfc1_block: " << i << ", " << wfc1_block(ii_loc, jj_loc)
                //                      << std::endl;
                //     }
                // }
                auto temp_nband_nao = multiply_scalapack(wfc1_block, desc_nband_Mu, C_nao_nao,
                                                         desc_Mu_nao, desc_nband_nao);
                // if (ik == 0 && mu == 1)
                // {
                //     auto i_3 = desc_nband_nao.indx_g2l_r(3);
                //     auto j_4 = desc_nband_nao.indx_g2l_c(4);
                //     if (i_3 >= 0 && j_4 >= 0)
                //         ofs_myid << "Cmn0: " << temp_nband_nao(i_3, j_4) << std::endl;
                // }
                ScalapackConnector::pgemm_f('N', 'T', n_states, n_states, n_basis, 1.0,
                                            temp_nband_nao.ptr(), 1, 1, desc_nband_nao.desc,
                                            wfc2_block.ptr(), 1, 1, desc_nband_nao.desc, 1.0,
                                            C_nband_nband.ptr(), 1, 1, desc_nband_nband.desc);
                auto C_tmp = init_local_mat<complex<double>>(desc_nband_nband, MAJOR::COL);
                ScalapackConnector::pgeadd_f('C', n_states, n_states, 1.0, C_nband_nband.ptr(), 1,
                                             1, desc_nband_nband.desc, 0.0, C_tmp.ptr(), 1, 1,
                                             desc_nband_nband.desc);
                ScalapackConnector::pgeadd_f('N', n_states, n_states, 1.0, C_tmp.ptr(), 1, 1,
                                             desc_nband_nband.desc, 1.0, C_nband_nband.ptr(), 1, 1,
                                             desc_nband_nband.desc);
                // if (ik == 0 && mu == 1)
                // {
                //     auto i_3 = desc_nband_nband.indx_g2l_r(3);
                //     auto j_4 = desc_nband_nband.indx_g2l_c(4);
                //     if (i_3 >= 0 && j_4 >= 0)
                //         ofs_myid << "Cmn1: " << C_nband_nband(i_3, j_4) << std::endl;
                // }
                // term2
                // wfc1_block = get_local_mat(wfc_Mu.c, MAJOR::ROW, desc_nband_Mu, MAJOR::COL);
                // wfc2_block =
                //     get_local_mat(conj(wfc_isp2_k).c, MAJOR::ROW, desc_nband_nao, MAJOR::COL);
                // auto tmp_band_Mu = init_local_mat<complex<double>>(desc_nband_Mu, MAJOR::COL);

                // ScalapackConnector::pgemm_f('N', 'C', n_states, n_ao_Mu, n_basis, 1.0,
                //                             wfc2_block.ptr(), 1, 1, desc_nband_nao.desc,
                //                             C_nao_nao.ptr(), 1, 1, desc_Mu_nao.desc, 0.0,
                //                             tmp_band_Mu.ptr(), 1, 1, desc_nband_Mu.desc);
                // ScalapackConnector::pgemm_f('N', 'T', n_states, n_states, n_ao_Mu, 1.0,
                //                             tmp_band_Mu.ptr(), 1, 1, desc_nband_Mu.desc,
                //                             wfc1_block.ptr(), 1, 1, desc_nband_Mu.desc, 1.0,
                //                             C_nband_nband.ptr(), 1, 1, desc_nband_nband.desc);
                // if (ik == 0 && mu == 1)
                // {
                //     auto i_3 = desc_nband_nband.indx_g2l_r(3);
                //     auto j_4 = desc_nband_nband.indx_g2l_c(4);
                //     if (i_3 >= 0 && j_4 >= 0)
                //         ofs_myid << "Cmn2: " << C_nband_nband(i_3, j_4) << std::endl;
                // }
            }
        }
    }
    // profiler.stop("scalapack_multiply");
    return std::make_pair(desc_nband_nband, C_nband_nband);
};

std::complex<double> diele_func::compute_wing(const int alpha, const int iomega, const int mu,
                                              const int ik, const int ispin,
                                              const ArrayDesc &desc_nband_nband,
                                              const matrix_m<complex<double>> &C_nband_nband)
{
    auto &velocity = this->meanfield_df.get_velocity();
    auto &eigenvalues = this->meanfield_df.get_eigenvals();

    double omega_ev = this->omega[iomega];  // * HA2EV;
    std::complex<double> wing_term = 0.0;

    auto &wg = this->meanfield_df.get_weight()[ispin];
    std::complex<double> test_tot = 0.0;
    for (int iocc = 0; iocc != n_states; iocc++)
    {
        for (int iunocc = 0; iunocc != n_states; iunocc++)
        {
            double egap =
                (eigenvalues[ispin](ik, iunocc) - eigenvalues[ispin](ik, iocc));  // * HA2EV;
            if (iocc < iunocc)
            {
                double factor1;
                double factor2;
                if (Params::use_soc)
                {
                    factor1 = wg.c[iocc] * (1.0 - wg.c[iunocc] * nk);
                    factor2 = wg.c[iunocc] * (1.0 - wg.c[iocc] * nk);
                }
                else
                {
                    factor1 = wg.c[iocc] / 2 * n_spin * (1.0 - wg.c[iunocc] / 2 * n_spin * nk);
                    factor2 = wg.c[iunocc] / 2 * n_spin * (1.0 - wg.c[iocc] / 2 * n_spin * nk);
                }
                if (factor1 > 1.e-8)
                {
                    auto loc_m = desc_nband_nband.indx_g2l_r(iocc);
                    auto loc_n = desc_nband_nband.indx_g2l_c(iunocc);
                    if (loc_m < 0 || loc_n < 0) continue;
                    auto tmp = C_nband_nband(loc_m, loc_n);
                    wing_term += factor1 * conj(tmp * velocity[ispin][ik][alpha](iunocc, iocc)) /
                                 (omega_ev * omega_ev + egap * egap);
                }
                if (factor2 > 1.e-8)
                {
                    auto loc_m = desc_nband_nband.indx_g2l_r(iocc);
                    auto loc_n = desc_nband_nband.indx_g2l_c(iunocc);
                    if (loc_m < 0 || loc_n < 0) continue;
                    auto tmp = C_nband_nband(loc_m, loc_n);
                    // for metal
                    wing_term += factor2 * tmp * velocity[ispin][ik][alpha](iunocc, iocc) /
                                 (omega_ev * omega_ev + egap * egap);
                }

                /*if (Params::debug)
                {
                    if (alpha == 0 && iomega == 0 && mu == 0)
                    {
                        std::complex<double> test =
                conj(Ctri_mn[mu][iocc][iunocc][kfrac_band[ik]] *
                                                         velocity[ispin][ik][alpha](iunocc,
                iocc)) / (omega_ev * omega_ev + egap * egap); if (iocc == 0 && iunocc == 10)
                        {
                            std::cout << "C,p: " << Ctri_mn[mu][iocc][iunocc][kfrac_band[ik]] <<
                ","
                                      << velocity[ispin][ik][alpha](iunocc, iocc) << std::endl;
                        }
                        test_tot += test;
                    }
                }*/
            }
        }
    }

    return wing_term;
};

void diele_func::wing_mu_to_lambda(matrix_m<std::complex<double>> &sqrtveig_blacs,
                                   ArrayDesc &desc_nabf_nabf_opt)
{
    using global::profiler;

    profiler.start("cal_wing");
    int n_lambda = this->n_nonsingular - 1;
    ArrayDesc desc_wing_mu(blacs_ctxt_global_h);
    desc_wing_mu.init_square_blk(n_abf, 3, 0, 0);
    ArrayDesc desc_wing(blacs_ctxt_global_h);
    desc_wing.init_square_blk(n_nonsingular - 1, 3, 0, 0);
    for (int iomega = 0; iomega != this->omega.size(); iomega++)
    {
        auto &wing_tmp = this->wing.at(iomega);
        wing_tmp = init_local_mat<complex<double>>(desc_wing, MAJOR::COL);
        // TODO: reconstruct wing_mu
        auto wing_mu_tmp = init_local_mat<complex<double>>(desc_wing_mu, MAJOR::COL);
        for (int alpha = 0; alpha != 3; alpha++)
        {
            auto loc_alpha = desc_wing_mu.indx_g2l_c(alpha);
            if (loc_alpha < 0) continue;
            for (int mu = 0; mu != n_abf; mu++)
            {
                auto loc_mu = desc_wing_mu.indx_g2l_r(mu);
                if (loc_mu < 0) continue;
                wing_mu_tmp(loc_mu, loc_alpha) = this->wing_mu.at(iomega)(mu, alpha);
            }
        }
        // this->wing.at(iomega)(lambda, alpha) +=
        // conj(sqrtveig_blacs(mu, lambda)) * this->wing_mu.at(iomega)(mu, alpha);
        // drop the first column of sqrtveig_blacs, the largest eigenvalue
        ScalapackConnector::pgemm_f('C', 'N', n_lambda, 3, n_abf, 1.0, sqrtveig_blacs.ptr(), 1, 2,
                                    desc_nabf_nabf_opt.desc, wing_mu_tmp.ptr(), 1, 1,
                                    desc_wing_mu.desc, 0.0, wing_tmp.ptr(), 1, 1, desc_wing.desc);
    }

    this->wing_mu.clear();
    profiler.stop("cal_wing");
};

// // real double diagonalization
// // Note complex diagonalization conserves symmetry much better
// void diele_func::get_Xv_real(double vq_threshold, const librpa_int::atpair_k_cplx_mat_t &Vq)
// {
//     using namespace librpa_int;
//     using RI::Tensor;
//     using RI::Communicate_Tensors_Map_Judge::comm_map2_first;
//     using librpa_int::global::mpi_comm_global_h;
// 
//     this->Coul_vector.clear();
//     this->Coul_value.clear();
//     const double CONE = 1.0;
//     std::array<double, 3> qa = {0.0, 0.0, 0.0};
//     Vector3_Order<double> q = {0.0, 0.0, 0.0};
//     size_t n_singular;
//     vec<double> eigenvalues(n_abf);
// 
//     const auto &comm_h = mpi_comm_global_h;
// 
//     comm_h.barrier();
//     BlacsCtxtHandler blacs_h(comm_h.comm);
// 
//     ArrayDesc desc_nabf_nabf(blacs_h);
//     desc_nabf_nabf.init_square_blk(n_abf, n_abf, 0, 0);
//     const auto set_IJ_nabf_nabf =
//         get_necessary_IJ_from_block_2D_sy('U', atomic_basis_abf_, desc_nabf_nabf);
//     const auto s0_s1 = get_s0_s1_for_comm_map2_first(set_IJ_nabf_nabf);
//     auto coul_eigen_block = init_local_mat<double>(desc_nabf_nabf, MAJOR::COL);
//     auto coulwc_block = init_local_mat<double>(desc_nabf_nabf, MAJOR::COL);
//     coulwc_block.zero_out();
//     std::map<int, std::map<std::pair<int, std::array<double, 3>>, RI::Tensor<double>>>
//         couleps_libri;
//     const int natom = atomic_basis_abf_.n_atoms;
//     const auto atpair_local = dispatch_upper_triangular_tasks(
//         natom, blacs_h.myid, blacs_h.nprows, blacs_h.npcols,
//         blacs_h.myprow, blacs_h.mypcol);
//     for (const auto &Mu_Nu : atpair_local)
//     {
//         const auto Mu = Mu_Nu.first;
//         const auto Nu = Mu_Nu.second;
//         // ofs_myid << "Mu " << Mu << " Nu " << Nu << endl;
//         if (Vq.count(Mu) == 0 || Vq.at(Mu).count(Nu) == 0 ||
//             Vq.at(Mu).at(Nu).count(q) == 0)
//             continue;
//         auto Vq_cpl = *(Vq.at(Mu).at(Nu).at(q));
//         // if (Vq.count(Mu) == 0 || Vq.at(Mu).count(Nu) == 0 || Vq.at(Mu).at(Nu).count(q) == 0)
//         //     continue;
//         // auto Vq_cpl = *(Vq.at(Mu).at(Nu).at(q));
//         const auto &Vq0 = std::make_shared<matrix>(Vq_cpl.real());
//         // const auto &Vq0 = Vq.at(Mu).at(Nu).at(q);
//         const auto n_mu = atomic_basis_abf_.get_atom_nb(Mu);
//         const auto n_nu = atomic_basis_abf_.get_atom_nb(Nu);
//         std::valarray<double> Vq_va(Vq0->c, Vq0->size);
//         auto pvq = std::make_shared<std::valarray<double>>();
//         *pvq = Vq_va;
//         couleps_libri[Mu][{Nu, qa}] = RI::Tensor<double>({n_mu, n_nu}, pvq);
//     }
//     const auto IJq_coul = RI::Communicate_Tensors_Map_Judge::comm_map2_first(
//         mpi_comm_global_h.comm, couleps_libri, s0_s1.first, s0_s1.second);
//     collect_block_from_ALL_IJ_Tensor(coulwc_block, desc_nabf_nabf, atomic_basis_abf_, qa,
//                                      true, CONE, IJq_coul, MAJOR::ROW);
//     power_hemat_blacs_real(coulwc_block, desc_nabf_nabf, coul_eigen_block, desc_nabf_nabf,
//                            n_singular, eigenvalues.c, 1.0, vq_threshold);
//     this->n_nonsingular = n_abf - n_singular;
//     for (int iv = 1; iv != n_nonsingular; iv++)
//     {
//         // Here eigen solved by Scalapack is ascending order,
//         // however, what we want is descending order.
//         this->Coul_value.push_back(eigenvalues.c[iv] + 0.0);  // throw away the largest one
//         std::vector<std::complex<double>> newRow;
// 
//         for (int jabf = 0; jabf != n_abf; jabf++)
//         {
//             newRow.push_back(coul_eigen_block(jabf, iv));
//         }
//         this->Coul_vector.push_back(newRow);
//     }
// 
//     if (mpi_comm_global_h.is_root())
//     {
//         std::cout << "The largest/smallest eigenvalue of Coulomb matrix(non-singular): "
//                   << this->Coul_value.front() << ", " << this->Coul_value.back() << std::endl;
//         std::cout << "The 1st/2nd/3rd/-1th eigenvalue of Coulomb matrix(Full): " << eigenvalues.c[0]
//                   << ", " << eigenvalues.c[1] << ", " << eigenvalues.c[2] << ", "
//                   << eigenvalues.c[n_abf - 1] << std::endl;
//         std::cout << "Dim of eigenvectors: " << coul_eigen_block.nr() << ", "
//                   << coul_eigen_block.nc() << std::endl;
//     }
//     /*std::cout << "Coulomb vector: lambda=-1" << std::endl;
//     for (int j = 0; j != n_abf; j++)
//     {
//         std::cout << j << "," << coul_eigen_block(j, 0) << std::endl;
//     }
//     std::cout << "Coulomb vector: lambda=0" << std::endl;
//     for (int j = 0; j != n_abf; j++)
//     {
//         std::cout << j << "," << Coul_vector[0][j] << std::endl;
//     }
//     std::cout << "Coulomb vector: lambda=1" << std::endl;
//     for (int j = 0; j != n_abf; j++)
//     {
//         std::cout << j << "," << Coul_vector[1][j] << std::endl;
//     }*/
//     if (mpi_comm_global_h.is_root())
//         std::cout << "* Success: diagonalize Coulomb matrix in the ABFs repre." << std::endl;
// };

void diele_func::get_Xv_cpl(double vq_threshold, const librpa_int::atpair_k_cplx_mat_t &Vq)
{
    using global::profiler;

    profiler.start("get_eigenvector_of_Coulomb_matrix");
    this->Coul_vector.clear();
    this->Coul_value.clear();
    const complex<double> CONE{1.0, 0.0};
    std::array<double, 3> qa = {0.0, 0.0, 0.0};
    Vector3_Order<double> q = {0.0, 0.0, 0.0};
    size_t n_singular;
    librpa_int::vec<double> eigenvalues(n_abf);

    const MpiCommHandler &comm_h = global::mpi_comm_global_h;  // TODO: replace with the actual comm

    comm_h.barrier();

    BlacsCtxtHandler blacs_ctxt_global_h;

    ArrayDesc desc_nabf_nabf(blacs_ctxt_global_h);
    desc_nabf_nabf.init_square_blk(n_abf, n_abf, 0, 0);
    const auto set_IJ_nabf_nabf =
        get_necessary_IJ_from_block_2D_sy('U', atomic_basis_abf_, desc_nabf_nabf);
    const auto s0_s1 = get_s0_s1_for_comm_map2_first(set_IJ_nabf_nabf);
    auto coul_eigen_block = init_local_mat<complex<double>>(desc_nabf_nabf, MAJOR::COL);
    auto coulwc_block = init_local_mat<complex<double>>(desc_nabf_nabf, MAJOR::COL);
    coulwc_block.zero_out();
    std::map<int, std::map<std::pair<int, std::array<double, 3>>, RI::Tensor<complex<double>>>>
        couleps_libri;
    const int natom = atomic_basis_abf_.n_atoms;
    const auto atpair_local = dispatch_upper_triangular_tasks(
        natom, blacs_ctxt_global_h.myid, blacs_ctxt_global_h.nprows, blacs_ctxt_global_h.npcols,
        blacs_ctxt_global_h.myprow, blacs_ctxt_global_h.mypcol);
    for (const auto &Mu_Nu : atpair_local)
    {
        const auto Mu = Mu_Nu.first;
        const auto Nu = Mu_Nu.second;
        // ofs_myid << "Mu " << Mu << " Nu " << Nu << endl;
        if (Vq.count(Mu) == 0 || Vq.at(Mu).count(Nu) == 0 ||
            Vq.at(Mu).at(Nu).count(q) == 0)
            continue;
        const auto &Vq0 = Vq.at(Mu).at(Nu).at(q);
        const auto n_mu = atomic_basis_abf_.get_atom_nb(Mu);
        const auto n_nu = atomic_basis_abf_.get_atom_nb(Nu);
        std::valarray<complex<double>> Vq_va(Vq0->c, Vq0->size);
        auto pvq = std::make_shared<std::valarray<complex<double>>>();
        *pvq = Vq_va;
        couleps_libri[Mu][{Nu, qa}] = RI::Tensor<complex<double>>({n_mu, n_nu}, pvq);
    }
    const auto IJq_coul = RI::Communicate_Tensors_Map_Judge::comm_map2_first(
        comm_h.comm, couleps_libri, s0_s1.first, s0_s1.second);
    collect_block_from_ALL_IJ_Tensor(coulwc_block, desc_nabf_nabf, atomic_basis_abf_, qa,
                                     true, CONE, IJq_coul, MAJOR::ROW);
    power_hemat_blacs_desc(coulwc_block, desc_nabf_nabf, coul_eigen_block, desc_nabf_nabf,
                           n_singular, eigenvalues.c, 1.0, vq_threshold);
    this->n_nonsingular = n_abf - n_singular;

    for (int iv = 1; iv != n_nonsingular; iv++)
    {
        // Here eigen solved by Scalapack is ascending order,
        // however, what we want is descending order.
        this->Coul_value.push_back(eigenvalues.c[iv]);  // throw away the largest one
        std::vector<std::complex<double>> newRow;

        for (int jabf = 0; jabf != n_abf; jabf++)
        {
            newRow.push_back(coul_eigen_block(jabf, iv));
        }
        this->Coul_vector.push_back(newRow);
        newRow.clear();
    }
    if (mpi_comm_global_h.is_root())
    {
        std::cout << "n_singular: " << n_singular << std::endl;
        std::cout << "The largest/smallest eigenvalue of Coulomb matrix(non-singular): "
                  << this->Coul_value.front() << ", " << this->Coul_value.back() << std::endl;
        std::cout << "The 1st/2nd/3rd/-1th eigenvalue of Coulomb matrix(Full): " << eigenvalues.c[0]
                  << ", " << eigenvalues.c[1] << ", " << eigenvalues.c[2] << ", "
                  << eigenvalues.c[n_abf - 1] << std::endl;
        // std::cout << "Dim of eigenvectors: " << coul_eigen_block.dataobj.nr() << ", "
        //           << coul_eigen_block.dataobj.nc() << std::endl;

        std::cout << "* Success: diagonalize Coulomb matrix in the ABFs repre.\n";
    }
    profiler.stop("get_eigenvector_of_Coulomb_matrix");
};

std::vector<double> diele_func::get_head_vec()
{
    std::vector<double> head_vec(this->omega.size(), 0.0);
    for (int iomega = 0; iomega != this->omega.size(); iomega++)
    {
        std::complex<double> df = 0;
        for (int alpha = 0; alpha != 3; alpha++)
        {
            df += this->head.at(iomega)(alpha, alpha);
        }
        head_vec[iomega] = df.real() / 3.0;
    }
    return head_vec;
};

void diele_func::test_head()
{
    using LIBRPA::utils::lib_printf;
    if (mpi_comm_global_h.is_root())
    {
        std::cout << "Freqency node & Head of dielectric function(Re, Im): " << std::endl;
        for (int iomega = 0; iomega != this->omega.size(); iomega++)
        {
            std::complex<double> df = 0;
            for (int alpha = 0; alpha != 3; alpha++)
            {
                df += this->head.at(iomega)(alpha, alpha);
            }
            lib_printf("%2d %15.8f %15.8f %15.8f\n", iomega, this->omega[iomega], df.real() / 3.0,
                       df.imag() / 3.0);
        }
    }
    // std::exit(0);
};

void diele_func::test_wing()
{
    using LIBRPA::utils::lib_printf;
    if (mpi_comm_global_h.is_root())
    {
        std::cout << "Index of abfs & wing(iomega=0, z) of dielectric function(Re, Im): "
                  << std::endl;
        for (int mu = 0; mu != n_abf; mu++)
        {
            std::complex<double> df = 0;
            // z direction
            df = this->wing_mu.at(0)(mu, 2);

            lib_printf("%2d %15.8f %15.8f\n", mu, df.real(), df.imag());
        }
    }
    // if (mpi_comm_global_h.myid == 1)
    // {
    //     std::cout << "id=1 Index of abfs & wing(iomega=0, z) of dielectric function(Re, Im): "
    //               << std::endl;
    //     for (int mu = 0; mu != n_abf; mu++)
    //     {
    //         std::complex<double> df = 0;
    //         // z direction
    //         df = this->wing_mu.at(0)(mu, 2);

    //         lib_printf("%2d %15.8f %15.8f\n", mu, df.real(), df.imag());
    //     }
    // }
    // std::exit(0);
};

ArrayDesc diele_func::get_body_inv(matrix_m<std::complex<double>> &chi0_block,
                                    ArrayDesc &desc_nabf_nabf_opt)
{
    using global::profiler;
    using global::mpi_comm_global_h;

    mpi_comm_global_h.barrier();
    BlacsCtxtHandler blacs_h(mpi_comm_global_h.comm);
    ArrayDesc desc_nabf_nabf(blacs_h);
    desc_nabf_nabf.init_square_blk(n_nonsingular - 1, n_nonsingular - 1, 0, 0);
    this->body_inv = init_local_mat<complex<double>>(desc_nabf_nabf, MAJOR::COL);
    profiler.start("get_inverse_body_of_chi0");
    mpi_comm_global_h.barrier();

    ArrayDesc desc_body(blacs_ctxt_global_h);
    desc_body.init_square_blk(n_nonsingular - 1, n_nonsingular - 1, 0, 0);
    this->body_inv = init_local_mat<complex<double>>(desc_body, MAJOR::COL);

    /* for (int ilambda = 0; ilambda < n_nonsingular - 1; ilambda++)
    {
        for (int jlambda = 0; jlambda < n_nonsingular - 1; jlambda++)
        {
            body_inv(ilambda, jlambda) = chi0_block(1 + ilambda, 1 + jlambda);
        }
    } */
    ScalapackConnector::pgemr2d_f(n_nonsingular - 1, n_nonsingular - 1, chi0_block.ptr(), 2, 2,
                                  desc_nabf_nabf_opt.desc, this->body_inv.ptr(), 1, 1,
                                  desc_body.desc, blacs_ctxt_global_h.ictxt);
    invert_scalapack(this->body_inv, desc_body);

    // std::cout << "* Success: get inverse body of chi0.\n";
    profiler.stop("get_inverse_body_of_chi0");
    return desc_body;
};

void diele_func::construct_L(const int ifreq, ArrayDesc &desc_body)
{
    profiler.start("cal_L");
    this->Lind.resize(3, 3, MAJOR::COL);
    this->bw.resize(n_nonsingular - 1, 3, MAJOR::COL);
    this->wb.resize(3, n_nonsingular - 1, MAJOR::COL);
    ArrayDesc desc_wing(blacs_ctxt_global_h);
    desc_wing.init_square_blk(n_nonsingular - 1, 3, 0, 0);

    ArrayDesc desc_lam_3(blacs_ctxt_global_h);
    desc_lam_3.init_square_blk(n_nonsingular - 1, 3, 0, 0);

    ArrayDesc desc_3_lam(blacs_ctxt_global_h);
    desc_3_lam.init_square_blk(3, n_nonsingular - 1, 0, 0);

    ArrayDesc desc_3_3(blacs_ctxt_global_h);
    desc_3_3.init_square_blk(3, 3, 0, 0);

    auto lam_3 = init_local_mat<complex<double>>(desc_lam_3, MAJOR::COL);
    auto _3_lam = init_local_mat<complex<double>>(desc_3_lam, MAJOR::COL);
    auto Lind_loc = init_local_mat<complex<double>>(desc_3_3, MAJOR::COL);
    // tmp = head.at(ifreq) - transpose(wing.at(ifreq), true) * body_inv * wing.at(ifreq);
    ScalapackConnector::pgemm_f('N', 'N', n_nonsingular - 1, 3, n_nonsingular - 1, 1.0,
                                body_inv.ptr(), 1, 1, desc_body.desc, wing.at(ifreq).ptr(), 1, 1,
                                desc_wing.desc, 0.0, lam_3.ptr(), 1, 1, desc_lam_3.desc);
    ScalapackConnector::pgemm_f('C', 'N', 3, 3, n_nonsingular - 1, 1.0, wing.at(ifreq).ptr(), 1, 1,
                                desc_wing.desc, lam_3.ptr(), 1, 1, desc_lam_3.desc, 0.0,
                                Lind_loc.ptr(), 1, 1, desc_3_3.desc);
    ScalapackConnector::pgemm_f('C', 'N', 3, n_nonsingular - 1, n_nonsingular - 1, 1.0,
                                wing.at(ifreq).ptr(), 1, 1, desc_wing.desc, body_inv.ptr(), 1, 1,
                                desc_body.desc, 0.0, _3_lam.ptr(), 1, 1, desc_3_lam.desc);

    for (int i = 0; i != 3; i++)
    {
        auto loc_i = desc_3_3.indx_g2l_r(i);
        for (int ilambda = 0; ilambda < n_nonsingular - 1; ilambda++)
        {
            auto loc_ilambda = desc_lam_3.indx_g2l_r(ilambda);
            auto loc_ibw = desc_lam_3.indx_g2l_c(i);
            if (loc_ibw >= 0 && loc_ilambda >= 0)
                this->bw(ilambda, i) = lam_3(loc_ilambda, loc_ibw);

            loc_ilambda = desc_3_lam.indx_g2l_c(ilambda);
            auto loc_iwb = desc_3_lam.indx_g2l_r(i);
            if (loc_iwb >= 0 && loc_ilambda >= 0)
                this->wb(i, ilambda) = _3_lam(loc_iwb, loc_ilambda);

            MPI_Allreduce(MPI_IN_PLACE, &bw(ilambda, i), 1, MPI_DOUBLE_COMPLEX, MPI_SUM,
                          mpi_comm_global_h.comm);
            MPI_Allreduce(MPI_IN_PLACE, &wb(i, ilambda), 1, MPI_DOUBLE_COMPLEX, MPI_SUM,
                          mpi_comm_global_h.comm);
        }

        for (int j = 0; j != 3; j++)
        {
            auto loc_j = desc_3_3.indx_g2l_c(j);
            if (loc_j >= 0 && loc_i >= 0)
                this->Lind(i, j) = head.at(ifreq)(i, j) - Lind_loc(loc_i, loc_j);
            MPI_Allreduce(MPI_IN_PLACE, &Lind(i, j), 1, MPI_DOUBLE_COMPLEX, MPI_SUM,
                          mpi_comm_global_h.comm);
        }
    }

    profiler.stop("cal_L");
};

void diele_func::get_Leb_points()
{
    auto quad_order = lebedev::QuadratureOrder::order_5810;
    // lebedev::QuadratureOrder::order_590;
    auto quad_points = lebedev::QuadraturePoints(quad_order);
    qx_leb = quad_points.get_x();
    qy_leb = quad_points.get_y();
    qz_leb = quad_points.get_z();
    qw_leb = quad_points.get_weights();
    for (int ileb = 0; ileb != qw_leb.size(); ileb++)
    {
        qw_leb[ileb] *= 2 * TWO_PI;
    }
};

void diele_func::get_g_enclosing_gamma()
{
    g_enclosing_gamma.clear();
    g_enclosing_gamma.resize(26);
    const auto &kv_nmp = pbc_.period_array;
    const auto &G = pbc_.G;
    int ik = 0;
    for (int a = -1; a != 2; a++)
    {
        for (int b = -1; b != 2; b++)
        {
            for (int c = -1; c != 2; c++)
            {
                if (a == 0 && b == 0 && c == 0) continue;
                g_enclosing_gamma.at(ik) = {
                    G.e11 * a / kv_nmp[0] + G.e12 * b / kv_nmp[1] + G.e13 * c / kv_nmp[2],
                    G.e21 * a / kv_nmp[0] + G.e22 * b / kv_nmp[1] + G.e23 * c / kv_nmp[2],
                    G.e31 * a / kv_nmp[0] + G.e32 * b / kv_nmp[1] + G.e33 * c / kv_nmp[2]};

                ik++;
            }
        }
    }
};

void diele_func::calculate_q_gamma()
{
    q_gamma.clear();
    q_gamma.resize(qw_leb.size());
#pragma omp parallel for schedule(dynamic)
    for (int ileb = 0; ileb != qw_leb.size(); ileb++)
    {
        double qmax = 1.0e10;
        Vector3_Order<double> q_quta = {qx_leb[ileb], qy_leb[ileb], qz_leb[ileb]};
        for (int ik = 0; ik != 26; ik++)
        {
            double denominator = q_quta * g_enclosing_gamma[ik];
            if (denominator > 1.0e-10)
            {
                double numerator = 0.5 * g_enclosing_gamma[ik] * g_enclosing_gamma[ik];
                double temp = numerator / denominator;
                qmax = std::min(qmax, temp);
            }
        }
        q_gamma[ileb] = qmax;
    }
};

void diele_func::cal_eps(const int ifreq, ArrayDesc &desc_nabf_nabf_opt, ArrayDesc &desc_body)
{
    using global::mpi_comm_global_h;
    using global::profiler;

    profiler.start("cal_inverse_dielectric_matrix");
    this->chi0 = init_local_mat<complex<double>>(desc_nabf_nabf_opt, MAJOR::COL);

    const double k_volume = std::abs(pbc_.G.Det());
    this->vol_gamma = k_volume / nk;
    double vol_gamma_numeric = 0.0;
    if (ifreq == 0 && mpi_comm_global_h.is_root())
    {
        for (int ileb = 0; ileb != qw_leb.size(); ileb++)
        {
            vol_gamma_numeric += qw_leb[ileb] * std::pow(q_gamma[ileb], 3) / 3.0;
        }
        std::cout << "Number of angular grids for average inverse dielectric matrix: "
                  << qw_leb.size() << std::endl;
        std::cout << "vol_gamma_numeric/vol_gamma: " << vol_gamma_numeric << ", " << vol_gamma
                  << std::endl;
        std::cout << "Angular quadrature accuracy for volume: " << vol_gamma_numeric / vol_gamma
                  << " (should be close to 1)" << std::endl;
    }
    /*std::cout << "major of Matz: " << wing[0].is_row_major() << "," << body_inv.is_row_major()
              << "," << transpose(wing.at(0), true).is_row_major() << "," << Lind.is_row_major()
              << std::endl;*/
    construct_L(ifreq, desc_body);

    profiler.start("precompute_q_data");

    const size_t nleb = qw_leb.size();
    std::vector<std::complex<double>> weights(nleb);
    std::vector<std::array<double, 3>> q_vectors(nleb);

    const auto L00 = Lind(0, 0), L01 = Lind(0, 1), L02 = Lind(0, 2);
    const auto L10 = Lind(1, 0), L11 = Lind(1, 1), L12 = Lind(1, 2);
    const auto L20 = Lind(2, 0), L21 = Lind(2, 1), L22 = Lind(2, 2);

#pragma omp parallel for schedule(static)
    for (int ileb = 0; ileb < nleb; ++ileb)
    {
        const double qx = qx_leb[ileb];
        const double qy = qy_leb[ileb];
        const double qz = qz_leb[ileb];

        q_vectors[ileb] = {qx, qy, qz};

        const auto qLq = qx * (qx * L00 + qy * L01 + qz * L02) +
                         qy * (qx * L10 + qy * L11 + qz * L12) +
                         qz * (qx * L20 + qy * L21 + qz * L22);

        weights[ileb] = qw_leb[ileb] * std::pow(q_gamma[ileb], 3) / (3.0 * vol_gamma) / qLq;
    }
    profiler.stop("precompute_q_data");

    profiler.start("cal_inverse_dielectric_matrix_ij");
    int i_start = 0, i_end = n_nonsingular;
    int j_start = 0, j_end = n_nonsingular;
#pragma omp parallel for schedule(dynamic, 4) collapse(2)
    for (int i = i_start; i != i_end; i++)
    {
        for (int j = j_start; j != j_end; j++)
        {
            const int ilo = desc_nabf_nabf_opt.indx_g2l_r(i);
            if (ilo < 0) continue;
            const int jlo = desc_nabf_nabf_opt.indx_g2l_c(j);
            if (jlo < 0) continue;

            complex<double> result = 0.0;

            if (i == 0 && j == 0)
            {
                for (int ileb = 0; ileb < nleb; ++ileb)
                {
                    result += weights[ileb];
                }
            }
            else if (i == 0 || j == 0)
            {
                result = 0.0;
            }
            else
            {
                const int idx_i = i - 1, idx_j = j - 1;

                const auto bw_i0 = bw(idx_i, 0), bw_i1 = bw(idx_i, 1), bw_i2 = bw(idx_i, 2);
                const auto wb_j0 = wb(0, idx_j), wb_j1 = wb(1, idx_j), wb_j2 = wb(2, idx_j);

                for (int ileb = 0; ileb < nleb; ++ileb)
                {
                    const auto &[qx, qy, qz] = q_vectors[ileb];
                    const auto bwq = bw_i0 * qx + bw_i1 * qy + bw_i2 * qz;
                    const auto qwb = qx * wb_j0 + qy * wb_j1 + qz * wb_j2;

                    result += weights[ileb] * bwq * qwb;
                }
            }
            chi0(ilo, jlo) = result;
        }
    }
    auto identity = init_local_mat<complex<double>>(desc_body, MAJOR::COL);
    for (int i = 0; i < n_nonsingular - 1; i++)
    {
        const int ilo = desc_body.indx_g2l_r(i);
        if (ilo < 0) continue;
        for (int j = 0; j < n_nonsingular - 1; j++)
        {
            const int jlo = desc_body.indx_g2l_c(j);
            if (jlo < 0) continue;
            if (i == j)
                identity(ilo, jlo) = 1.0;
            else
                identity(ilo, jlo) = 0.0;
        }
    }
    ScalapackConnector::pgemm_f('N', 'N', n_nonsingular - 1, n_nonsingular - 1, n_nonsingular - 1,
                                1.0, body_inv.ptr(), 1, 1, desc_body.desc, identity.ptr(), 1, 1,
                                desc_body.desc, 1.0, chi0.ptr(), 2, 2, desc_nabf_nabf_opt.desc);
    profiler.stop("cal_inverse_dielectric_matrix_ij");
    if (mpi_comm_global_h.is_root())
        std::cout << "* Success: calculate average inverse dielectric matrix no." << ifreq + 1
                  << "." << std::endl;
    profiler.stop("cal_inverse_dielectric_matrix");
};

/*std::complex<double> diele_func::compute_chi0_inv_00(const int ifreq)
{
    std::complex<double> total = 0.0;
    std::vector<std::complex<double>> partial_sum(qw_leb.size(), 0.0);
#pragma omp parallel for schedule(dynamic)
    for (int ileb = 0; ileb != qw_leb.size(); ileb++)
    {
        matrix_m<std::complex<double>> q_unit(3, 1, MAJOR::COL);
        q_unit(0, 0) = qx_leb[ileb];
        q_unit(1, 0) = qy_leb[ileb];
        q_unit(2, 0) = qz_leb[ileb];

        auto den = transpose(q_unit, false) * Lind * q_unit;
        // total += qw_leb[ileb] * std::pow(q_gamma[ileb], 3) / den(0, 0);
        partial_sum[ileb] = qw_leb[ileb] * std::pow(q_gamma[ileb], 3) / den(0, 0);
    }
    total = std::accumulate(partial_sum.begin(), partial_sum.end(), std::complex<double>(0.0,
0.0)); total *= 1.0 / 3.0 / vol_gamma;

    return total;
};

std::complex<double> diele_func::compute_chi0_inv_ij(const int ifreq, int i, int j)
{
    const std::complex<double> bw_i0 = this->bw(i, 0);
    const std::complex<double> bw_i1 = this->bw(i, 1);
    const std::complex<double> bw_i2 = this->bw(i, 2);
    const std::complex<double> wb_j0 = this->wb(0, j);
    const std::complex<double> wb_j1 = this->wb(1, j);
    const std::complex<double> wb_j2 = this->wb(2, j);

    const std::complex<double> L00 = Lind(0, 0);
    const std::complex<double> L01 = Lind(0, 1);
    const std::complex<double> L02 = Lind(0, 2);
    const std::complex<double> L10 = Lind(1, 0);
    const std::complex<double> L11 = Lind(1, 1);
    const std::complex<double> L12 = Lind(1, 2);
    const std::complex<double> L20 = Lind(2, 0);
    const std::complex<double> L21 = Lind(2, 1);
    const std::complex<double> L22 = Lind(2, 2);

    std::complex<double> total = 0.0;

    const size_t nleb = qw_leb.size();

#pragma omp parallel for reduction(+ : total)
    for (int ileb = 0; ileb < nleb; ++ileb)
    {
        const double qx = qx_leb[ileb];
        const double qy = qy_leb[ileb];
        const double qz = qz_leb[ileb];

        const std::complex<double> qLq = qx * (qx * L00 + qy * L01 + qz * L02) +
                                         qy * (qx * L10 + qy * L11 + qz * L12) +
                                         qz * (qx * L20 + qy * L21 + qz * L22);

        const std::complex<double> bwq = bw_i0 * qx + bw_i1 * qy + bw_i2 * qz;
        const std::complex<double> qwb = qx * wb_j0 + qy * wb_j1 + qz * wb_j2;

        total += qw_leb[ileb] * std::pow(q_gamma[ileb], 3) * bwq * qwb / qLq;
    }

    return total * (1.0 / (3.0 * vol_gamma));
}*/

void diele_func::assign_chi0(matrix_m<std::complex<double>> &chi0_block,
                             ArrayDesc &desc_nabf_nabf_opt)
{
    profiler.start("assign_chi0");
    mpi_comm_global_h.barrier();

    ScalapackConnector::pgemr2d_f(n_abf, n_abf, this->chi0.ptr(), 1, 1, desc_nabf_nabf_opt.desc,
                                  chi0_block.ptr(), 1, 1, desc_nabf_nabf_opt.desc,
                                  blacs_ctxt_global_h.ictxt);

    profiler.stop("assign_chi0");
}

void diele_func::rewrite_eps(matrix_m<std::complex<double>> &chi0_block, const int ifreq,
                             ArrayDesc &desc_nabf_nabf_opt)
{
    auto desc_body = get_body_inv(chi0_block, desc_nabf_nabf_opt);
    cal_eps(ifreq, desc_nabf_nabf_opt, desc_body);
    assign_chi0(chi0_block, desc_nabf_nabf_opt);
    // this->chi0.clear();
    this->Lind.clear();
    this->body_inv.clear();
};

}
