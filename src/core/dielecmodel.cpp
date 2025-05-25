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
#include "atomic_basis.h"
#include "pbc.h"
#include "ri.h"
// #include "utils_atomic_basis_blacs.h"
#include "../lebedev-quadrature/lebedev_quadrature.hpp"
#ifdef LIBRPA_USE_LIBRI
#include <RI/comm/mix/Communicate_Tensors_Map_Judge.h>
#include <RI/global/Tensor.h>
#else
#include "../utils/libri_stub.h"
#endif

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
    for (int iomega = 0; iomega != n_omega; iomega++)
    {
        wing_mu[iomega].resize(n_abf, 3, MAJOR::COL);
        wing[iomega].resize(n_nonsingular - 1, 3, MAJOR::COL);
    }
    get_Leb_points();
    get_g_enclosing_gamma();
    calculate_q_gamma();
    std::cout << "* Success: initalize and calculate lebdev points and q_gamma." << std::endl;
}

void diele_func::cal_head()
{
    using global::profiler;

    profiler.start("cal_head");
    //! spin = 1 only
    auto &wg = this->meanfield_df.get_weight()[n_spin - 1];
    auto &eigenvalues = this->meanfield_df.get_eigenvals()[n_spin - 1];
    auto &velocity = this->meanfield_df.get_velocity()[n_spin - 1];
    std::complex<double> tmp;
    int nocc = 0;

    double dielectric_unit = cal_factor("head");
    for (int i = 0; i != wg.size; i++)
    {
        if (wg.c[i] == 0.)
        {
            nocc = i;
            break;
        }
    }
    for (int ik = 0; ik != nk; ik++)
    {
        for (int iocc = 0; iocc != nocc; iocc++)
        {
            for (int iunocc = nocc; iunocc != n_states; iunocc++)
            {
                double egap = (eigenvalues(ik, iocc) - eigenvalues(ik, iunocc));  // * HA2EV;
                for (int alpha = 0; alpha != 3; alpha++)
                {
                    for (int beta = 0; beta != 3; beta++)
                    {
                        for (int iomega = 0; iomega != this->omega.size(); iomega++)
                        {
                            double omega_ev = this->omega[iomega];  // * HA2EV;
                            tmp = 2.0 * velocity[ik][alpha](iunocc, iocc) *
                                  velocity[ik][beta](iocc, iunocc) /
                                  (egap * egap + omega_ev * omega_ev) / egap;
                            this->head.at(iomega)(alpha, beta) -= tmp;
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
                if (n_spin == 1)
                {
                    if (Params::use_soc)
                        this->head.at(iomega)(alpha, beta) *= dielectric_unit;
                    else
                        this->head.at(iomega)(alpha, beta) *= dielectric_unit * 2;
                }
                else if (n_spin == 4)
                    this->head.at(iomega)(alpha, beta) *= dielectric_unit;
                if (alpha == beta)
                {
                    this->head.at(iomega)(alpha, beta) += std::complex<double>(1.0, 0.0);
                }
            }
        }
    }
    std::cout << "* Success: calculate head term." << std::endl;
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
        dielectric_unit = 2 * TWO_PI / primitive_cell_volume / this->meanfield_df.get_n_kpoints();
    }
    else if (name == "wing")
    {
        // abacus
        /* dielectric_unit = TWO_PI * hbar / h_divide_e2 * sqrt(2 * TWO_PI / primitive_cell_volume)
         * 1.0e15 / this->meanfield_df.get_n_kpoints() / epsilon0 / TWO_PI / eV; */
        // aims
        dielectric_unit = 2 * sqrt(2 * TWO_PI / primitive_cell_volume) /
                          this->meanfield_df.get_n_kpoints();  // bohr
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
            for (int lambda = 0; lambda != n_lambda; lambda++)
            {
                this->wing.at(iomega)(lambda, alpha) = 0.0;
            }
        }
    }
};

void diele_func::cal_wing(const librpa_int::Cs_LRI &Cs_data)
{
    profiler.start("cal_wing_mu");
    init_wing();
    int n_lambda = this->n_nonsingular - 1;
    init_Cs(Cs_data);
    FT_R2k(Cs_data);
    Cs_ij2mn();

#pragma omp parallel for schedule(dynamic) collapse(3)
    for (int iomega = 0; iomega != this->omega.size(); iomega++)
    {
        for (int alpha = 0; alpha != 3; alpha++)
        {
            for (int mu = 0; mu != n_abf; mu++)
            {
                this->wing_mu.at(iomega)(mu, alpha) = compute_wing(alpha, iomega, mu);
            }
        }
    }
    double dielectric_unit = cal_factor("wing");

    for (int alpha = 0; alpha != 3; alpha++)
    {
        for (int mu = 0; mu != n_abf; mu++)
        {
            for (int iomega = 0; iomega != this->omega.size(); iomega++)
            {
                if (n_spin == 1)
                {
                    if (Params::use_soc)
                        this->wing_mu.at(iomega)(mu, alpha) *= -dielectric_unit;
                    else
                        this->wing_mu.at(iomega)(mu, alpha) *= -dielectric_unit * 2.0;
                }
                else if (n_spin == 4)
                    this->wing_mu.at(iomega)(mu, alpha) *= -dielectric_unit;
            }
        }
    }
    // transform_mu_to_lambda();
    std::cout << "* Success: calculate wing term." << std::endl;

    // if (Params::debug)
    // {
    //     df_headwing.test_head();
    //     df_headwing.test_wing();
    // }
    // this->wing_mu.clear();
    this->Coul_vector.clear();
    this->Coul_value.clear();
    this->Ctri_mn.clear();
    this->Ctri_ij.clear();
    profiler.stop("cal_wing_mu");
};

void diele_func::wing_mu_to_lambda(matrix_m<std::complex<double>> &sqrtveig_blacs)
{
    int n_lambda = this->n_nonsingular - 1;
    for (int alpha = 0; alpha != 3; alpha++)
    {
        for (int iomega = 0; iomega != this->omega.size(); iomega++)
        {
            for (int lambda = 0; lambda != n_lambda; lambda++)
            {
                for (int mu = 0; mu != n_abf; mu++)
                {
                    this->wing.at(iomega)(lambda, alpha) +=
                        conj(sqrtveig_blacs(mu, lambda)) * this->wing_mu.at(iomega)(mu, alpha);
                }
            }
        }
    }
    this->wing_mu.clear();
};

std::complex<double> diele_func::compute_wing(int alpha, int iomega, int mu)
{
    const bool debug = false;
    auto &wg = this->meanfield_df.get_weight()[n_spin - 1];
    auto &velocity = this->meanfield_df.get_velocity();
    auto &eigenvalues = this->meanfield_df.get_eigenvals();
    int nocc = 0;
    for (int i = 0; i != wg.size; i++)
    {
        if (wg.c[i] == 0.)
        {
            nocc = i;
            break;
        }
    }
    double omega_ev = this->omega[iomega];  // * HA2EV;
    std::complex<double> wing_term = 0.0;

    // std::complex<double> tmp = 0.0;
    for (int ispin = 0; ispin != n_spin; ispin++)
    {
        for (int ik = 0; ik != nk; ik++)
        {
            std::complex<double> test_tot = 0.0;
            for (int iocc = 0; iocc != n_states; iocc++)
            {
                for (int iunocc = iocc; iunocc != n_states; iunocc++)
                {
                    double egap = (eigenvalues[ispin](ik, iunocc) -
                                   eigenvalues[ispin](ik, iocc));  // * HA2EV;
                    if (iocc < nocc && iunocc >= nocc)
                    {
                        /*tmp += conj(this->Ctri_mn[mu][iocc][iunocc][kfrac_band[ik]] *
                                    velocity[ispin][ik][alpha](iunocc, iocc)) /
                               (omega_ev * omega_ev + egap * egap);*/

                        wing_term += conj(this->Ctri_mn[mu][iocc][iunocc][kfrac_band[ik]] *
                                          velocity[ispin][ik][alpha](iunocc, iocc)) /
                                     (omega_ev * omega_ev + egap * egap);

                        /*if (iocc == 5 && iunocc == 40 && alpha == 0 && mu == 0 && iomega == 0)
                        {
                            std::cout << "mu, ik: " << mu << "," << ik << std::endl;
                            std::cout << "C: " << std::scientific << std::setprecision(8)
                                      << conj(this->Ctri_mn[mu][iocc][iunocc][kfrac_band[ik]])
                                      << std::endl;
                            std::cout << "p: " << std::scientific << std::setprecision(8)
                                      << conj(velocity[ispin][ik][alpha](iunocc, iocc))
                                      << std::endl;
                            std::cout << "E_m, E_n: " << std::scientific << std::setprecision(8)
                                      << eigenvalues[ispin](ik, iunocc) << "," << std::scientific
                                      << std::setprecision(8) << eigenvalues[ispin](ik, iocc)
                                      << std::endl;
                            std::cout << "C*p: "
                                      << conj(
                                              this->Ctri_mn[mu][iocc][iunocc][kfrac_band[ik]] *
                                              velocity[ispin][ik][alpha](iunocc, iocc)) /
                                             (omega_ev * omega_ev + egap * egap)
                                      << std::endl;
                        }*/
                        if (debug)
                        {
                            if (alpha == 0 && iomega == 0 && mu == 0)
                            {
                                std::complex<double> test =
                                    conj(Ctri_mn[mu][iocc][iunocc][kfrac_band[ik]] *
                                         velocity[ispin][ik][alpha](iunocc, iocc)) /
                                    (omega_ev * omega_ev + egap * egap);
                                if (iocc == 0 && iunocc == 10)
                                {
                                    std::cout
                                        << "C,p: " << Ctri_mn[mu][iocc][iunocc][kfrac_band[ik]]
                                        << "," << velocity[ispin][ik][alpha](iunocc, iocc)
                                        << std::endl;
                                }
                                test_tot += test;
                            }
                        }
                    }
                    else if (iunocc < nocc && iocc >= nocc)
                    {
                        // for metal
                        wing_term += 0.0 * this->Ctri_mn[mu][iocc][iunocc][kfrac_band[ik]] *
                                     velocity[ispin][ik][alpha](iunocc, iocc) /
                                     (omega_ev * omega_ev + egap * egap);
                    }
                }
            }
            if (debug)
            {
                if (alpha == 0 && iomega == 0 && mu == 0)
                {
                    std::cout << "sum over mn: " << ik << "," << test_tot << std::endl;
                }
            }
        }

        /*if (alpha == 0 && lambda == 0 && iomega == 0)
        {
            std::cout << "wing(x,l=0,w=0,mu): " << mu << ", " << wing_term << std::endl;
        }*/
        // wing_term += conj(this->Coul_vector.at(lambda).at(mu)) * tmp;
    }
    return wing_term;
};

void diele_func::init_Cs(const librpa_int::Cs_LRI &Cs_data)
{
    using global::profiler;

    profiler.start("init_Cs");
    using RI::Tensor;
    const int n_atom = Cs_data.data_libri.size();

    for (int ik = 0; ik != nk; ik++)
    {
        for (int I = 0; I != n_atom; I++)
        {
            for (int J = 0; J != n_atom; J++)
            {
                int n_mu_I = atomic_basis_abf_.get_atom_nb(I);
                int n_ao_I = atomic_basis_wfc_.get_atom_nb(I);
                int n_ao_J = atomic_basis_wfc_.get_atom_nb(J);

                Vector3_Order<double> k_frac = kfrac_band[ik];
                const std::array<double, 3> k_array = {k_frac.x, k_frac.y, k_frac.z};
                size_t total = n_ao_I * n_ao_J * n_mu_I;
                std::complex<double> *Cs_in = new std::complex<double>[total]();
                matrix_m<std::complex<double>> mat(n_ao_I * n_ao_J, n_mu_I, Cs_in, MAJOR::ROW,
                                                   MAJOR::COL);
                const std::initializer_list<std::size_t> shape{static_cast<std::size_t>(n_mu_I),
                                                               static_cast<std::size_t>(n_ao_I),
                                                               static_cast<std::size_t>(n_ao_J)};
                this->Ctri_ij.data_libri[I][{J, k_array}] =
                    RI::Tensor<std::complex<double>>(shape, mat.sptr());
                delete[] Cs_in;
            }
        }
    }

    this->Ctri_mn.resize(n_abf);
    for (int mu = 0; mu != n_abf; mu++)
    {
        this->Ctri_mn.at(mu).resize(n_states);
        for (int m = 0; m != n_states; m++)
        {
            this->Ctri_mn.at(mu).at(m).resize(n_states);
            for (int n = 0; n != n_states; n++)
            {
                for (int ik = 0; ik != nk; ik++)
                {
                    this->Ctri_mn.at(mu).at(m).at(n).insert(std::make_pair(kfrac_band[ik], 0.0));
                }
            }
        }
    }
    // std::cout << "* Success: Initialize Ctri_ij and Ctri_mn.\n";
    profiler.stop("init_Cs");
};

void diele_func::FT_R2k(const librpa_int::Cs_LRI &Cs_data)
{
    using global::profiler;

    profiler.start("fourier_r2k");
    const int n_atom = Cs_data.data_libri.size();
    const bool use_shrink_abfs = false;  // TODO: replace with the actual shrink tag
    // std::cout << "Number of atom: " << n_atom << std::endl;

    for (int ik = 0; ik != nk; ik++)
    {
        for (int I = 0; I != n_atom; I++)
        {
            for (int J = 0; J != n_atom; J++)
            {
                int n_mu_I = atomic_basis_abf_.get_atom_nb(I);
                int n_ao_I = atomic_basis_wfc_.get_atom_nb(I);
                int n_ao_J = atomic_basis_wfc_.get_atom_nb(J);
#pragma omp parallel for schedule(dynamic) collapse(3)
                for (int mu = 0; mu != n_mu_I; mu++)
                {
                    for (int i = 0; i != n_ao_I; i++)
                    {
                        for (int j = 0; j != n_ao_J; j++)
                        {
                            Vector3_Order<double> k_frac = kfrac_band[ik];
                            const std::array<double, 3> k_array = {k_frac.x, k_frac.y, k_frac.z};
                            if (use_shrink_abfs)
                            {
                                this->Ctri_ij.data_libri[I][{J, k_array}](mu, i, j) =
                                    compute_Cijk(Cs_shrinked_data, mu, I, i, J, j, ik);
                            }
                            else
                            {
                                this->Ctri_ij.data_libri[I][{J, k_array}](mu, i, j) =
                                    compute_Cijk(Cs_data, mu, I, i, J, j, ik);
                            }

                            /*if (ik == 19 && I == 1 && J == 0 && i == 0 && j == 1)
                            {
                                std::cout << "Cij: " << mu << ", "
                                          << this->Ctri_ij.data_libri[I][{J, k_array}](mu, i, j)
                                          << std::endl;
                            }*/
                        }
                    }
                }
            }
        }
    }
    /* Vector3_Order<int> period{kv_nmp[0], kv_nmp[1], kv_nmp[2]};
    auto Rlist = construct_R_grid(period);
    std::cout << "Number of Bvk cell: " << Rlist.size() << std::endl; */
    std::cout << "* Success: Fourier transform from Cs(R) to Cs(k)." << std::endl;
    profiler.stop("fourier_r2k");
};

std::complex<double> diele_func::compute_Cijk(const librpa_int::Cs_LRI &Cs_data, int mu, int I, int i, int J, int j, int ik)
{
    using librpa_int::TWO_PI;

    std::complex<double> Cijk = 0.0;
    Vector3_Order<int> period{3, 3, 3}; // temp
    auto Rlist = construct_R_grid(period);
    Vector3_Order<double> k_frac = kfrac_band[ik];
    for (auto outer : Cs_data.data_libri.at(I))
    {
        auto J_Ra = outer.first;
        auto Ra = J_Ra.second;
        Vector3_Order<int> R = {Ra[0], Ra[1], Ra[2]};
        double ang = k_frac * R * TWO_PI;
        std::complex<double> kphase = std::complex<double>(cos(ang), sin(ang));
        if (J_Ra.first == J)
        {
            // std::cout << I << "," << J << "," << Ra[0] << "," << Ra[1] << "," << Ra[2] <<
            // std::endl;
            Cijk += kphase * Cs_data.data_libri.at(I).at({J, Ra})(mu, i, j);
        }
    }
    return Cijk;
};

void diele_func::Cs_ij2mn()
{
    profiler.start("transform_Cs_NAO_to_KS");
#pragma omp parallel for schedule(dynamic) collapse(4)
    for (int ik = 0; ik != nk; ik++)
    {
        for (int m = 0; m != n_states; m++)
        {
            for (int n = 0; n != n_states; n++)
            {
                for (int mu = 0; mu != n_abf; mu++)
                {
                    this->Ctri_mn.at(mu).at(m).at(n).at(kfrac_band[ik]) =
                        compute_Cs_ij2mn(mu, m, n, ik);
                    /*if (ik == 26 && m == 5 && n == 40)
                    {
                        lib_printf("Cmn: %5d, %15.5e, %15.5e\n", mu,
                                   Ctri_mn.at(mu).at(m).at(n).at(kfrac_band[ik]).real(),
                                   Ctri_mn.at(mu).at(m).at(n).at(kfrac_band[ik]).imag());
                    }*/
                }
            }
        }
    }

    std::cout << "* Success: transform of Cs from NAO to KS." << std::endl;
    profiler.stop("transform_Cs_NAO_to_KS");
};

std::complex<double> diele_func::compute_Cs_ij2mn(int mu, int m, int n, int ik)
{
    const std::array<double, 3> k_array = {kfrac_band[ik].x, kfrac_band[ik].y, kfrac_band[ik].z};
    int Mu = atomic_basis_abf_.get_i_atom(mu);
    int mu_local = atomic_basis_abf_.get_local_index(mu, Mu);
    std::complex<double> total = 0.0;
    const int n_atom = Ctri_ij.data_libri.size();
    const int n_ao_Mu = atomic_basis_wfc_.get_atom_nb(Mu);
    const int n_spins = meanfield_df.get_n_spins();
    const int n_soc = meanfield_df.get_n_soc();

    // #pragma omp parallel for schedule(dynamic) collapse(2)
    for (int isp = 0; isp != n_spins; isp++)
    {
        for (int is1 = 0; is1 != n_soc; is1++)
        {
            const ComplexMatrix eigenvectors1 = meanfield_df.get_eigenvectors()[isp][is1][ik];
            for (int is2 = 0; is2 != n_soc; is2++)
            {
                const ComplexMatrix eigenvectors2 = meanfield_df.get_eigenvectors()[isp][is2][ik];
                for (int i = 0; i != n_ao_Mu; i++)
                {
                    for (int J = 0; J != n_atom; J++)
                    {
                        int n_ao_J = atomic_basis_wfc_.get_atom_nb(J);
                        for (int j = 0; j != n_ao_J; j++)
                        {
                            std::complex<double> term1 =
                                conj(eigenvectors1(m, atomic_basis_wfc_.get_global_index(Mu, i))) *
                                Ctri_ij.data_libri[Mu][{J, k_array}](mu_local, i, j) *
                                eigenvectors2(n, atomic_basis_wfc_.get_global_index(J, j));
                            std::complex<double> term2 =
                                eigenvectors1(n, atomic_basis_wfc_.get_global_index(Mu, i)) *
                                conj(Ctri_ij.data_libri[Mu][{J, k_array}](mu_local, i, j)) *
                                conj(eigenvectors2(m, atomic_basis_wfc_.get_global_index(J, j)));
                            // #pragma omp critical
                            total += term1 + term2;
                        }
                    }
                }
            }
        }
    }

    return total;
};

// real double diagonalization
// Note complex diagonalization conserves symmetry much better
void diele_func::get_Xv_real(double vq_threshold, const librpa_int::atpair_k_cplx_mat_t &Vq)
{
    using namespace librpa_int;
    using RI::Tensor;
    using RI::Communicate_Tensors_Map_Judge::comm_map2_first;
    using librpa_int::global::mpi_comm_global_h;

    this->Coul_vector.clear();
    this->Coul_value.clear();
    const double CONE = 1.0;
    std::array<double, 3> qa = {0.0, 0.0, 0.0};
    Vector3_Order<double> q = {0.0, 0.0, 0.0};
    size_t n_singular;
    vec<double> eigenvalues(n_abf);

    const auto &comm_h = mpi_comm_global_h;

    comm_h.barrier();
    BlacsCtxtHandler blacs_h(comm_h.comm);

    ArrayDesc desc_nabf_nabf(blacs_h);
    desc_nabf_nabf.init_square_blk(n_abf, n_abf, 0, 0);
    const auto set_IJ_nabf_nabf =
        get_necessary_IJ_from_block_2D_sy('U', atomic_basis_abf_, desc_nabf_nabf);
    const auto s0_s1 = get_s0_s1_for_comm_map2_first(set_IJ_nabf_nabf);
    auto coul_eigen_block = init_local_mat<double>(desc_nabf_nabf, MAJOR::COL);
    auto coulwc_block = init_local_mat<double>(desc_nabf_nabf, MAJOR::COL);
    coulwc_block.zero_out();
    std::map<int, std::map<std::pair<int, std::array<double, 3>>, RI::Tensor<double>>>
        couleps_libri;
    const int natom = atomic_basis_abf_.n_atoms;
    const auto atpair_local = dispatch_upper_triangular_tasks(
        natom, blacs_h.myid, blacs_h.nprows, blacs_h.npcols,
        blacs_h.myprow, blacs_h.mypcol);
    for (const auto &Mu_Nu : atpair_local)
    {
        const auto Mu = Mu_Nu.first;
        const auto Nu = Mu_Nu.second;
        // ofs_myid << "Mu " << Mu << " Nu " << Nu << endl;
        if (Vq.count(Mu) == 0 || Vq.at(Mu).count(Nu) == 0 ||
            Vq.at(Mu).at(Nu).count(q) == 0)
            continue;
        auto Vq_cpl = *(Vq.at(Mu).at(Nu).at(q));
        // if (Vq.count(Mu) == 0 || Vq.at(Mu).count(Nu) == 0 || Vq.at(Mu).at(Nu).count(q) == 0)
        //     continue;
        // auto Vq_cpl = *(Vq.at(Mu).at(Nu).at(q));
        const auto &Vq0 = std::make_shared<matrix>(Vq_cpl.real());
        // const auto &Vq0 = Vq.at(Mu).at(Nu).at(q);
        const auto n_mu = atomic_basis_abf_.get_atom_nb(Mu);
        const auto n_nu = atomic_basis_abf_.get_atom_nb(Nu);
        std::valarray<double> Vq_va(Vq0->c, Vq0->size);
        auto pvq = std::make_shared<std::valarray<double>>();
        *pvq = Vq_va;
        couleps_libri[Mu][{Nu, qa}] = RI::Tensor<double>({n_mu, n_nu}, pvq);
    }
    const auto IJq_coul = RI::Communicate_Tensors_Map_Judge::comm_map2_first(
        mpi_comm_global_h.comm, couleps_libri, s0_s1.first, s0_s1.second);
    collect_block_from_ALL_IJ_Tensor(coulwc_block, desc_nabf_nabf, atomic_basis_abf_, qa,
                                     true, CONE, IJq_coul, MAJOR::ROW);
    power_hemat_blacs_real(coulwc_block, desc_nabf_nabf, coul_eigen_block, desc_nabf_nabf,
                           n_singular, eigenvalues.c, 1.0, vq_threshold);
    this->n_nonsingular = n_abf - n_singular;
    for (int iv = 1; iv != n_nonsingular; iv++)
    {
        // Here eigen solved by Scalapack is ascending order,
        // however, what we want is descending order.
        this->Coul_value.push_back(eigenvalues.c[iv] + 0.0);  // throw away the largest one
        std::vector<std::complex<double>> newRow;

        for (int jabf = 0; jabf != n_abf; jabf++)
        {
            newRow.push_back(coul_eigen_block(jabf, iv));
        }
        this->Coul_vector.push_back(newRow);
    }

    std::cout << "The largest/smallest eigenvalue of Coulomb matrix(non-singular): "
              << this->Coul_value.front() << ", " << this->Coul_value.back() << std::endl;
    std::cout << "The 1st/2nd/3rd/-1th eigenvalue of Coulomb matrix(Full): " << eigenvalues.c[0]
              << ", " << eigenvalues.c[1] << ", " << eigenvalues.c[2] << ", "
              << eigenvalues.c[n_abf - 1] << std::endl;
    std::cout << "Dim of eigenvectors: " << coul_eigen_block.nr() << ", "
              << coul_eigen_block.nc() << std::endl;
    /*std::cout << "Coulomb vector: lambda=-1" << std::endl;
    for (int j = 0; j != n_abf; j++)
    {
        std::cout << j << "," << coul_eigen_block(j, 0) << std::endl;
    }
    std::cout << "Coulomb vector: lambda=0" << std::endl;
    for (int j = 0; j != n_abf; j++)
    {
        std::cout << j << "," << Coul_vector[0][j] << std::endl;
    }
    std::cout << "Coulomb vector: lambda=1" << std::endl;
    for (int j = 0; j != n_abf; j++)
    {
        std::cout << j << "," << Coul_vector[1][j] << std::endl;
    }*/
    std::cout << "* Success: diagonalize Coulomb matrix in the ABFs repre." << std::endl;
};

// complex diagonalization
void diele_func::get_Xv_cpl(double vq_threshold, const librpa_int::atpair_k_cplx_mat_t &Vq)
{
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
    std::cout << "n_singular: " << n_singular << std::endl;
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

    std::cout << "The largest/smallest eigenvalue of Coulomb matrix(non-singular): "
              << this->Coul_value.front() << ", " << this->Coul_value.back() << std::endl;
    std::cout << "The 1st/2nd/3rd/-1th eigenvalue of Coulomb matrix(Full): "
              << eigenvalues.c[n_abf - 1] << ", " << eigenvalues.c[n_abf - 2] << ", "
              << eigenvalues.c[n_abf - 3] << ", " << eigenvalues.c[0] << std::endl;
    std::cout << "Dim of eigenvectors: " << coul_eigen_block.nr() << ", "
              << coul_eigen_block.nc() << std::endl;
    std::cout << "Coulomb vector: lambda=0" << std::endl;
    std::cout << "Coulomb vector: lambda=-1" << std::endl;
    for (int j = 0; j != n_abf; j++)
    {
        std::cout << j << "," << coul_eigen_block(j, 0) << std::endl;
    }
    std::cout << "Coulomb vector: lambda=0" << std::endl;
    for (int j = 0; j != n_abf; j++)
    {
        std::cout << j << "," << Coul_vector[0][j] << std::endl;
    }
    std::cout << "Coulomb vector: lambda=1" << std::endl;
    for (int j = 0; j != n_abf; j++)
    {
        std::cout << j << "," << Coul_vector[1][j] << std::endl;
    }

    std::cout << "* Success: diagonalize Coulomb matrix in the ABFs repre.\n";
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
    std::cout << "BEGIN test head !!!!!!!!!!" << std::endl;
    for (int iomega = 0; iomega != this->omega.size(); iomega++)
    {
        std::complex<double> df = 0;
        for (int alpha = 0; alpha != 3; alpha++)
        {
            df += this->head.at(iomega)(alpha, alpha);
        }
        std::cout << this->omega[iomega] << " " << df.real() / 3.0 << " " << df.imag() / 3.0
                  << std::endl;
    }
    std::cout << "END test head !!!!!!!!!!" << std::endl;
    // std::exit(0);
};

void diele_func::test_wing()
{
    std::cout << "BEGIN test wing !!!!!!!!!!" << std::endl;
    /*std::cout << "wing(z, mu=0) vs omega" << std::endl;
    for (int iomega = 0; iomega != this->omega.size(); iomega++)
    {
        std::complex<double> df = 0;
        // z direction and the first lambda
        df = this->wing_mu.at(2).at(0).at(iomega);

        std::cout << this->omega[iomega] << " " << df.real() << " " << df.imag() << std::endl;
    }*/
    std::cout << "wing_mu(z, iomega=0) vs mu" << std::endl;
    for (int mu = 0; mu != n_abf; mu++)
    {
        std::complex<double> df = 0;
        // z direction
        df = this->wing_mu.at(0)(mu, 2);

        std::cout << mu << " " << df.real() << " " << df.imag() << std::endl;
    }
    /*std::cout << "wing_lambda(z, iomega=0) vs lambda" << std::endl;
    for (int lambda = 0; lambda != n_nonsingular - 1; lambda++)
    {
        std::complex<double> df = 0;
        // z direction
        df = this->wing.at(0)(lambda, 2);

        std::cout << lambda << " " << df.real() << " " << df.imag() << std::endl;
    }
    std::cout << "END test wing !!!!!!!!!!" << std::endl;
    std::exit(0);
}
    std::cout << "END test wing !!!!!!!!!!" << std::endl;*/
    // std::exit(0);
};

// M
void diele_func::get_body_inv(matrix_m<std::complex<double>> &chi0_block)
{
    using librpa_int::ArrayDesc;
    using librpa_int::global::mpi_comm_global_h;

    mpi_comm_global_h.barrier();
    BlacsCtxtHandler blacs_h(mpi_comm_global_h.comm);
    ArrayDesc desc_nabf_nabf(blacs_h);
    desc_nabf_nabf.init_square_blk(n_nonsingular - 1, n_nonsingular - 1, 0, 0);
    this->body_inv = init_local_mat<complex<double>>(desc_nabf_nabf, MAJOR::COL);

    for (int ilambda = 0; ilambda < n_nonsingular - 1; ilambda++)
    {
        for (int jlambda = 0; jlambda < n_nonsingular - 1; jlambda++)
        {
            body_inv(ilambda, jlambda) = chi0_block(1 + ilambda, 1 + jlambda);
        }
    }
    /*std::cout << "Here body: " << std::endl;
    for (int i = 0; i < 5; i++)
    {
        std::cout << i << "," << body_inv(i, i) << std::endl;
    }*/
    invert_scalapack(body_inv, desc_nabf_nabf);
    /*std::cout << "Here binv: " << std::endl;
    for (int i = 0; i < 5; i++)
    {
        std::cout << i << "," << body_inv(i, i) << std::endl;
    }*/

    // test invert
    /*matrix_m<std::complex<double>> test(2, 2, MAJOR::COL);
    Array_Desc desc_2_2(blacs_ctxt_global_h);
    desc_2_2.init_square_blk(2, 2, 0, 0);
    test = init_local_mat<complex<double>>(desc_2_2, MAJOR::COL);
    test(0, 0) = 1;
    test(0, 1) = 2;
    test(1, 0) = 3;
    test(1, 1) = 4;
    invert_scalapack(test, desc_2_2);
    std::cout << test(0, 0) << "," << test(0, 1) << "," << test(1, 0) << "," << test(1, 1)
              << std::endl;*/
    // test invert

    // std::cout << "* Success: get inverse body of chi0.\n";
};

void diele_func::construct_L(const int ifreq)
{
    profiler.start("cal_L");
    matrix_m<std::complex<double>> tmp(3, 3, MAJOR::COL);
    tmp = head.at(ifreq) - transpose(wing.at(ifreq), true) * body_inv * wing.at(ifreq);
    this->Lind = tmp;
    profiler.stop("cal_L");
};

void diele_func::get_Leb_points()
{
    auto quad_order = lebedev::QuadratureOrder::order_590;
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

void diele_func::cal_eps(const int ifreq)
{
    using global::mpi_comm_global_h;
    using global::profiler;

    profiler.start("cal_inverse_dielectric_matrix");
    const auto &comm_h = mpi_comm_global_h;
    comm_h.barrier();

    BlacsCtxtHandler blacs_h(comm_h.comm);

    ArrayDesc desc_nabf_nabf(blacs_h);
    desc_nabf_nabf.init_square_blk(n_abf, n_abf, 0, 0);
    this->chi0 = init_local_mat<complex<double>>(desc_nabf_nabf, MAJOR::COL);

    const double k_volume = std::abs(pbc_.G.Det());
    this->vol_gamma = k_volume / nk;
    double vol_gamma_numeric = 0.0;
    if (ifreq == 0)
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
    construct_L(ifreq);
    profiler.start("cal_inverse_dielectric_matrix_ij");
#pragma omp parallel for schedule(dynamic) collapse(2)
    for (int i = 0; i != n_nonsingular; i++)
    {
        for (int j = 0; j != n_nonsingular; j++)
        {
            if (i == 0 || j == 0)
            {
                if (i == 0 && j == 0)
                {
                    chi0(i, j) = compute_chi0_inv_00(ifreq);
                }
                else
                    chi0(i, j) = 0.0;
            }
            else
                chi0(i, j) = compute_chi0_inv_ij(ifreq, i - 1, j - 1);
            /*if (i == j)
            {
                chi0(i, j) -= 1.0;
            }*/
        }
    }
    profiler.stop("cal_inverse_dielectric_matrix_ij");
    std::cout << "* Success: calculate average inverse dielectric matrix no." << ifreq + 1 << "."
              << std::endl;
    profiler.stop("cal_inverse_dielectric_matrix");
};

std::complex<double> diele_func::compute_chi0_inv_00(const int ifreq)
{
    std::complex<double> total = 0.0;

    for (int ileb = 0; ileb != qw_leb.size(); ileb++)
    {
        matrix_m<std::complex<double>> q_unit(3, 1, MAJOR::COL);
        q_unit(0, 0) = qx_leb[ileb];
        q_unit(1, 0) = qy_leb[ileb];
        q_unit(2, 0) = qz_leb[ileb];

        auto den = transpose(q_unit, false) * Lind * q_unit;
        total += qw_leb[ileb] * std::pow(q_gamma[ileb], 3) / den(0, 0);
    }
    total *= 1.0 / 3.0 / vol_gamma;

    return total;
};

std::complex<double> diele_func::compute_chi0_inv_ij(const int ifreq, int i, int j)
{
    std::complex<double> total = 0.0;
    std::vector<std::complex<double>> partial_sum(qw_leb.size(), 0.0);
    matrix_m<std::complex<double>> body_inv_i(1, n_nonsingular - 1, MAJOR::COL);
    matrix_m<std::complex<double>> body_inv_j(n_nonsingular - 1, 1, MAJOR::COL);
    for (int ii = 0; ii != n_nonsingular - 1; ii++)
    {
        body_inv_i(0, ii) = body_inv(i, ii);
        body_inv_j(ii, 0) = body_inv(ii, j);
    }
#pragma omp parallel for schedule(dynamic)
    for (int ileb = 0; ileb != qw_leb.size(); ileb++)
    {
        matrix_m<std::complex<double>> q_unit(3, 1, MAJOR::COL);
        q_unit(0, 0) = qx_leb[ileb];
        q_unit(1, 0) = qy_leb[ileb];
        q_unit(2, 0) = qz_leb[ileb];
        auto den = transpose(q_unit, false) * Lind * q_unit;
        auto bwq_i = body_inv_i * wing.at(ifreq) * q_unit;
        auto qwb_j = transpose(q_unit, false) * transpose(wing.at(ifreq), true) * body_inv_j;
        partial_sum[ileb] =
            qw_leb[ileb] * std::pow(q_gamma[ileb], 3) * bwq_i(0, 0) * qwb_j(0, 0) / den(0, 0);
    }
    total = std::accumulate(partial_sum.begin(), partial_sum.end(), std::complex<double>(0.0, 0.0));
    total *= 1.0 / 3.0 / vol_gamma;
    total += body_inv(i, j);

    // if (ifreq == 0 && i == 1 && j == 1) std::cout << "* Success: calculate epsilon_11.\n";
    return total;
};

void diele_func::rewrite_eps(matrix_m<std::complex<double>> &chi0_block, const int ifreq)
{
    get_body_inv(chi0_block);
    cal_eps(ifreq);
    chi0_block = this->chi0;

    // this->chi0.clear(); // quote by chi0_block, should not be clear
    this->body_inv.clear();
    this->Lind.clear();
    /*if (ifreq == 0)
        std::cout << "* Success: replace average inverse dielectric matrix with head and wing.\n";*/
};

}
