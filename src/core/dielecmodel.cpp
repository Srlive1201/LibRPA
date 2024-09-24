#include "dielecmodel.h"

#include <cassert>
#include <cmath>

#include "../math/fitting.h"
#include "../math/interpolate.h"
#include "../utils/error.h"
#include "../utils/constants.h"

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

}

void diele_func::cal_head()
{
    using librpa_int::HA2EV;

    int nk = this->meanfield_df.get_n_kpoints();
    int nbands = this->meanfield_df.get_n_bands();
    int nspin = this->meanfield_df.get_n_spins();  //! spin = 1 only
    auto &wg = this->meanfield_df.get_weight()[nspin - 1];
    auto &eigenvalues = this->meanfield_df.get_eigenvals()[nspin - 1];
    auto &velocity = this->meanfield_df.get_velocity()[nspin - 1];
    std::complex<double> tmp;
    int nocc = 0;

    double dielectric_unit = cal_factor();
    for (int i = 0; i != wg.size; i++)
    {
        if (wg.c[i] == 0.)
        {
            nocc = i;
            break;
        }
    }
    init_head();

    for (int ik = 0; ik != nk; ik++)
    {
        for (int iocc = 0; iocc != nocc; iocc++)
        {
            for (int iunocc = nocc; iunocc != nbands; iunocc++)
            {
                double egap = (eigenvalues(ik, iocc) - eigenvalues(ik, iunocc)) * HA2EV;
                for (int alpha = 0; alpha != 3; alpha++)
                {
                    for (int beta = 0; beta != 3; beta++)
                    {
                        for (int iomega = 0; iomega != this->omega.size(); iomega++)
                        {
                            double omega_ev = this->omega[iomega] * HA2EV;
                            tmp = 2.0 * velocity[ik][alpha](iunocc, iocc) *
                                  velocity[ik][beta](iocc, iunocc) /
                                  (egap * egap + omega_ev * omega_ev) / egap;
                            this->head.at(alpha).at(beta).at(iomega) -= tmp;
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
                if (nspin == 1)
                {
                    this->head.at(alpha).at(beta).at(iomega) *= dielectric_unit * 2;
                }
                else if (nspin == 4)
                    this->head.at(alpha).at(beta).at(iomega) *= dielectric_unit;
                if (alpha == beta)
                {
                    this->head.at(alpha).at(beta).at(iomega) += std::complex<double>(1.0, 0.0);
                }
            }
        }
    }
};

void diele_func::cal_wing() {

};

double diele_func::cal_factor()
{
    using librpa_int::TWO_PI;
    using librpa_int::BOHR2ANG;

    const double h_divide_e2 = 25812.80745;
    const double epsilon0 = 8.854187817e-12;
    const double hbar = 1.05457182e-34;
    const double eV = 1.60217662e-19;
    //! Bohr to A
    const double primitive_cell_volume = latvec.Det() * BOHR2ANG * BOHR2ANG * BOHR2ANG;
    // latvec.print();
    double dielectric_unit = TWO_PI / h_divide_e2 / primitive_cell_volume /
                             this->meanfield_df.get_n_kpoints() * 1.0e10 * hbar / epsilon0 / eV;
    return dielectric_unit;
};

void diele_func::init_head()
{
    head.clear();
    wing.clear();
    head.resize(3);
    wing.resize(3);
    for (int ia = 0; ia < 3; ia++)
    {
        head[ia].resize(3);
        wing[ia].resize(3);
        for (int ib = 0; ib < 3; ib++)
        {
            head[ia][ib].resize(this->omega.size());
            wing[ia][ib].resize(this->omega.size());
        }
    }

    for (int alpha = 0; alpha != 3; alpha++)
    {
        for (int beta = 0; beta != 3; beta++)
        {
            for (int iomega = 0; iomega != this->omega.size(); iomega++)
            {
                this->head.at(alpha).at(beta).at(iomega) = std::complex<double>(0.0, 0.0);
            }
        }
    }
}

void diele_func::test_head()
{
    std::cout << "BEGIN test head !!!!!!!!!!" << std::endl;
    for (int iomega = 0; iomega != this->omega.size(); iomega++)
    {
        std::complex<double> df = 0;
        for (int alpha = 0; alpha != 3; alpha++)
        {
            df += this->head.at(alpha).at(alpha).at(iomega);
        }
        std::cout << this->omega[iomega] << " " << df.real() / 3.0 << " " << df.imag() / 3.0
                  << std::endl;
    }
    std::cout << "END test head !!!!!!!!!!" << std::endl;
    std::exit(0);
}
