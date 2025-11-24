/*!
 * @author    Min-Ye Zhang
 * @date      2024-04-25
 */
#include "../utils/constants.h"
#include "analycont.h"

#include <cassert>
// #include <iostream>

#include "../math/complexmatrix.h"
// #include "../io/stl_io_helper.h"

namespace librpa_int
{

AnalyContPade::AnalyContPade(int n_pars_in, const std::vector<cplxdb> &xs, const std::vector<cplxdb> &data)
    : n_pars(n_pars_in)
{
    const int n_data = data.size();
    std::vector<cplxdb> data_npar;

    assert (n_pars > 0);

    if (n_data <= n_pars)
    {
        // Use all data points
        n_pars = n_data;
        par_x = xs;
        data_npar = data;
    }
    else
    {
        // Select the data points evenly, when number of parameters are fewer than data points
        par_x.resize(n_pars);
        data_npar.resize(n_pars);
        int step = n_data / (n_pars - 1);
        for (int ipar = 0; ipar < n_pars - 1; ipar++)
        {
            par_x[ipar] = xs[ipar * step];
            data_npar[ipar] = data[ipar * step];
        }
        data_npar[n_pars-1] = data[n_data-1];
    }

    // Calculate the continuation coefficients, using Thiel's reciprocal difference method
    ComplexMatrix g(n_pars, n_pars);
    for (int i_par = 0; i_par < n_pars; i_par++)
    {
        g(i_par, 0) = data_npar[i_par];
    }

    for (int i_par = 1; i_par < n_pars; i_par++)
    {
        for (int i = i_par; i < n_pars; i++)
        {
            g(i, i_par) = 
                (g(i_par-1, i_par-1) - g(i, i_par-1)) / ((par_x[i] - par_x[i_par-1]) * g(i, i_par-1));
        }
    }

    par_y.resize(n_pars);
    for (int i_par = 0; i_par < n_pars; i_par++)
    {
        par_y[i_par] = g(i_par, i_par);
    }
}

cplxdb
AnalyContPade::get(const cplxdb &x) const
{
    cplxdb tmp = {1.0, 0.0};

    for (int i_par = n_pars - 1; i_par > 0; i_par--)
    {
        tmp = 1.0 + par_y[i_par] * (x - par_x[i_par-1]) / tmp;
    }
    return par_y[0] / tmp;
}

const std::vector<double> get_specfunc(const AnalyCont &ac, const std::vector<cplxdb> omegas,
                                       const double &ref, const double &e_ks, const double &v_xc,
                                       const double &v_exx,
                                       const double &sigc_omega_imag_shift,
                                       const double &gf_omega_imag_shift)
{
    const int n_freq = omegas.size();
    std::vector<double> sf(n_freq, 0.0);
    for (int i = 0; i < n_freq; i++)
    {
        const auto omega_gf = omegas[i] + cplxdb(0.0, gf_omega_imag_shift);
        const auto omega_sigc = omegas[i] + cplxdb(0.0, sigc_omega_imag_shift);
        cplxdb sf_c = omega_gf - e_ks + v_xc - v_exx - ac.get(omega_sigc - ref);
        sf[i] = - ((1.0 / PI) / sf_c).imag();
    }
    return sf;
}

}
