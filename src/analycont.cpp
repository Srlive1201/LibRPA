/*!
 * @author    Min-Ye Zhang
 * @date      2024-04-25
 */
#include "analycont.h"

#include <cassert>
// #include <iostream>

#include "complexmatrix.h"
// #include "stl_io_helper.h"

namespace LIBRPA
{

AnalyContPade::AnalyContPade(int n_pars_in, const std::vector<cplxdb> &xs, const std::vector<cplxdb> &data)
    : n_pars(n_pars_in)
{
    int n_data = data.size();
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

}
