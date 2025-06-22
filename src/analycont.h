/*!
 * @file      analycont.h
 * @brief     Utilities for analytic continuation
 * @author    Min-Ye Zhang
 * @date      2024-04-23
 */
#pragma once
#include <vector>

#include "base_utility.h"

namespace LIBRPA
{

class AnalyCont
{
public:
    virtual cplxdb get(const cplxdb &x) const = 0;
};

class AnalyContPade: public AnalyCont
{
private:
    int n_pars;
    std::vector<cplxdb> par_x;
    std::vector<cplxdb> par_y;

public:
    AnalyContPade(int n_pars_in,
                  const std::vector<cplxdb> &xs,
                  const std::vector<cplxdb> &data);

    /*!
     * @brief get the value of continued function at complex number
     *
     * @param [in]    x    complex argument of function
     *
     * @return    a complex double, the value of function at x
     */
    cplxdb get(const cplxdb &x) const;
};

/*!
 * @brief get the spectral function from Pade approximant
 */
const std::vector<double> get_specfunc(const AnalyCont &ac, const std::vector<cplxdb> omegas,
                                       const double &ref, const double &e_ks, const double &v_xc,
                                       const double &v_exx);
}
