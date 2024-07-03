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

class AnalyContPade
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

}
