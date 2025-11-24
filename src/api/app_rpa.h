#pragma once
#include <complex>
#include <vector>

namespace librpa_int
{

namespace app
{

void get_rpa_correlation_energy_(
        std::complex<double> &rpa_corr, 
        std::vector<std::complex<double>> &rpa_corr_k_contrib 
        );


} /* end of namespace app */

} /* end of namespace librpa_int */
