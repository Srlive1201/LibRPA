#pragma once
#include <complex>
#include <map>
#include <vector>

#include "complexmatrix.h"
#include "vector3_order.h"

namespace LIBRPA
{

namespace app
{

void get_rpa_correlation_energy_(std::complex<double> &rpa_corr,
                                 std::vector<std::complex<double>> &rpa_corr_k_contrib,
                                 std::map<Vector3_Order<double>, ComplexMatrix> &sinvS,
                                 const std::string &input_dir, const bool use_shrink_abfs = 0);

} /* end of namespace app */

} /* end of namespace LIBRPA */
