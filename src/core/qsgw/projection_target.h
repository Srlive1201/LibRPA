#pragma once

#include <string>
#include <vector>

#include "../../math/vector3_order.h"
#include "../meanfield.h"

namespace librpa_int
{
namespace qsgw
{

void validate_projection_target(const MeanField& reference,
                                const std::vector<Vector3_Order<double>>& kpoints,
                                int expected_n_spins, int expected_n_spinors, int expected_n_aos,
                                const std::string& label);

} // namespace qsgw
} // namespace librpa_int
