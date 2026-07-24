#pragma once

#include "../core/meanfield.h"
#include "../math/vector3_order.h"

#include <string>
#include <vector>

namespace librpa_int
{
namespace qsgw
{

struct ProjectionTargetShape
{
    int n_spins = 0;
    int n_spinors = 0;
    int n_kpoints = 0;
    int n_bands = 0;
    int n_aos = 0;
};

ProjectionTargetShape validate_projection_target(
    const MeanField& reference,
    const std::vector<Vector3_Order<double>>& kpoints,
    int expected_n_spins,
    int expected_n_spinors,
    int expected_n_aos,
    const std::string& label);

} // namespace qsgw
} // namespace librpa_int
