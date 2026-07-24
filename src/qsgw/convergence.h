#pragma once

#include "../core/meanfield.h"

#include <vector>

namespace librpa_int
{
namespace qsgw
{

using EigenvalueSnapshot = std::vector<matrix>;

EigenvalueSnapshot eigenvalue_snapshot(const MeanField& meanfield);

double max_eigenvalue_change(
    const MeanField& meanfield,
    const EigenvalueSnapshot& previous);

bool qsgw_iteration_converged(int iteration,
                              int minimum_iterations,
                              double maximum_change,
                              double tolerance);

} // namespace qsgw
} // namespace librpa_int
