#include "convergence.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace librpa_int
{
namespace qsgw
{

EigenvalueSnapshot eigenvalue_snapshot(const MeanField& meanfield)
{
    if (!meanfield.initialized())
    {
        throw std::invalid_argument(
            "QSGW convergence mean field is not initialized");
    }
    EigenvalueSnapshot result;
    result.reserve(static_cast<std::size_t>(meanfield.get_n_spins()));
    for (const matrix& eigenvalues : meanfield.get_eigenvals())
    {
        matrix copy(eigenvalues.nr, eigenvalues.nc, false);
        for (int index = 0; index < eigenvalues.size; ++index)
        {
            if (!std::isfinite(eigenvalues.c[index]))
            {
                throw std::invalid_argument(
                    "QSGW convergence eigenvalue is non-finite");
            }
            copy.c[index] = eigenvalues.c[index];
        }
        result.push_back(std::move(copy));
    }
    return result;
}

double max_eigenvalue_change(
    const MeanField& meanfield,
    const EigenvalueSnapshot& previous)
{
    if (!meanfield.initialized() ||
        previous.size() !=
            static_cast<std::size_t>(meanfield.get_n_spins()))
    {
        throw std::invalid_argument(
            "QSGW convergence snapshot has an incompatible spin layout");
    }
    double result = 0.0;
    for (int spin = 0; spin < meanfield.get_n_spins(); ++spin)
    {
        const matrix& current = meanfield.get_eigenvals()[spin];
        const matrix& old = previous[static_cast<std::size_t>(spin)];
        if (current.nr != old.nr || current.nc != old.nc)
        {
            throw std::invalid_argument(
                "QSGW convergence snapshot has incompatible dimensions");
        }
        for (int index = 0; index < current.size; ++index)
        {
            if (!std::isfinite(current.c[index]) ||
                !std::isfinite(old.c[index]))
            {
                throw std::invalid_argument(
                    "QSGW convergence eigenvalue is non-finite");
            }
            result = std::max(
                result, std::abs(current.c[index] - old.c[index]));
        }
    }
    return result;
}

bool qsgw_iteration_converged(const int iteration,
                              const int minimum_iterations,
                              const double maximum_change,
                              const double tolerance)
{
    if (iteration < 1 || minimum_iterations < 1 ||
        maximum_change < 0.0 || !std::isfinite(maximum_change) ||
        !(tolerance > 0.0) || !std::isfinite(tolerance))
    {
        throw std::invalid_argument(
            "QSGW convergence inputs are invalid");
    }
    return iteration >= minimum_iterations && maximum_change < tolerance;
}

} // namespace qsgw
} // namespace librpa_int
