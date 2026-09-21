#include "hamiltonian_mixing.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace librpa_int
{
namespace qsgw
{
namespace
{

void require_matching_hamiltonians(const SpinKMatrixMap& input,
                                   const SpinKMatrixMap& output)
{
    if (input.empty() || input.size() != output.size())
        throw std::invalid_argument("QSGW Hamiltonian spin maps do not match");
    for (const auto& [spin, by_kpoint] : input)
    {
        const auto out_spin = output.find(spin);
        if (by_kpoint.empty() || out_spin == output.end() ||
            by_kpoint.size() != out_spin->second.size())
            throw std::invalid_argument("QSGW Hamiltonian k-point maps do not match");
        for (const auto& [kpoint, value] : by_kpoint)
        {
            const auto out_k = out_spin->second.find(kpoint);
            if (out_k == out_spin->second.end() || value.nr() <= 0 ||
                value.nr() != value.nc() ||
                value.nr() != out_k->second.nr() ||
                value.nc() != out_k->second.nc())
                throw std::invalid_argument("QSGW Hamiltonian matrix dimensions do not match");
        }
    }
}

} // namespace

SpinKHamiltonianResidual measure_spin_k_hamiltonian_residual(
    const SpinKMatrixMap& output,
    const SpinKMatrixMap& input)
{
    require_matching_hamiltonians(input, output);
    double square_sum = 0.0;
    double maximum = 0.0;
    for (const auto& [spin, by_kpoint] : input)
        for (const auto& [kpoint, value] : by_kpoint)
        {
            const Matz& target = output.at(spin).at(kpoint);
            for (int row = 0; row < value.nr(); ++row)
                for (int column = 0; column < value.nc(); ++column)
                {
                    const double magnitude = std::abs(
                        target(row, column) - value(row, column));
                    square_sum += magnitude * magnitude;
                    maximum = std::max(maximum, magnitude);
                }
        }
    if (!std::isfinite(square_sum))
        throw std::invalid_argument("QSGW Hamiltonian residual is non-finite");
    return {std::sqrt(square_sum), maximum};
}

SpinKMatrixMap mix_spin_k_hamiltonian(
    const SpinKMatrixMap& input,
    const SpinKMatrixMap& output,
    const double beta)
{
    if (!(beta > 0.0 && beta <= 1.0))
        throw std::invalid_argument("QSGW mixing beta must be in (0, 1]");
    require_matching_hamiltonians(input, output);
    if (beta == 1.0) return output;
    SpinKMatrixMap result;
    for (const auto& [spin, by_kpoint] : input)
        for (const auto& [kpoint, value] : by_kpoint)
        {
            const Matz& target = output.at(spin).at(kpoint);
            // Matz copies share storage. Allocate a new matrix so the immutable
            // reference and the driver's current input cannot be changed here.
            Matz mixed(value.nr(), value.nc(), value.major());
            for (int row = 0; row < value.nr(); ++row)
                for (int column = 0; column < value.nc(); ++column)
                    mixed(row, column) = value(row, column) + beta *
                        (target(row, column) - value(row, column));
            result[spin][kpoint] = std::move(mixed);
        }
    return result;
}

} // namespace qsgw
} // namespace librpa_int
