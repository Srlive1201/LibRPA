#include "mixing.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace librpa_int
{
namespace qsgw
{
namespace
{

void require_same_shape(const matrix& lhs,
                        const matrix& rhs,
                        const char* channel)
{
    if (lhs.nr != rhs.nr || lhs.nc != rhs.nc)
    {
        throw std::invalid_argument(
            std::string("QSGW ") + channel +
            " mixing matrix dimensions do not match");
    }
}

void require_finite_matrix(const matrix& value, const char* channel)
{
    if (value.nr <= 0 || value.nc <= 0)
    {
        throw std::invalid_argument(
            std::string("QSGW ") + channel +
            " mixing matrix must be non-empty");
    }
    for (int index = 0; index < value.size; ++index)
    {
        if (!std::isfinite(value.c[index]))
        {
            throw std::invalid_argument(
                std::string("QSGW ") + channel +
                " mixing matrix contains non-finite data");
        }
    }
}

double residual_l2_norm(const matrix& residual)
{
    double squared_norm = 0.0;
    for (int index = 0; index < residual.size; ++index)
    {
        squared_norm += residual.c[index] * residual.c[index];
    }
    return std::sqrt(squared_norm);
}

double residual_max_norm(const matrix& residual)
{
    double result = 0.0;
    for (int index = 0; index < residual.size; ++index)
    {
        result = std::max(result, std::abs(residual.c[index]));
    }
    return result;
}

matrix linear_mix(const matrix& input,
                  const matrix& output,
                  const double beta)
{
    matrix residual = output - input;
    residual *= beta;
    return input + residual;
}

} // namespace

HamiltonianMixer::HamiltonianMixer(MixingOptions options)
    : options_(options)
{
    if (options_.mode != MixingMode::Linear)
    {
        throw std::invalid_argument(
            "QSGW supports only linear Hamiltonian mixing");
    }
    if (!(options_.beta > 0.0 && options_.beta <= 1.0) ||
        !std::isfinite(options_.beta))
    {
        throw std::invalid_argument(
            "QSGW mixing beta must be finite and in (0, 1]");
    }
}

void HamiltonianMixer::initialize(const matrix& grid_input)
{
    require_finite_matrix(grid_input, "grid input");
    grid_input_ = grid_input;
    band_input_.reset();
}

void HamiltonianMixer::initialize(const matrix& grid_input,
                                  const matrix& band_input)
{
    require_finite_matrix(grid_input, "grid input");
    require_finite_matrix(band_input, "band input");
    grid_input_ = grid_input;
    band_input_ = band_input;
}

HamiltonianMixResult HamiltonianMixer::mix(const matrix& grid_output)
{
    return mix_impl(grid_output, std::nullopt);
}

HamiltonianMixResult HamiltonianMixer::mix(
    const matrix& grid_output,
    const matrix& band_output)
{
    return mix_impl(grid_output, band_output);
}

HamiltonianMixResult HamiltonianMixer::mix_impl(
    const matrix& grid_output,
    const std::optional<matrix>& band_output)
{
    if (!grid_input_)
    {
        throw std::logic_error(
            "QSGW Hamiltonian mixer must be initialized before use");
    }
    require_same_shape(*grid_input_, grid_output, "grid");
    require_finite_matrix(grid_output, "grid output");
    if (band_input_.has_value() != band_output.has_value())
    {
        throw std::invalid_argument(
            "QSGW band mixing input/output presence does not match");
    }
    if (band_input_)
    {
        require_same_shape(*band_input_, *band_output, "band");
        require_finite_matrix(*band_output, "band output");
    }

    const matrix grid_residual = grid_output - *grid_input_;
    require_finite_matrix(grid_residual, "grid residual");
    const double residual_l2 = residual_l2_norm(grid_residual);
    const double residual_max = residual_max_norm(grid_residual);
    if (!std::isfinite(residual_l2) || !std::isfinite(residual_max))
    {
        throw std::invalid_argument(
            "QSGW grid residual norm is not finite");
    }

    matrix mixed_grid =
        linear_mix(*grid_input_, grid_output, options_.beta);
    require_finite_matrix(mixed_grid, "mixed grid");
    std::optional<matrix> mixed_band;
    if (band_input_)
    {
        mixed_band =
            linear_mix(*band_input_, *band_output, options_.beta);
        require_finite_matrix(*mixed_band, "mixed band");
    }

    MixingDecision decision;
    decision.requested_mode = MixingMode::Linear;
    decision.applied_mode = MixingMode::Linear;
    decision.beta = options_.beta;
    decision.coefficients = {1.0};

    grid_input_ = mixed_grid;
    band_input_ = mixed_band;
    return {mixed_grid, mixed_band, decision, residual_l2, residual_max};
}

} // namespace qsgw
} // namespace librpa_int
