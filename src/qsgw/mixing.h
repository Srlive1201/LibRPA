#pragma once

#include "../math/matrix.h"

#include <optional>
#include <string>
#include <vector>

namespace librpa_int
{
namespace qsgw
{

enum class MixingMode
{
    Linear,
};

struct MixingOptions
{
    MixingMode mode = MixingMode::Linear;
    double beta = 0.2;
};

struct MixingDecision
{
    MixingMode requested_mode = MixingMode::Linear;
    MixingMode applied_mode = MixingMode::Linear;
    double beta = 0.2;
    bool fell_back = false;
    std::string fallback_reason;
    double reciprocal_condition = 1.0;
    std::vector<double> coefficients;
};

struct HamiltonianMixResult
{
    matrix grid;
    std::optional<matrix> band;
    MixingDecision decision;
    double residual_l2 = 0.0;
    double residual_max = 0.0;
};

class HamiltonianMixer
{
public:
    explicit HamiltonianMixer(MixingOptions options = {});

    void initialize(const matrix& grid_input);
    void initialize(const matrix& grid_input, const matrix& band_input);

    HamiltonianMixResult mix(const matrix& grid_output);
    HamiltonianMixResult mix(const matrix& grid_output,
                             const matrix& band_output);

private:
    MixingOptions options_;
    std::optional<matrix> grid_input_;
    std::optional<matrix> band_input_;

    HamiltonianMixResult mix_impl(
        const matrix& grid_output,
        const std::optional<matrix>& band_output);
};

} // namespace qsgw
} // namespace librpa_int
