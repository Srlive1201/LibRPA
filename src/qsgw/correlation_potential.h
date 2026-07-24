#pragma once

#include "../core/meanfield.h"
#include "../math/matrix_m.h"

#include <map>
#include <vector>

namespace librpa_int
{
namespace qsgw
{

enum class CorrelationPotentialMode
{
    ModeA,
    ModeB,
};

struct CorrelationPotentialSettings
{
    CorrelationPotentialMode mode = CorrelationPotentialMode::ModeB;
    int n_params_anacon = -1;
    int n_params_anacon_resample = -1;
    std::vector<cplxdb> resample_imagfreqs;
};

Matz build_qsgw_correlation_potential(
    const MeanField& meanfield,
    const std::vector<double>& source_frequencies,
    const std::map<double, Matz>& sigma_imaginary_axis,
    int spin,
    int kpoint,
    const CorrelationPotentialSettings& settings = {});

} // namespace qsgw
} // namespace librpa_int
