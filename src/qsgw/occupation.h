#pragma once

#include "../core/meanfield.h"

#include <vector>

namespace librpa_int
{
namespace qsgw
{

struct OccupationSettings
{
    double temperature_kelvin = 0.0;
    double degeneracy_tolerance_ha = 1.0e-10;
    double electron_tolerance = 1.0e-12;
};

struct OccupationResult
{
    double chemical_potential = 0.0;
    double electron_count = 0.0;
    double vbm = 0.0;
    double cbm = 0.0;
    double gap = 0.0;
    bool metallic = false;
};

double physical_electron_count(
    const MeanField& meanfield,
    const std::vector<double>& kpoint_weights,
    double tolerance = 1.0e-12);

OccupationResult analyze_qsgw_occupations(
    const MeanField& meanfield,
    const std::vector<double>& kpoint_weights,
    double total_electrons,
    const OccupationSettings& settings = {});

OccupationResult update_qsgw_occupations(
    MeanField& live_meanfield,
    const MeanField& reference_meanfield,
    const std::vector<double>& kpoint_weights,
    double total_electrons,
    const OccupationSettings& settings = {});

} // namespace qsgw
} // namespace librpa_int
