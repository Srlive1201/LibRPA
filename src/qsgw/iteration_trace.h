#pragma once

#include "matrix_map.h"
#include "mixing.h"

#include "../core/meanfield.h"

#include <iosfwd>
#include <string>
#include <vector>

namespace librpa_int
{
namespace qsgw
{

enum class IterationChannel
{
    Grid = 0,
    Band = 1,
    Headwing = 2,
};

struct IterationSummary
{
    int iteration = 0;
    double maximum_eigenvalue_change_ev = 0.0;
    double residual_l2_ha = 0.0;
    double residual_max_ha = 0.0;
    double fermi_energy_ev = 0.0;
    double gap_ev = 0.0;
    double electron_count = 0.0;
    bool has_mixing_decision = true;
    MixingMode requested_mode = MixingMode::Linear;
    MixingMode applied_mode = MixingMode::Linear;
    double beta = 0.2;
    bool fell_back = false;
    double reciprocal_condition = 1.0;
    std::vector<double> coefficients;
    bool converged = false;
    std::string fallback_reason;
};

void write_iteration_summary_header(std::ostream& output);
void write_iteration_summary(std::ostream& output,
                             const IterationSummary& summary);

void write_eigenvalue_trace_header(std::ostream& output);
void write_eigenvalue_trace(
    std::ostream& output,
    int iteration,
    IterationChannel channel,
    const MeanField& meanfield,
    const std::vector<Vector3_Order<double>>& kpoints);

void write_matrix_trace_header(std::ostream& output);
void write_matrix_component_trace(
    std::ostream& output,
    int iteration,
    IterationChannel channel,
    const std::string& component,
    const SpinKMatrixMap& matrices);
void write_frequency_matrix_component_trace(
    std::ostream& output,
    int iteration,
    IterationChannel channel,
    const std::string& component,
    const SpinKFrequencyMatrixMap& matrices);
void write_occupation_trace(
    std::ostream& output,
    int iteration,
    IterationChannel channel,
    const MeanField& meanfield);
void write_scalar_component_trace(
    std::ostream& output,
    int iteration,
    IterationChannel channel,
    const std::string& component,
    double value);
void write_wavefunction_trace(
    std::ostream& output,
    int iteration,
    IterationChannel channel,
    const MeanField& meanfield,
    const std::string& component_prefix = "wfc");
void write_velocity_trace(
    std::ostream& output,
    int iteration,
    IterationChannel channel,
    const std::vector<std::vector<std::vector<ComplexMatrix>>>& velocity,
    const std::string& component_prefix = "velocity");

} // namespace qsgw
} // namespace librpa_int
