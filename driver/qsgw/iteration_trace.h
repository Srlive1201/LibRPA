#pragma once

#include "../../src/qsgw/matrix_map.h"

#include "../../src/core/meanfield.h"

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
    double beta = 1.0;
    bool converged = false;
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
void write_wavefunction_trace(
    std::ostream& output,
    int iteration,
    IterationChannel channel,
    const MeanField& meanfield,
    const std::string& component_prefix = "wfc");
} // namespace qsgw
} // namespace librpa_int
