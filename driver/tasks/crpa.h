#pragma once

#include <string>
#include <vector>

#include "librpa_crpa.h"

namespace driver
{
//! Static species/metric input shared by standalone and iterative calculations.
struct CrpaDriverInput
{
    std::vector<std::string> species;
    std::vector<const char*> labels;
    std::vector<double> overlap_packed;
    int n_kpoints = 0, n_aos = 0, n_spinor = 0;
    LibrpaCrpaInput input();
    std::string tracking_definition() const;
};

CrpaDriverInput read_crpa_driver_input();

//! Standalone file/native velocity preset; QSGW already supplies live velocity.
void prepare_crpa_optical_velocity();

//! Collective analysis and root output with propagated root failures.
void write_crpa_state_outputs(const LibrpaCrpaResult* result, CrpaDriverInput& input,
                              const std::string& directory, int state_index = -1,
                              int gw_update = -1, bool final_state = false,
                              bool write_analysis = true,
                              const LibrpaOrbitalResult* supplied_analysis = nullptr);
}  // namespace driver
