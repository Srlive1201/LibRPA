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
    std::vector<double> previous_parent, parent_capture;
    int n_parent = 0;

    LibrpaCrpaInput input();
    void advance_tracking(const LibrpaCrpaResult* result);
    //! Transport on the current E/C without a response or self-energy calculation.
    void update_tracking();
    //! Root writes, all ranks may read; the caller propagates errors collectively.
    void write_tracking_checkpoint(const std::string& filename, int iteration,
                                   const std::vector<double>& checkpoint_energies) const;
    void read_tracking_checkpoint(const std::string& filename, int iteration,
                                  const std::vector<double>& checkpoint_energies);
};

CrpaDriverInput read_crpa_driver_input();

//! Standalone file/native velocity preset; QSGW already supplies live velocity.
void prepare_crpa_optical_velocity();

//! Collective analysis and root output with propagated root failures.
void write_crpa_state_outputs(const LibrpaCrpaResult* result, CrpaDriverInput& input,
                              const std::string& directory, int state_index = -1,
                              int gw_update = -1, bool final_state = false,
                              bool write_analysis = true);
}  // namespace driver
