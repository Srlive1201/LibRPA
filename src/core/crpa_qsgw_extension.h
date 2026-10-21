#pragma once

#include <functional>

#include "crpa_task.h"
#include "orbital_weights.h"
#include "qsgw/task_extension.h"

namespace librpa_int
{
struct CrpaQsgwControls
{
    int interval = 1;
    bool track_orbitals = false, output_analysis = true;
    std::string tracking_definition;
};

using CrpaQsgwOutputSink = std::function<void(const CrpaTaskResult&, const OrbitalWeightResult*,
                                              int state, int gw_update, bool final)>;

//! Shares the live GW response while transporting the cRPA parent frame.
class CrpaQsgwExtension final : public qsgw::TaskExtension
{
public:
    CrpaQsgwExtension(CrpaTaskInput input, CrpaQsgwControls controls, CrpaQsgwOutputSink sink);
    void initialize(Dataset&, const LibrpaOptions&, int completed_iteration) override;
    bool wants_response(int entering_state) const override;
    void evaluate_response(Dataset&, const LibrpaOptions&, const Chi0QMap&, int entering_state,
                           int gw_update) override;
    void state_updated(Dataset&, int completed_iteration) override;
    void read_checkpoint(Dataset&, const std::string& filename, int iteration) override;
    void write_checkpoint(const Dataset&, const std::string& filename,
                          int iteration) const override;
    void evaluate_final(Dataset&, const LibrpaOptions&, int completed_iteration) override;

private:
    void update_tracking(const CrpaTaskResult&);
    void emit(Dataset&, CrpaTaskResult result, int state, int gw_update, bool final);

    CrpaTaskInput input_;
    CrpaQsgwControls controls_;
    CrpaQsgwOutputSink sink_;
    LibrpaOptions options_{};
    CrpaTrackingState tracking_;
};
}  // namespace librpa_int
