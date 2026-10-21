#pragma once
#include <string>

#include "../dataset.h"
#include "librpa_options.h"

namespace librpa_int::qsgw
{
// Numerical companions share the response and live state; QSGW itself has no
// dependency on any particular companion calculation.
class TaskExtension
{
public:
    virtual ~TaskExtension() = default;
    virtual void initialize(Dataset&, const LibrpaOptions&, int completed_iteration) = 0;
    virtual bool wants_response(int entering_state) const = 0;
    virtual void evaluate_response(Dataset&, const LibrpaOptions&, const Chi0QMap&,
                                   int entering_state, int gw_update) = 0;
    virtual void state_updated(Dataset&, int completed_iteration) = 0;
    virtual void read_checkpoint(Dataset&, const std::string& qsgw_filename, int iteration) = 0;
    virtual void write_checkpoint(const Dataset&, const std::string& qsgw_filename,
                                  int iteration) const = 0;
    virtual void evaluate_final(Dataset&, const LibrpaOptions&, int completed_iteration) = 0;
};
}  // namespace librpa_int::qsgw
