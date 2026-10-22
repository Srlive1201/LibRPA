#pragma once
#include <string>

#include "task_context.h"
#include "librpa_options.h"

namespace librpa_int::qsgw
{
// Numerical companions share the response and live state; QSGW itself has no
// dependency on any particular companion calculation.
class TaskExtension
{
public:
    virtual ~TaskExtension() = default;
    virtual void initialize(TaskContext&, const LibrpaOptions&, int completed_iteration) = 0;
    virtual bool wants_response(int entering_state) const = 0;
    virtual void evaluate_response(TaskContext&, const LibrpaOptions&, const Chi0QMap&,
                                   int entering_state, int gw_update) = 0;
    virtual void state_updated(TaskContext&, int completed_iteration) = 0;
    virtual void read_checkpoint(TaskContext&, const std::string& qsgw_filename, int iteration) = 0;
    virtual void write_checkpoint(const TaskContext&, const std::string& qsgw_filename,
                                  int iteration) const = 0;
    virtual void evaluate_final(TaskContext&, const LibrpaOptions&, int completed_iteration) = 0;
};
}  // namespace librpa_int::qsgw
