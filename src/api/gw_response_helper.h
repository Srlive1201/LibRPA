#pragma once

#include "dataset.h"
#include "../core/qsgw/task_context.h"

namespace librpa_int::api
{
//! Reuse the existing API GW build for a core task, with an optional live-response consumer.
void build_gw_sigma_for_task(Dataset& dataset, const LibrpaOptions& options,
                             const qsgw::ResponseObserver& on_response = {});
}  // namespace librpa_int::api
