#pragma once

#include "../core/crpa_task.h"
#include "dataset.h"
#include "librpa_crpa.h"

namespace librpa_int::api
{
std::shared_ptr<CrpaTaskContext> make_crpa_task_context(const dataset_ptr_t& dataset);

//! Translate caller-owned C buffers; numerical validation belongs to the task.
CrpaTaskInput unpack_crpa_input(Dataset& dataset, const LibrpaCrpaInput* input);
//! Own the numerical snapshot at the public result-lifetime boundary.
LibrpaCrpaResult* own_crpa_result(CrpaTaskResult result);
}  // namespace librpa_int::api
