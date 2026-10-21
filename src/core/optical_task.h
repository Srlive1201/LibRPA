#pragma once

#include "dataset.h"
#include "librpa_optical.h"

namespace librpa_int
{

void prepare_native_velocity(Dataset &dataset, const LibrpaOptions *options,
                             const LibrpaNativeVelocityInput *root_input, int prepare_optical);
void prepare_supplied_velocity(Dataset &dataset, const LibrpaOptions *options,
                               const LibrpaSuppliedVelocityInput *input, int need_wing);

}  // namespace librpa_int
