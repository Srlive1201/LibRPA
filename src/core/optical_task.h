#pragma once
#include "dataset.h"
#include "librpa_optical.h"
namespace librpa_int
{
void prepare_native_velocity(Dataset&, const LibrpaOptions*, const LibrpaNativeVelocityInput*, int);
void prepare_supplied_velocity(Dataset&, const LibrpaOptions*, const LibrpaSuppliedVelocityInput*,
                               int);
}  // namespace librpa_int
