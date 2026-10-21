#pragma once
#include <string>

#include "librpa_crpa.h"

namespace driver
{
//! Root-only export; caller propagates failures collectively.
void write_crpa_model(const LibrpaCrpaResult* result, const std::string& directory, int state_index,
                      int gw_update, bool final_state);
}  // namespace driver
