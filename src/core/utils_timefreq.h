#pragma once

#include "timefreq.h"

#include <string>

#include "meanfield.h"

namespace librpa_int
{

namespace utils
{

TFGrids generate_timefreq_grids(unsigned ngrids, const std::string &grid_type_str, const MeanField &mf);

} /* end of namespace utils */

} /* end of namespace librpa_int */
