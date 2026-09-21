#pragma once

#include <map>

#include "../../math/matrix_m.h"

namespace librpa_int
{
namespace qsgw
{

using SpinKMatrixMap = std::map<int, std::map<int, Matz>>;
using SpinKFrequencyMatrixMap =
    std::map<int, std::map<int, std::map<double, Matz>>>;

} // namespace qsgw
} // namespace librpa_int
