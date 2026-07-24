#pragma once

#include "../math/matrix_m.h"
#include "../math/vector3_order.h"

#include <map>

namespace librpa_int
{
namespace qsgw
{

using SpinKMatrixMap = std::map<int, std::map<int, Matz>>;
using SpinKFrequencyMatrixMap =
    std::map<int, std::map<int, std::map<double, Matz>>>;
using RealSpaceMatrixMap = std::map<Vector3_Order<int>, Matz>;
using SpinRMatrixMap = std::map<int, RealSpaceMatrixMap>;

} // namespace qsgw
} // namespace librpa_int
