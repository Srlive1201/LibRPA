#pragma once

#include <map>
#include <vector>

#include "../math/vector3_order.h"

namespace librpa_int
{

enum class SymmetryQPointRestoreMode
{
    NONE,
    TIME_REVERSAL,
    FULL_CRYSTAL,
};

struct SymmetryQPointView
{
    SymmetryQPointRestoreMode restore_mode = SymmetryQPointRestoreMode::NONE;
    std::vector<Vector3_Order<double>> representatives;
    std::map<Vector3_Order<double>, std::vector<Vector3_Order<double>>> members;
    std::map<Vector3_Order<double>, double> weights;
};

}  // namespace librpa_int
