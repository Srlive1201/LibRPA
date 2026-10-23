#pragma once

#include <array>
#include <map>

#include "matrix_m.h"
#include "vector3_order.h"

namespace librpa_int
{

// Dense blocks with independent lattice-translation supports. Storage order,
// physical units and Fourier-completion conventions belong to the caller.
using LocalizedMatrix = std::map<Vector3_Order<int>, Matz>;

struct LocalizedAoOperators
{
    int dimension = 0;
    LocalizedMatrix overlap;
    std::array<LocalizedMatrix, 3> position;
};

}  // namespace librpa_int
