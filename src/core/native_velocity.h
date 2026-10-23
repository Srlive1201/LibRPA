#pragma once

#include <array>
#include <vector>

#include "../mpi/base_mpi.h"
#include "librpa_optical.h"
#include "meanfield.h"
#include "pbc.h"
#include "symmetry_context.h"

namespace librpa_int
{
struct NativeVelocityResult
{
    std::vector<cplxdb> velocity;
    std::array<double, 4> residuals{};
    // Root-owned normalization supplied alongside H/S/r. The API distributes
    // it only if analytic head/wing preparation is requested.
    std::array<double, 2> coulomb_normalization{};
};

NativeVelocityResult build_native_velocity(const MeanField &mf, const PeriodicBoundaryData &pbc,
                                           const MpiCommHandler &comm,
                                           const LibrpaNativeVelocityInput *root_input,
                                           const SymmetryContext *symmetry = nullptr);
}  // namespace librpa_int
