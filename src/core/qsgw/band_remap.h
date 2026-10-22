#pragma once

#include "../atom.h"
#include "../geometry.h"
#include "../pbc.h"

namespace librpa_int::qsgw
{
AtomPairBvKRemap<atom_t> build_qsgw_band_bvk_remap(
    const Atoms& atoms, const PeriodicBoundaryData& pbc, int remap_convention);
}  // namespace librpa_int::qsgw
