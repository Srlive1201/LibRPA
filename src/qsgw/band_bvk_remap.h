#pragma once

#include "../core/atom.h"
#include "../core/geometry.h"
#include "../core/pbc.h"

namespace librpa_int
{
namespace qsgw
{

AtomPairBvKRemap<atom_t> build_legacy_band_bvk_remap(
    const Atoms& atoms,
    const PeriodicBoundaryData& pbc,
    int remap_convention);

} // namespace qsgw
} // namespace librpa_int
