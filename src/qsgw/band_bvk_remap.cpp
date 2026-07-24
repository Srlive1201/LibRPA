#include "band_bvk_remap.h"

#include "../api/compute_helper.h"
#include "../io/global_io.h"

#include <cmath>

namespace librpa_int
{
namespace qsgw
{

AtomPairBvKRemap<atom_t> build_legacy_band_bvk_remap(
    const Atoms& atoms,
    const PeriodicBoundaryData& pbc,
    const int remap_convention)
{
    auto remap =
        api::build_band_bvk_remap(atoms, pbc, remap_convention);
    if (remap_convention != 0 || pbc.Rlist.empty())
    {
        return remap;
    }

    const Vector3_Order<int> positive_corner{
        (pbc.period.x - 1) / 2,
        (pbc.period.y - 1) / 2,
        (pbc.period.z - 1) / 2};
    int legacy_overrides = 0;
    for (const auto& atom_coord : atoms.coords_frac)
    {
        const auto atom = atom_coord.first;
        const auto& coord = atom_coord.second;
        const bool displaced =
            std::abs(coord.x) + std::abs(coord.y) + std::abs(coord.z) >
            1.0e-12;
        if (displaced &&
            remap.erase_mapping({atom, atom}, positive_corner))
        {
            ++legacy_overrides;
        }
    }
    global::ofs_myid
        << "QSGW legacy band BvK identity overrides: "
        << legacy_overrides << "\n";
    return remap;
}

} // namespace qsgw
} // namespace librpa_int
