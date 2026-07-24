#ifdef NDEBUG
#undef NDEBUG
#endif

#include "../api/compute_helper.h"
#include "../qsgw/band_bvk_remap.h"

#include <cassert>
#include <iostream>
#include <vector>

using librpa_int::Atoms;
using librpa_int::Matrix3;
using librpa_int::PeriodicBoundaryData;
using librpa_int::Vector3;
using librpa_int::Vector3_Order;
using librpa_int::atom_t;
using librpa_int::qsgw::build_legacy_band_bvk_remap;

namespace
{

void test_displaced_same_atom_positive_corner_matches_legacy_mapping()
{
    constexpr double a = 5.1315533365962613;
    const Matrix3 lattice(
        0.0, a, a,
        a, 0.0, a,
        a, a, 0.0);

    PeriodicBoundaryData pbc;
    pbc.set_latvec({
        0.0, a, a,
        a, 0.0, a,
        a, a, 0.0});
    pbc.set_period(4, 4, 4);

    Atoms atoms;
    atoms.set(
        {14, 14},
        {Vector3<double>{0.0, 0.0, 0.0},
         Vector3<double>{0.5 * a, 0.5 * a, 0.5 * a}},
        lattice);

    const Vector3_Order<int> corner{1, 1, 1};
    const auto upstream =
        librpa_int::api::build_band_bvk_remap(atoms, pbc, 0);
    const auto* upstream_origin =
        upstream.find_R_bvk({atom_t{0}, atom_t{0}}, corner);
    const auto* upstream_displaced =
        upstream.find_R_bvk({atom_t{1}, atom_t{1}}, corner);
    assert(upstream_origin != nullptr);
    assert(upstream_displaced != nullptr);
    assert(upstream_origin->front() == Vector3_Order<int>(-3, 1, 1));
    assert(upstream_displaced->front() == Vector3_Order<int>(-3, 1, 1));

    const auto legacy = build_legacy_band_bvk_remap(atoms, pbc, 0);
    const auto* legacy_origin =
        legacy.find_R_bvk({atom_t{0}, atom_t{0}}, corner);
    const auto* legacy_displaced =
        legacy.find_R_bvk({atom_t{1}, atom_t{1}}, corner);
    assert(legacy_origin != nullptr);
    assert(legacy_origin->front() == Vector3_Order<int>(-3, 1, 1));
    assert(legacy_displaced == nullptr);
}

} // namespace

int main()
{
    test_displaced_same_atom_positive_corner_matches_legacy_mapping();
    std::cout << "test_qsgw_band_bvk_remap: all tests passed\n";
    return 0;
}
