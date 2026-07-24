#ifdef NDEBUG
#undef NDEBUG
#endif

#include "../qsgw/projection_target.h"

#include <cassert>
#include <complex>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

using librpa_int::MeanField;
using librpa_int::Vector3_Order;
using librpa_int::qsgw::validate_projection_target;

namespace
{

template <typename Function>
void assert_throws(Function&& function)
{
    bool threw = false;
    try
    {
        function();
    }
    catch (const std::exception&)
    {
        threw = true;
    }
    assert(threw);
}

MeanField make_meanfield(const int n_kpoints, const int n_bands,
                         const int n_aos = 4, const int n_spinors = 1)
{
    MeanField result(1, n_kpoints, n_bands, n_aos, n_spinors);
    for (int spinor = 0; spinor < n_spinors; ++spinor)
    {
        for (int kpoint = 0; kpoint < n_kpoints; ++kpoint)
        {
            auto& wfc = result.get_eigenvectors()[0][spinor][kpoint];
            wfc.create(n_bands, n_aos);
            for (int band = 0; band < n_bands; ++band)
            {
                for (int ao = 0; ao < n_aos; ++ao)
                {
                    wfc(band, ao) = std::complex<double>(
                        100.0 * spinor + 10.0 * kpoint + band + 0.1 * ao,
                        -0.01 * (spinor + band + ao));
                }
            }
        }
    }
    return result;
}

std::vector<Vector3_Order<double>> make_kpoints(const int count)
{
    std::vector<Vector3_Order<double>> result;
    for (int index = 0; index < count; ++index)
    {
        result.push_back({static_cast<double>(index) / count, 0.0, 0.0});
    }
    return result;
}

void test_grid_and_band_targets_may_have_independent_layouts()
{
    const MeanField grid = make_meanfield(2, 3);
    const MeanField band = make_meanfield(5, 6);

    const auto grid_shape = validate_projection_target(
        grid, make_kpoints(2), 1, 1, 4, "grid");
    const auto band_shape = validate_projection_target(
        band, make_kpoints(5), 1, 1, 4, "band");

    assert(grid_shape.n_kpoints == 2);
    assert(grid_shape.n_bands == 3);
    assert(band_shape.n_kpoints == 5);
    assert(band_shape.n_bands == 6);
}

void test_missing_or_mismatched_target_data_is_rejected()
{
    MeanField target = make_meanfield(2, 3);
    target.get_eigenvectors()[0][0].erase(1);
    assert_throws([&] {
        validate_projection_target(
            target, make_kpoints(2), 1, 1, 4, "grid");
    });

    const MeanField valid = make_meanfield(2, 3);
    assert_throws([&] {
        validate_projection_target(
            valid, make_kpoints(1), 1, 1, 4, "grid");
    });
    assert_throws([&] {
        validate_projection_target(
            valid, make_kpoints(2), 2, 1, 4, "grid");
    });
    assert_throws([&] {
        validate_projection_target(
            valid, make_kpoints(2), 1, 1, 5, "grid");
    });

    const MeanField spinor_target = make_meanfield(2, 3, 4, 2);
    assert_throws([&] {
        validate_projection_target(
            spinor_target, make_kpoints(2), 1, 1, 4, "band");
    });
}

void test_nonfinite_wavefunction_is_rejected()
{
    MeanField target = make_meanfield(1, 2);
    target.get_eigenvectors()[0][0][0](0, 0) =
        std::complex<double>(std::numeric_limits<double>::quiet_NaN(), 0.0);
    assert_throws([&] {
        validate_projection_target(
            target, make_kpoints(1), 1, 1, 4, "grid");
    });
}

} // namespace

int main()
{
    test_grid_and_band_targets_may_have_independent_layouts();
    test_missing_or_mismatched_target_data_is_rejected();
    test_nonfinite_wavefunction_is_rejected();
    std::cout << "test_qsgw_projection_target: all tests passed\n";
    return 0;
}
