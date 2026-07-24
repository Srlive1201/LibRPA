#ifdef NDEBUG
#undef NDEBUG
#endif

#include "../qsgw/hamiltonian_mixing.h"
#include "../qsgw/hamiltonian_cut.h"

#include <cassert>
#include <cmath>
#include <complex>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

using librpa_int::Matz;
using librpa_int::MeanField;
using librpa_int::cplxdb;
using librpa_int::qsgw::HamiltonianCutMode;
using librpa_int::qsgw::HamiltonianCutOptions;
using librpa_int::qsgw::MixingMode;
using librpa_int::qsgw::MixingOptions;
using librpa_int::qsgw::SpinKHamiltonianMixer;
using librpa_int::qsgw::SpinKMatrixMap;
using librpa_int::qsgw::apply_hamiltonian_cut;
using librpa_int::qsgw::measure_spin_k_hamiltonian_residual;

namespace
{

void assert_close(const cplxdb actual, const cplxdb expected,
                  const double tolerance = 1.0e-13)
{
    assert(std::abs(actual - expected) < tolerance);
}

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

Matz hermitian_2x2(const double d0, const double d1,
                   const cplxdb off_diagonal,
                   const librpa_int::MAJOR major = librpa_int::MAJOR::ROW)
{
    Matz result(2, 2, major);
    result(0, 0) = d0;
    result(1, 1) = d1;
    result(0, 1) = off_diagonal;
    result(1, 0) = std::conj(off_diagonal);
    return result;
}

SpinKMatrixMap deep_copy_map(const SpinKMatrixMap& input)
{
    SpinKMatrixMap result;
    for (const auto& [spin, by_kpoint] : input)
    {
        for (const auto& [kpoint, value] : by_kpoint)
        {
            result[spin][kpoint] = value.copy();
        }
    }
    return result;
}

SpinKMatrixMap make_grid(const double shift)
{
    SpinKMatrixMap result;
    result[0][0] = hermitian_2x2(
        1.0 + shift, 2.0 + shift, {0.25 + shift, -0.5});
    result[0][1] = hermitian_2x2(
        3.0 + shift, 4.0 + shift, {-0.75, 0.2 + shift});
    result[1][0] = hermitian_2x2(
        5.0 + shift, 6.0 + shift, {0.1, 0.3 + shift});
    result[1][1] = hermitian_2x2(
        7.0 + shift, 8.0 + shift, {-0.4 + shift, -0.2});
    return result;
}

SpinKMatrixMap make_band(const double shift)
{
    SpinKMatrixMap result;
    for (int spin = 0; spin < 2; ++spin)
    {
        for (int k = 0; k < 3; ++k)
        {
            Matz value(1, 1);
            value(0, 0) = 10.0 * spin + k + shift;
            result[spin][k] = value;
        }
    }
    return result;
}

void assert_linear_map(const SpinKMatrixMap& mixed,
                       const SpinKMatrixMap& input,
                       const SpinKMatrixMap& output,
                       const double beta)
{
    assert(mixed.size() == input.size());
    for (const auto& [spin, by_k] : input)
    {
        for (const auto& [k, matrix] : by_k)
        {
            const auto& got = mixed.at(spin).at(k);
            const auto& target = output.at(spin).at(k);
            assert(got.nr() == matrix.nr());
            assert(got.nc() == matrix.nc());
            for (int row = 0; row < matrix.nr(); ++row)
            {
                for (int column = 0; column < matrix.nc(); ++column)
                {
                    assert_close(
                        got(row, column),
                        (1.0 - beta) * matrix(row, column) +
                            beta * target(row, column));
                }
            }
        }
    }
}

void assert_maps_equal(const SpinKMatrixMap& actual,
                       const SpinKMatrixMap& expected)
{
    assert(actual.size() == expected.size());
    for (const auto& [spin, by_k] : expected)
    {
        assert(actual.count(spin) == 1);
        assert(actual.at(spin).size() == by_k.size());
        for (const auto& [kpoint, matrix] : by_k)
        {
            const auto& got = actual.at(spin).at(kpoint);
            assert(got.nr() == matrix.nr());
            assert(got.nc() == matrix.nc());
            assert(got.major() == matrix.major());
            for (int row = 0; row < matrix.nr(); ++row)
            {
                for (int column = 0; column < matrix.nc(); ++column)
                {
                    assert_close(got(row, column), matrix(row, column));
                }
            }
        }
    }
}

void test_grid_and_band_storage_major_are_independent()
{
    MixingOptions options;
    options.beta = 0.25;
    SpinKHamiltonianMixer mixer(options);

    SpinKMatrixMap grid0;
    grid0[0][0] = hermitian_2x2(
        1.0, 2.0, {0.25, -0.5}, librpa_int::MAJOR::COL);
    SpinKMatrixMap band0;
    band0[0][0] = hermitian_2x2(
        3.0, 4.0, {-0.75, 0.2}, librpa_int::MAJOR::ROW);
    SpinKMatrixMap grid_target;
    grid_target[0][0] = hermitian_2x2(
        5.0, 6.0, {0.5, 0.25}, librpa_int::MAJOR::ROW);
    SpinKMatrixMap band_target;
    band_target[0][0] = hermitian_2x2(
        7.0, 8.0, {-0.1, -0.4}, librpa_int::MAJOR::COL);

    mixer.initialize(grid0, band0);
    const auto result = mixer.mix(grid_target, band_target);

    assert(result.grid.at(0).at(0).major() == librpa_int::MAJOR::COL);
    assert(result.band->at(0).at(0).major() == librpa_int::MAJOR::ROW);
    assert_linear_map(result.grid, grid0, grid_target, 0.25);
    assert_linear_map(*result.band, band0, band_target, 0.25);
}

void test_complex_grid_and_different_band_layout_mix_together()
{
    MixingOptions options;
    options.beta = 0.2;
    SpinKHamiltonianMixer mixer(options);

    const auto grid0 = make_grid(0.0);
    const auto band0 = make_band(0.0);
    const auto grid_target = make_grid(1.5);
    const auto band_target = make_band(5.0);
    mixer.initialize(grid0, band0);

    const auto result = mixer.mix(grid_target, band_target);
    assert(result.decision.applied_mode == MixingMode::Linear);
    assert(result.band.has_value());
    assert_linear_map(result.grid, grid0, grid_target, 0.2);
    assert_linear_map(*result.band, band0, band_target, 0.2);

    for (const auto& [spin, by_k] : result.grid)
    {
        (void)spin;
        for (const auto& [k, matrix] : by_k)
        {
            (void)k;
            assert_close(matrix(0, 1), std::conj(matrix(1, 0)));
            assert_close(matrix(0, 0).imag(), 0.0);
            assert_close(matrix(1, 1).imag(), 0.0);
        }
    }
}

void test_band_channel_never_changes_grid_mixing_history_or_decision()
{
    MixingOptions options;
    options.mode = MixingMode::Linear;
    options.beta = 0.2;
    SpinKHamiltonianMixer grid_only(options);
    SpinKHamiltonianMixer grid_and_band(options);
    const auto grid0 = make_grid(0.0);
    const auto band0 = make_band(0.0);
    grid_only.initialize(grid0);
    grid_and_band.initialize(grid0, band0);

    std::vector<SpinKMatrixMap> grid_targets{
        make_grid(0.7), make_grid(-0.2), make_grid(1.1),
    };
    grid_targets[1].at(0).at(0)(0, 1) += cplxdb(0.15, -0.05);
    grid_targets[1].at(0).at(0)(1, 0) =
        std::conj(grid_targets[1].at(0).at(0)(0, 1));
    grid_targets[2].at(1).at(1)(0, 0) += 0.35;

    for (std::size_t step = 0; step < grid_targets.size(); ++step)
    {
        const auto grid_result = grid_only.mix(grid_targets[step]);
        const auto combined_result = grid_and_band.mix(
            grid_targets[step], make_band(2.0 + 3.0 * step));
        assert_maps_equal(combined_result.grid, grid_result.grid);
        assert_close(combined_result.residual_l2,
                     grid_result.residual_l2);
        assert_close(combined_result.residual_max,
                     grid_result.residual_max);
        assert(combined_result.decision.requested_mode ==
               grid_result.decision.requested_mode);
        assert(combined_result.decision.applied_mode ==
               grid_result.decision.applied_mode);
        assert(combined_result.decision.fell_back ==
               grid_result.decision.fell_back);
        assert_close(combined_result.decision.reciprocal_condition,
                     grid_result.decision.reciprocal_condition);
        assert(combined_result.decision.coefficients.size() ==
               grid_result.decision.coefficients.size());
        for (std::size_t index = 0;
             index < grid_result.decision.coefficients.size(); ++index)
        {
            assert_close(combined_result.decision.coefficients[index],
                         grid_result.decision.coefficients[index]);
        }
    }
}

void test_invalid_band_call_is_transactional()
{
    MixingOptions options;
    options.beta = 0.2;
    const auto grid0 = make_grid(0.0);
    const auto band0 = make_band(0.0);
    const auto grid_target = make_grid(1.0);
    const auto band_target = make_band(2.0);

    SpinKHamiltonianMixer tested(options);
    SpinKHamiltonianMixer fresh(options);
    tested.initialize(grid0, band0);
    fresh.initialize(grid0, band0);

    auto missing_band_k = band_target;
    missing_band_k.at(0).erase(1);
    assert_throws([&] { tested.mix(grid_target, missing_band_k); });

    const auto after_rejection = tested.mix(grid_target, band_target);
    const auto expected = fresh.mix(grid_target, band_target);
    assert_linear_map(after_rejection.grid, grid0, grid_target, 0.2);
    assert_linear_map(*after_rejection.band, band0, band_target, 0.2);
    assert_close(after_rejection.residual_l2, expected.residual_l2);
    assert_close(after_rejection.residual_max, expected.residual_max);
}

void test_nonhermitian_and_nonfinite_maps_are_rejected()
{
    SpinKHamiltonianMixer mixer;
    const auto grid0 = make_grid(0.0);
    mixer.initialize(grid0);

    auto nonhermitian = make_grid(1.0);
    nonhermitian.at(0).at(0)(1, 0) += cplxdb(0.0, 0.25);
    assert_throws([&] { mixer.mix(nonhermitian); });

    auto nonfinite = make_grid(1.0);
    nonfinite.at(0).at(0)(0, 0) =
        std::numeric_limits<double>::quiet_NaN();
    assert_throws([&] { mixer.mix(nonfinite); });
}

void test_residual_uses_complex_frobenius_norm()
{
    SpinKMatrixMap input;
    input[0][0] = Matz(2, 2);
    SpinKMatrixMap output = deep_copy_map(input);
    output[0][0](0, 0) = 3.0;
    output[0][0](0, 1) = {0.0, 4.0};
    output[0][0](1, 0) = {0.0, -4.0};
    const auto residual = measure_spin_k_hamiltonian_residual(output, input);
    assert_close(residual.l2, std::sqrt(41.0));
    assert_close(residual.maximum, 4.0);
}

void test_exact_cut_region_does_not_enter_linear_mixing_residual()
{
    MeanField live(1, 1, 5, 1, 1);
    live.get_efermi() = 0.0;
    live.get_eigenvals()[0](0, 0) = -2.0;
    live.get_eigenvals()[0](0, 1) = -0.5;
    live.get_eigenvals()[0](0, 2) = 0.0;
    live.get_eigenvals()[0](0, 3) = 0.4;
    live.get_eigenvals()[0](0, 4) = 0.8;

    SpinKMatrixMap reference;
    reference[0][0] = Matz(5, 5);
    SpinKMatrixMap raw = deep_copy_map(reference);
    for (int band = 0; band < 5; ++band)
    {
        reference[0][0](band, band) = -2.0 + band;
        raw[0][0](band, band) = 10.0 + band;
    }
    raw[0][0](0, 1) = cplxdb(1.0, 2.0);
    raw[0][0](1, 0) = std::conj(raw[0][0](0, 1));
    raw[0][0](3, 4) = cplxdb(100.0, -50.0);
    raw[0][0](4, 3) = std::conj(raw[0][0](3, 4));

    HamiltonianCutOptions cut_options;
    cut_options.mode = HamiltonianCutMode::ShiftedReferenceDiagonal;
    cut_options.unoccupied_keep = 1;
    cut_options.shift_ha = 20.0;
    const auto cut_reference = apply_hamiltonian_cut(
        reference, reference, live, cut_options);
    const auto cut_raw = apply_hamiltonian_cut(
        raw, reference, live, cut_options);

    MixingOptions mixing_options;
    mixing_options.mode = MixingMode::Linear;
    mixing_options.beta = 0.2;
    SpinKHamiltonianMixer mixer(mixing_options);
    mixer.initialize(cut_reference);
    const auto mixed = mixer.mix(cut_raw);
    const Matz& matrix = mixed.grid.at(0).at(0);

    assert_close(matrix(0, 0), 0.8 * reference.at(0).at(0)(0, 0) +
                                      0.2 * raw.at(0).at(0)(0, 0));
    assert_close(matrix(0, 1), 0.2 * raw.at(0).at(0)(0, 1));
    for (int band = 3; band < 5; ++band)
    {
        assert_close(matrix(band, band),
                     reference.at(0).at(0)(band, band) + 20.0);
    }
    assert_close(matrix(3, 4), 0.0);
    assert_close(matrix(4, 3), 0.0);

    const auto projected = apply_hamiltonian_cut(
        mixed.grid, reference, live, cut_options);
    mixer.initialize(projected);
    const auto repeated = mixer.mix(projected);
    assert_close(repeated.residual_l2, 0.0);
    assert_close(repeated.residual_max, 0.0);
}

} // namespace

int main()
{
    test_complex_grid_and_different_band_layout_mix_together();
    test_grid_and_band_storage_major_are_independent();
    test_band_channel_never_changes_grid_mixing_history_or_decision();
    test_invalid_band_call_is_transactional();
    test_nonhermitian_and_nonfinite_maps_are_rejected();
    test_residual_uses_complex_frobenius_norm();
    test_exact_cut_region_does_not_enter_linear_mixing_residual();
    std::cout << "test_qsgw_hamiltonian_mixing: all tests passed\n";
    return 0;
}
