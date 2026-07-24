// These tests intentionally use assert; keep it active in Release test builds.
#ifdef NDEBUG
#undef NDEBUG
#endif

#include "../qsgw/mixing.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

using librpa_int::matrix;
using librpa_int::qsgw::HamiltonianMixer;
using librpa_int::qsgw::MixingMode;
using librpa_int::qsgw::MixingOptions;

namespace
{

void assert_close(const double actual, const double expected)
{
    assert(std::abs(actual - expected) < 1.0e-14);
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

void test_default_linear_mixing_updates_grid_and_band_together()
{
    HamiltonianMixer mixer;

    matrix grid0(2, 1, true);
    grid0(0, 0) = 1.0;
    grid0(1, 0) = 2.0;
    matrix band0(1, 1, true);
    band0(0, 0) = 10.0;
    mixer.initialize(grid0, band0);

    matrix grid_output(2, 1, true);
    grid_output(0, 0) = 3.0;
    grid_output(1, 0) = 6.0;
    matrix band_output(1, 1, true);
    band_output(0, 0) = 20.0;

    const auto result = mixer.mix(grid_output, band_output);
    assert(result.decision.requested_mode == MixingMode::Linear);
    assert(result.decision.applied_mode == MixingMode::Linear);
    assert_close(result.decision.beta, 0.2);
    assert(!result.decision.fell_back);
    assert(result.decision.coefficients.size() == 1);
    assert_close(result.decision.coefficients.front(), 1.0);
    assert(result.band.has_value());
    assert_close(result.grid(0, 0), 1.4);
    assert_close(result.grid(1, 0), 2.8);
    assert_close(result.band->operator()(0, 0), 12.0);
}

void test_reported_residual_is_unmixed_even_when_beta_is_one()
{
    MixingOptions options;
    options.beta = 1.0;
    HamiltonianMixer mixer(options);

    matrix grid0(2, 1, true);
    mixer.initialize(grid0);
    matrix output(2, 1, true);
    output(0, 0) = 3.0;
    output(1, 0) = 4.0;

    const auto result = mixer.mix(output);
    assert_close(result.residual_l2, 5.0);
    assert_close(result.residual_max, 4.0);
    assert_close(result.grid(0, 0), 3.0);
    assert_close(result.grid(1, 0), 4.0);
}

void test_rejected_band_shape_does_not_change_linear_state()
{
    HamiltonianMixer tested;
    HamiltonianMixer fresh;

    matrix grid0(2, 1, true);
    matrix band0(1, 1, true);
    tested.initialize(grid0, band0);
    fresh.initialize(grid0, band0);

    matrix rejected_grid(2, 1, true);
    rejected_grid(0, 0) = 1.0;
    matrix rejected_band(2, 1, true);
    assert_throws([&] { tested.mix(rejected_grid, rejected_band); });

    matrix accepted_grid(2, 1, true);
    accepted_grid(1, 0) = 1.0;
    matrix accepted_band(1, 1, true);
    accepted_band(0, 0) = 2.0;
    const auto actual = tested.mix(accepted_grid, accepted_band);
    const auto expected = fresh.mix(accepted_grid, accepted_band);

    assert_close(actual.grid(0, 0), expected.grid(0, 0));
    assert_close(actual.grid(1, 0), expected.grid(1, 0));
    assert_close(actual.band->operator()(0, 0),
                 expected.band->operator()(0, 0));
}

void test_nonfinite_input_or_output_is_rejected_before_state_mutation()
{
    HamiltonianMixer tested;
    HamiltonianMixer fresh;

    matrix grid0(2, 1, true);
    matrix band0(1, 1, true);
    tested.initialize(grid0, band0);
    fresh.initialize(grid0, band0);

    matrix nonfinite_grid(2, 1, true);
    nonfinite_grid(0, 0) = std::numeric_limits<double>::quiet_NaN();
    assert_throws([&] { tested.mix(nonfinite_grid, band0); });

    matrix accepted_grid(2, 1, true);
    accepted_grid(0, 0) = 1.0;
    matrix accepted_band(1, 1, true);
    accepted_band(0, 0) = 2.0;
    const auto actual = tested.mix(accepted_grid, accepted_band);
    const auto expected = fresh.mix(accepted_grid, accepted_band);
    assert_close(actual.grid(0, 0), expected.grid(0, 0));
    assert_close(actual.band->operator()(0, 0),
                 expected.band->operator()(0, 0));
}

void test_overflowing_residual_is_rejected_before_state_mutation()
{
    HamiltonianMixer tested;
    HamiltonianMixer fresh;

    matrix grid0(1, 1, true);
    grid0(0, 0) = std::numeric_limits<double>::max();
    tested.initialize(grid0);
    fresh.initialize(grid0);

    matrix overflowing_output(1, 1, true);
    overflowing_output(0, 0) = -std::numeric_limits<double>::max();
    assert_throws([&] { tested.mix(overflowing_output); });

    const auto actual = tested.mix(grid0);
    const auto expected = fresh.mix(grid0);
    assert_close(actual.grid(0, 0), expected.grid(0, 0));
    assert_close(actual.residual_l2, expected.residual_l2);
    assert_close(actual.residual_max, expected.residual_max);
}

} // namespace

int main()
{
    test_default_linear_mixing_updates_grid_and_band_together();
    test_reported_residual_is_unmixed_even_when_beta_is_one();
    test_rejected_band_shape_does_not_change_linear_state();
    test_nonfinite_input_or_output_is_rejected_before_state_mutation();
    test_overflowing_residual_is_rejected_before_state_mutation();
    std::cout << "test_qsgw_mixing: all tests passed\n";
    return 0;
}
