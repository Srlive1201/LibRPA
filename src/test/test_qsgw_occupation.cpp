// These tests intentionally use assert; keep it active in Release test builds.
#ifdef NDEBUG
#undef NDEBUG
#endif

#include "../qsgw/occupation.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

using librpa_int::MeanField;
using librpa_int::qsgw::OccupationSettings;
using librpa_int::qsgw::analyze_qsgw_occupations;
using librpa_int::qsgw::physical_electron_count;
using librpa_int::qsgw::update_qsgw_occupations;

namespace
{

void assert_close(const double actual, const double expected, const double tolerance = 1.0e-12)
{
    assert(std::abs(actual - expected) < tolerance);
}

template <typename Function>
void assert_invalid_argument(Function&& function)
{
    bool threw = false;
    try
    {
        function();
    }
    catch (const std::invalid_argument&)
    {
        threw = true;
    }
    assert(threw);
}

double stored_total_weight(const MeanField& meanfield)
{
    double result = 0.0;
    for (int spin = 0; spin < meanfield.get_n_spins(); ++spin)
    {
        for (int kpoint = 0; kpoint < meanfield.get_n_kpoints(); ++kpoint)
        {
            for (int band = 0; band < meanfield.get_n_bands(); ++band)
            {
                result += meanfield.get_weight()[spin](kpoint, band);
            }
        }
    }
    return result;
}

void test_analysis_preserves_input_meanfield()
{
    MeanField meanfield(1, 1, 2, 2, 1);
    meanfield.get_weight()[0].zero_out();
    meanfield.get_weight()[0](0, 0) = 2.0;
    meanfield.get_eigenvals()[0](0, 0) = -1.0;
    meanfield.get_eigenvals()[0](0, 1) = 2.0;
    meanfield.get_efermi() = 0.25;

    const auto result = analyze_qsgw_occupations(
        meanfield, {1.0}, 2.0, OccupationSettings{});

    assert_close(meanfield.get_weight()[0](0, 0), 2.0);
    assert_close(meanfield.get_weight()[0](0, 1), 0.0);
    assert_close(meanfield.get_efermi(), 0.25);
    assert_close(result.chemical_potential, 0.5);
    assert_close(result.electron_count, 2.0);
    assert_close(result.gap, 3.0);
}

void test_global_filling_preserves_nonuniform_kpoint_weights()
{
    MeanField reference(1, 2, 2, 2, 1);
    reference.get_weight()[0].zero_out();
    reference.get_weight()[0](0, 0) = 1.5;
    reference.get_weight()[0](1, 0) = 0.5;

    MeanField live = reference;
    live.get_eigenvals()[0](0, 0) = -1.0;
    live.get_eigenvals()[0](0, 1) = 2.0;
    live.get_eigenvals()[0](1, 0) = -0.5;
    live.get_eigenvals()[0](1, 1) = 3.0;

    const std::vector<double> kpoint_weights{0.75, 0.25};
    const auto result = update_qsgw_occupations(
        live, reference, kpoint_weights, 2.0, OccupationSettings{});

    assert_close(live.get_weight()[0](0, 0), 1.5);
    assert_close(live.get_weight()[0](1, 0), 0.5);
    assert_close(live.get_weight()[0](0, 1), 0.0);
    assert_close(live.get_weight()[0](1, 1), 0.0);
    assert_close(stored_total_weight(live), 2.0);
    assert_close(physical_electron_count(live, kpoint_weights), 2.0);
    assert_close(result.electron_count, 2.0);
    assert(result.chemical_potential > -0.5);
    assert(result.chemical_potential < 2.0);
    assert(!result.metallic);
}

void test_physical_electron_count_uses_stored_geometric_weights_once()
{
    MeanField reference(1, 2, 1, 1, 1);
    reference.get_weight()[0].zero_out();
    reference.get_weight()[0](0, 0) = 1.5;

    const std::vector<double> kpoint_weights{0.75, 0.25};
    assert_close(stored_total_weight(reference), 1.5);
    assert_close(physical_electron_count(reference, kpoint_weights), 1.5);
}

void test_global_spin_filling_does_not_fill_one_electron_per_spin()
{
    MeanField reference(2, 1, 2, 2, 1);
    reference.get_weight()[0].zero_out();
    reference.get_weight()[1].zero_out();
    reference.get_weight()[0](0, 0) = 1.0;

    MeanField live = reference;
    live.get_eigenvals()[0](0, 0) = -1.0;
    live.get_eigenvals()[0](0, 1) = 10.0;
    live.get_eigenvals()[1](0, 0) = 0.0;
    live.get_eigenvals()[1](0, 1) = 20.0;

    const auto result = update_qsgw_occupations(
        live, reference, {1.0}, 1.0, OccupationSettings{});

    assert_close(live.get_weight()[0](0, 0), 1.0);
    assert_close(live.get_weight()[0](0, 1), 0.0);
    assert_close(live.get_weight()[1](0, 0), 0.0);
    assert_close(live.get_weight()[1](0, 1), 0.0);
    assert_close(stored_total_weight(live), 1.0);
    assert_close(result.electron_count, 1.0);
    assert(result.chemical_potential > -1.0);
    assert(result.chemical_potential < 0.0);
}

void test_degenerate_frontier_uses_one_capacity_fraction()
{
    MeanField reference(1, 2, 1, 1, 1);
    reference.get_weight()[0].zero_out();
    reference.get_weight()[0](0, 0) = 0.5;
    reference.get_weight()[0](1, 0) = 0.5;

    MeanField live = reference;
    live.get_eigenvals()[0](0, 0) = 0.0;
    live.get_eigenvals()[0](1, 0) = 0.0;

    const auto result = update_qsgw_occupations(
        live, reference, {0.75, 0.25}, 1.0, OccupationSettings{});

    assert_close(live.get_weight()[0](0, 0), 0.75);
    assert_close(live.get_weight()[0](1, 0), 0.25);
    assert_close(result.electron_count, 1.0);
    assert_close(result.chemical_potential, 0.0);
    assert_close(result.gap, 0.0);
    assert(result.metallic);
}

void test_overlapping_reference_manifolds_report_zero_gap()
{
    MeanField reference(1, 2, 2, 2, 1);
    reference.get_weight()[0].zero_out();
    reference.get_weight()[0](0, 0) = 1.0;
    reference.get_weight()[0](1, 0) = 1.0;

    MeanField live = reference;
    live.get_eigenvals()[0](0, 0) = 1.0;
    live.get_eigenvals()[0](0, 1) = 2.0;
    live.get_eigenvals()[0](1, 0) = 0.0;
    live.get_eigenvals()[0](1, 1) = 0.5;

    const auto result = update_qsgw_occupations(
        live, reference, {0.5, 0.5}, 2.0, OccupationSettings{});

    assert_close(result.electron_count, 2.0);
    assert_close(result.vbm, 1.0);
    assert_close(result.cbm, 0.5);
    assert_close(result.gap, 0.0);
    assert_close(result.chemical_potential, 0.75);
    assert(result.metallic);
}

void test_spinor_capacity_is_one_electron_per_state()
{
    MeanField reference(1, 1, 2, 2, 2);
    reference.get_weight()[0].zero_out();
    reference.get_weight()[0](0, 0) = 1.0;

    MeanField live = reference;
    live.get_eigenvals()[0](0, 0) = -1.0;
    live.get_eigenvals()[0](0, 1) = 1.0;

    const auto result = update_qsgw_occupations(
        live, reference, {1.0}, 1.0, OccupationSettings{});

    assert_close(live.get_weight()[0](0, 0), 1.0);
    assert_close(live.get_weight()[0](0, 1), 0.0);
    assert_close(result.electron_count, 1.0);
    assert_close(result.chemical_potential, 0.0);
    assert(!result.metallic);
}

void test_zero_weight_kpoint_has_zero_storage_capacity()
{
    MeanField reference(1, 2, 1, 1, 1);
    reference.get_weight()[0].zero_out();
    reference.get_weight()[0](0, 0) = 2.0;

    MeanField live = reference;
    live.get_eigenvals()[0](0, 0) = 0.0;
    live.get_eigenvals()[0](1, 0) = 0.0;

    const auto result = update_qsgw_occupations(
        live, reference, {1.0, 0.0}, 2.0, OccupationSettings{});

    assert_close(live.get_weight()[0](0, 0), 2.0);
    assert_close(live.get_weight()[0](1, 0), 0.0);
    assert_close(result.electron_count, 2.0);
}

void test_nonfinite_energy_rejection_preserves_live_state()
{
    MeanField reference(1, 1, 2, 2, 1);
    reference.get_weight()[0].zero_out();
    reference.get_weight()[0](0, 0) = 2.0;

    MeanField live = reference;
    live.get_weight()[0](0, 0) = 0.25;
    live.get_weight()[0](0, 1) = 0.75;
    live.get_efermi() = 0.125;
    live.get_eigenvals()[0](0, 0) =
        std::numeric_limits<double>::quiet_NaN();
    live.get_eigenvals()[0](0, 1) = 1.0;

    assert_invalid_argument([&]() {
        update_qsgw_occupations(
            live, reference, {1.0}, 2.0, OccupationSettings{});
    });
    assert_close(live.get_weight()[0](0, 0), 0.25);
    assert_close(live.get_weight()[0](0, 1), 0.75);
    assert_close(live.get_efermi(), 0.125);
}

void test_invalid_reference_rejection_preserves_live_state()
{
    MeanField reference(1, 1, 2, 2, 1);
    reference.get_weight()[0].zero_out();
    reference.get_weight()[0](0, 0) = 2.5;

    MeanField live = reference;
    live.get_weight()[0](0, 0) = 0.5;
    live.get_weight()[0](0, 1) = 0.25;
    live.get_efermi() = -0.25;

    assert_invalid_argument([&]() {
        update_qsgw_occupations(
            live, reference, {1.0}, 2.5, OccupationSettings{});
    });
    assert_close(live.get_weight()[0](0, 0), 0.5);
    assert_close(live.get_weight()[0](0, 1), 0.25);
    assert_close(live.get_efermi(), -0.25);
}

void test_live_reference_alias_is_rejected_without_mutation()
{
    MeanField live(1, 1, 2, 2, 1);
    live.get_weight()[0].zero_out();
    live.get_weight()[0](0, 0) = 2.0;
    live.get_eigenvals()[0](0, 0) = -1.0;
    live.get_eigenvals()[0](0, 1) = 1.0;
    live.get_efermi() = 0.25;

    assert_invalid_argument([&]() {
        update_qsgw_occupations(
            live, live, {1.0}, 2.0, OccupationSettings{});
    });
    assert_close(live.get_weight()[0](0, 0), 2.0);
    assert_close(live.get_weight()[0](0, 1), 0.0);
    assert_close(live.get_efermi(), 0.25);
}

void test_finite_temperature_rejection_preserves_live_state()
{
    MeanField reference(1, 1, 1, 1, 1);
    reference.get_weight()[0](0, 0) = 2.0;
    MeanField live = reference;
    live.get_weight()[0](0, 0) = 0.5;
    live.get_efermi() = -0.5;

    OccupationSettings settings;
    settings.temperature_kelvin = 300.0;
    assert_invalid_argument([&]() {
        update_qsgw_occupations(live, reference, {1.0}, 2.0, settings);
    });
    assert_close(live.get_weight()[0](0, 0), 0.5);
    assert_close(live.get_efermi(), -0.5);
}

} // namespace

int main()
{
    test_analysis_preserves_input_meanfield();
    test_global_filling_preserves_nonuniform_kpoint_weights();
    test_physical_electron_count_uses_stored_geometric_weights_once();
    test_global_spin_filling_does_not_fill_one_electron_per_spin();
    test_degenerate_frontier_uses_one_capacity_fraction();
    test_overlapping_reference_manifolds_report_zero_gap();
    test_spinor_capacity_is_one_electron_per_state();
    test_zero_weight_kpoint_has_zero_storage_capacity();
    test_nonfinite_energy_rejection_preserves_live_state();
    test_invalid_reference_rejection_preserves_live_state();
    test_live_reference_alias_is_rejected_without_mutation();
    test_finite_temperature_rejection_preserves_live_state();
    std::cout << "test_qsgw_occupation: all tests passed\n";
    return 0;
}
