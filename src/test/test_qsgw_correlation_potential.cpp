#ifdef NDEBUG
#undef NDEBUG
#endif

#include "../qsgw/correlation_potential.h"

#include "../core/analycont.h"

#include <cassert>
#include <cmath>
#include <complex>
#include <limits>
#include <map>
#include <stdexcept>
#include <vector>

using librpa_int::AnalyContPade;
using librpa_int::Matz;
using librpa_int::MeanField;
using librpa_int::cplxdb;
using librpa_int::qsgw::CorrelationPotentialMode;
using librpa_int::qsgw::CorrelationPotentialSettings;
using librpa_int::qsgw::build_qsgw_correlation_potential;

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

std::vector<cplxdb> make_imagfreqs(const std::vector<double>& frequencies)
{
    std::vector<cplxdb> result;
    for (const double frequency: frequencies)
    {
        result.emplace_back(0.0, frequency);
    }
    return result;
}

std::map<double, Matz> make_sigma(const std::vector<double>& frequencies)
{
    Matz constant(2, 2);
    constant(0, 0) = 2.123456789e-7;
    constant(0, 1) = cplxdb(1.234567891e-7, 0.456789123e-7);
    constant(1, 0) = std::conj(constant(0, 1));
    constant(1, 1) = 4.234567891e-7;

    Matz slope(2, 2);
    slope(0, 0) = 0.712345678e-7;
    slope(0, 1) = cplxdb(-0.312345678e-7, 0.212345678e-7);
    slope(1, 0) = std::conj(slope(0, 1));
    slope(1, 1) = -0.423456789e-7;

    std::map<double, Matz> sigma;
    for (const double frequency: frequencies)
    {
        sigma[frequency] = constant + cplxdb(0.0, frequency) * slope;
    }
    return sigma;
}

cplxdb evaluate_element_with_upstream_pade(
    const std::vector<double>& source_frequencies,
    const std::map<double, Matz>& sigma,
    const int row,
    const int column,
    const cplxdb target,
    const std::vector<cplxdb>& resample_frequencies,
    const int n_params_anacon,
    const int n_params_anacon_resample)
{
    const auto source_imagfreqs = make_imagfreqs(source_frequencies);
    std::vector<cplxdb> source_data;
    for (const double frequency: source_frequencies)
    {
        source_data.push_back(sigma.at(frequency)(row, column));
    }

    if (resample_frequencies.empty())
    {
        return AnalyContPade(
                   n_params_anacon, source_imagfreqs, source_data).get(target);
    }

    const AnalyContPade source_pade(
        n_params_anacon_resample, source_imagfreqs, source_data);
    std::vector<cplxdb> resampled_data;
    for (const cplxdb frequency: resample_frequencies)
    {
        resampled_data.push_back(source_pade.get(frequency));
    }
    return AnalyContPade(
               n_params_anacon, resample_frequencies, resampled_data).get(target);
}

Matz evaluate_matrix(
    const std::vector<double>& source_frequencies,
    const std::map<double, Matz>& sigma,
    const double energy,
    const std::vector<cplxdb>& resample_frequencies,
    const int n_params_anacon,
    const int n_params_anacon_resample)
{
    Matz result(2, 2);
    for (int row = 0; row < 2; ++row)
    {
        for (int column = 0; column < 2; ++column)
        {
            result(row, column) = evaluate_element_with_upstream_pade(
                source_frequencies, sigma, row, column,
                cplxdb(energy, 0.0), resample_frequencies,
                n_params_anacon, n_params_anacon_resample);
        }
    }
    return 0.5 * (result + transpose(result, true));
}

Matz expected_potential(
    const MeanField& meanfield,
    const std::vector<double>& source_frequencies,
    const std::map<double, Matz>& sigma,
    const CorrelationPotentialMode mode,
    const std::vector<cplxdb>& resample_frequencies,
    const int n_params_anacon,
    const int n_params_anacon_resample)
{
    const double efermi = meanfield.get_efermi();
    const double e0 = meanfield.get_eigenvals()[0](0, 0) - efermi;
    const double e1 = meanfield.get_eigenvals()[0](0, 1) - efermi;
    const Matz sigma0 = evaluate_matrix(
        source_frequencies, sigma, e0, resample_frequencies,
        n_params_anacon, n_params_anacon_resample);
    const Matz sigma1 = evaluate_matrix(
        source_frequencies, sigma, e1, resample_frequencies,
        n_params_anacon, n_params_anacon_resample);
    const Matz sigma_fermi = evaluate_matrix(
        source_frequencies, sigma, 0.0, resample_frequencies,
        n_params_anacon, n_params_anacon_resample);

    Matz expected(2, 2);
    if (mode == CorrelationPotentialMode::ModeB)
    {
        expected(0, 0) = sigma0(0, 0);
        expected(1, 1) = sigma1(1, 1);
        expected(0, 1) = sigma_fermi(0, 1);
        expected(1, 0) = sigma_fermi(1, 0);
        return expected;
    }

    expected(0, 0) = sigma0(0, 0);
    expected(1, 1) = sigma1(1, 1);
    expected(0, 1) = 0.5 * (sigma0(0, 1) + sigma1(0, 1));
    expected(1, 0) = 0.5 * (sigma1(1, 0) + sigma0(1, 0));
    return expected;
}

void run_mode_test(
    const CorrelationPotentialMode mode,
    const std::vector<cplxdb>& resample_frequencies,
    const int n_params_anacon = -1,
    const int n_params_anacon_resample = -1)
{
    MeanField meanfield(1, 1, 2, 2, 1);
    meanfield.get_eigenvals()[0](0, 0) = -0.4;
    meanfield.get_eigenvals()[0](0, 1) = 0.6;
    meanfield.get_efermi() = 0.1;
    const std::vector<double> source_frequencies{0.2, 1.0};
    const auto sigma = make_sigma(source_frequencies);

    CorrelationPotentialSettings settings;
    settings.mode = mode;
    settings.n_params_anacon = n_params_anacon;
    settings.n_params_anacon_resample = n_params_anacon_resample;
    settings.resample_imagfreqs = resample_frequencies;
    const Matz actual = build_qsgw_correlation_potential(
        meanfield, source_frequencies, sigma, 0, 0, settings);
    const Matz expected = expected_potential(
        meanfield, source_frequencies, sigma, mode, resample_frequencies,
        n_params_anacon, n_params_anacon_resample);

    for (int row = 0; row < 2; ++row)
    {
        for (int column = 0; column < 2; ++column)
        {
            assert_close(actual(row, column), expected(row, column));
        }
    }
    assert(std::abs(actual(0, 0)) > 1.0e-8);
}

void test_mode_a_and_b_match_upstream_pade_without_qsgw_rounding()
{
    run_mode_test(CorrelationPotentialMode::ModeA, {});
    run_mode_test(CorrelationPotentialMode::ModeB, {});
}

void test_optional_resampling_matches_the_upstream_two_stage_pade()
{
    run_mode_test(
        CorrelationPotentialMode::ModeB,
        {cplxdb(0.0, 0.3), cplxdb(0.0, 0.8)});
}

void test_two_stage_parameter_counts_follow_upstream_order()
{
    const std::vector<cplxdb> resample_frequencies{
        cplxdb(0.0, 0.3), cplxdb(0.0, 0.65),
        cplxdb(0.0, 0.9)};

    MeanField meanfield(1, 1, 2, 2, 1);
    meanfield.get_eigenvals()[0](0, 0) = -0.4;
    meanfield.get_eigenvals()[0](0, 1) = 0.6;
    meanfield.get_efermi() = 0.1;
    const std::vector<double> source_frequencies{0.2, 0.6, 1.0};
    auto sigma = make_sigma(source_frequencies);
    sigma.at(0.6)(0, 0) += 1.1e-7;
    sigma.at(0.6)(0, 1) += cplxdb(0.4e-7, -0.2e-7);
    sigma.at(0.6)(1, 0) += cplxdb(-0.3e-7, 0.1e-7);
    sigma.at(0.6)(1, 1) -= 0.7e-7;

    CorrelationPotentialSettings settings;
    settings.mode = CorrelationPotentialMode::ModeB;
    settings.n_params_anacon = 2;
    settings.n_params_anacon_resample = 3;
    settings.resample_imagfreqs = resample_frequencies;
    const Matz actual = build_qsgw_correlation_potential(
        meanfield, source_frequencies, sigma, 0, 0, settings);
    const Matz correct = expected_potential(
        meanfield, source_frequencies, sigma, CorrelationPotentialMode::ModeB,
        resample_frequencies, 2, 3);
    const Matz swapped = expected_potential(
        meanfield, source_frequencies, sigma, CorrelationPotentialMode::ModeB,
        resample_frequencies, 3, 2);
    for (int row = 0; row < 2; ++row)
    {
        for (int column = 0; column < 2; ++column)
        {
            assert_close(actual(row, column), correct(row, column));
        }
    }
    assert(std::abs(correct(0, 0) - swapped(0, 0)) > 1.0e-12);
}

void test_nonfinite_matrix_data_is_rejected_instead_of_zeroed()
{
    MeanField meanfield(1, 1, 2, 2, 1);
    meanfield.get_eigenvals()[0](0, 0) = -0.4;
    meanfield.get_eigenvals()[0](0, 1) = 0.6;
    meanfield.get_efermi() = 0.1;
    const std::vector<double> frequencies{0.2, 1.0};
    auto sigma = make_sigma(frequencies);
    sigma.at(0.2)(0, 0) =
        std::numeric_limits<double>::quiet_NaN();
    assert_throws([&] {
        (void)build_qsgw_correlation_potential(
            meanfield, frequencies, sigma, 0, 0, {});
    });
}

void test_invalid_settings_and_frequency_contracts_are_rejected()
{
    MeanField meanfield(1, 1, 2, 2, 1);
    meanfield.get_eigenvals()[0](0, 0) = -0.4;
    meanfield.get_eigenvals()[0](0, 1) = 0.6;
    meanfield.get_efermi() = 0.1;
    const std::vector<double> frequencies{0.2, 1.0};
    const auto sigma = make_sigma(frequencies);

    CorrelationPotentialSettings settings;
    settings.n_params_anacon = 0;
    assert_throws([&] {
        (void)build_qsgw_correlation_potential(
            meanfield, frequencies, sigma, 0, 0, settings);
    });

    settings = {};
    settings.resample_imagfreqs = {cplxdb(0.0, 0.3), cplxdb(0.0, 0.8)};
    settings.n_params_anacon_resample = 0;
    assert_throws([&] {
        (void)build_qsgw_correlation_potential(
            meanfield, frequencies, sigma, 0, 0, settings);
    });

    settings = {};
    settings.resample_imagfreqs = {cplxdb(0.1, 0.3)};
    assert_throws([&] {
        (void)build_qsgw_correlation_potential(
            meanfield, frequencies, sigma, 0, 0, settings);
    });

    assert_throws([&] {
        (void)build_qsgw_correlation_potential(
            meanfield, {1.0, 0.2}, sigma, 0, 0, {});
    });
    assert_throws([&] {
        (void)build_qsgw_correlation_potential(
            meanfield, {0.2, 0.2}, sigma, 0, 0, {});
    });

    auto missing_frequency = sigma;
    missing_frequency.erase(1.0);
    assert_throws([&] {
        (void)build_qsgw_correlation_potential(
            meanfield, frequencies, missing_frequency, 0, 0, {});
    });

    auto wrong_shape = sigma;
    wrong_shape.at(1.0) = Matz(1, 1);
    wrong_shape.at(1.0)(0, 0) = 1.0;
    assert_throws([&] {
        (void)build_qsgw_correlation_potential(
            meanfield, frequencies, wrong_shape, 0, 0, {});
    });

    MeanField nonfinite_meanfield = meanfield;
    nonfinite_meanfield.get_eigenvals()[0](0, 1) =
        std::numeric_limits<double>::infinity();
    assert_throws([&] {
        (void)build_qsgw_correlation_potential(
            nonfinite_meanfield, frequencies, sigma, 0, 0, {});
    });
}

} // namespace

int main()
{
    test_mode_a_and_b_match_upstream_pade_without_qsgw_rounding();
    test_optional_resampling_matches_the_upstream_two_stage_pade();
    test_two_stage_parameter_counts_follow_upstream_order();
    test_nonfinite_matrix_data_is_rejected_instead_of_zeroed();
    test_invalid_settings_and_frequency_contracts_are_rejected();
    return 0;
}
