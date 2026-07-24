#ifdef NDEBUG
#undef NDEBUG
#endif

#include "../qsgw/iteration_trace.h"

#include <cassert>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using librpa_int::MeanField;
using librpa_int::Vector3_Order;
using librpa_int::qsgw::IterationChannel;
using librpa_int::qsgw::IterationSummary;
using librpa_int::qsgw::MixingMode;
using librpa_int::qsgw::SpinKFrequencyMatrixMap;
using librpa_int::qsgw::SpinKMatrixMap;
using librpa_int::qsgw::write_eigenvalue_trace;
using librpa_int::qsgw::write_eigenvalue_trace_header;
using librpa_int::qsgw::write_frequency_matrix_component_trace;
using librpa_int::qsgw::write_iteration_summary;
using librpa_int::qsgw::write_iteration_summary_header;
using librpa_int::qsgw::write_matrix_component_trace;
using librpa_int::qsgw::write_matrix_trace_header;
using librpa_int::qsgw::write_occupation_trace;
using librpa_int::qsgw::write_scalar_component_trace;
using librpa_int::qsgw::write_velocity_trace;
using librpa_int::qsgw::write_wavefunction_trace;

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

void test_summary_is_numeric_through_the_convergence_column()
{
    IterationSummary summary;
    summary.iteration = 3;
    summary.maximum_eigenvalue_change_ev = 0.125;
    summary.residual_l2_ha = 0.25;
    summary.residual_max_ha = 0.5;
    summary.fermi_energy_ev = 1.5;
    summary.gap_ev = 2.5;
    summary.electron_count = 8.0;
    summary.requested_mode = MixingMode::Linear;
    summary.applied_mode = MixingMode::Linear;
    summary.beta = 0.2;
    summary.fell_back = true;
    summary.reciprocal_condition = 1.0e-14;
    summary.coefficients = {0.25, -0.75};
    summary.converged = false;
    summary.fallback_reason = "ill conditioned residual system";

    std::ostringstream output;
    write_iteration_summary_header(output);
    write_iteration_summary(output, summary);

    std::istringstream lines(output.str());
    std::string header;
    std::string row;
    std::getline(lines, header);
    std::getline(lines, row);
    assert(header.find("requested_mode") != std::string::npos);
    assert(header.find("fallback_reason") != std::string::npos);

    std::istringstream values(row);
    int iteration = -1;
    double max_delta = 0.0;
    double residual_l2 = 0.0;
    double residual_max = 0.0;
    double efermi = 0.0;
    double gap = 0.0;
    double electron_count = 0.0;
    int requested = -9;
    int applied = -9;
    double beta = 0.0;
    int fallback = 0;
    double rcond = 0.0;
    double coefficient_l1 = 0.0;
    int coefficient_count = 0;
    int converged = 1;
    std::string coefficients;
    std::string reason;
    values >> iteration >> max_delta >> residual_l2 >> residual_max >>
        efermi >> gap >> electron_count >> requested >> applied >> beta >> fallback >>
        rcond >> coefficient_l1 >> coefficient_count >> converged >>
        coefficients >> reason;
    assert(values);
    assert(iteration == 3);
    assert(std::abs(max_delta - 0.125) < 1.0e-14);
    assert(std::abs(electron_count - 8.0) < 1.0e-14);
    assert(requested == 0);
    assert(applied == 0);
    assert(std::abs(beta - 0.2) < 1.0e-14);
    assert(fallback == 1);
    assert(std::abs(rcond - 1.0e-14) < 1.0e-25);
    assert(std::abs(coefficient_l1 - 1.0) < 1.0e-14);
    assert(coefficient_count == 2);
    assert(converged == 0);
    assert(coefficients.find(',') != std::string::npos);
    assert(reason == "ill_conditioned_residual_system");
}

void test_iteration_zero_and_band_channel_eigenvalues_are_explicit()
{
    MeanField meanfield(1, 2, 2, 1, 1);
    meanfield.get_eigenvals()[0](0, 0) = -0.5;
    meanfield.get_eigenvals()[0](0, 1) = 0.25;
    meanfield.get_eigenvals()[0](1, 0) = -0.4;
    meanfield.get_eigenvals()[0](1, 1) = 0.35;
    const std::vector<Vector3_Order<double>> kpoints{
        {0.0, 0.0, 0.0},
        {0.5, 0.0, 0.0},
    };

    std::ostringstream output;
    write_eigenvalue_trace_header(output);
    write_eigenvalue_trace(
        output, 0, IterationChannel::Band, meanfield, kpoints);

    std::istringstream lines(output.str());
    std::string header;
    std::getline(lines, header);
    assert(header.find("energy_eV") != std::string::npos);

    int row_count = 0;
    std::string row;
    while (std::getline(lines, row))
    {
        if (row.empty()) continue;
        std::istringstream values(row);
        int iteration = -1;
        int channel = -1;
        int spin = -1;
        int kpoint = -1;
        double kx = 0.0;
        double ky = 0.0;
        double kz = 0.0;
        int band = -1;
        double energy_ev = 0.0;
        values >> iteration >> channel >> spin >> kpoint >>
            kx >> ky >> kz >> band >> energy_ev;
        assert(values);
        assert(iteration == 0);
        assert(channel == 1);
        assert(spin == 0);
        assert(kpoint == row_count / 2);
        assert(band == row_count % 2);
        ++row_count;
    }
    assert(row_count == 4);

    std::ostringstream headwing_output;
    write_eigenvalue_trace(
        headwing_output, 3, IterationChannel::Headwing, meanfield, kpoints);
    std::istringstream headwing_rows(headwing_output.str());
    int iteration = -1;
    int channel = -1;
    headwing_rows >> iteration >> channel;
    assert(headwing_rows);
    assert(iteration == 3);
    assert(channel == 2);
}

void test_invalid_trace_inputs_are_rejected()
{
    MeanField meanfield(1, 1, 1, 1, 1);
    std::ostringstream output;
    assert_throws([&] {
        write_eigenvalue_trace(
            output, -1, IterationChannel::Grid, meanfield,
            {{0.0, 0.0, 0.0}});
    });
    assert_throws([&] {
        write_eigenvalue_trace(
            output, 0, IterationChannel::Grid, meanfield, {});
    });
}

void test_matrix_components_include_full_complex_data_and_frequency_index()
{
    SpinKMatrixMap exchange;
    exchange[0][2] = librpa_int::Matz(2, 2, librpa_int::MAJOR::ROW);
    const double roundtrip_value = std::nextafter(1.0, 2.0);
    exchange[0][2](0, 0) = {roundtrip_value, 0.0};
    exchange[0][2](0, 1) = {2.0, 3.0};
    exchange[0][2](1, 0) = {2.0, -3.0};
    exchange[0][2](1, 1) = {4.0, 0.0};

    SpinKFrequencyMatrixMap sigma;
    sigma[0][2][0.25] = exchange[0][2].copy();
    sigma[0][2][1.50] = exchange[0][2].copy();

    std::ostringstream output;
    write_matrix_trace_header(output);
    write_matrix_component_trace(
        output, 1, IterationChannel::Grid, "exx", exchange);
    write_frequency_matrix_component_trace(
        output, 1, IterationChannel::Band, "sigma_c_iw", sigma);

    std::istringstream lines(output.str());
    std::string header;
    std::getline(lines, header);
    assert(header.find("frequency_index") != std::string::npos);
    assert(header.find("imag_value") != std::string::npos);

    int rows = 0;
    std::string row;
    while (std::getline(lines, row))
    {
        if (row.empty()) continue;
        std::istringstream values(row);
        int iteration = -1;
        int channel = -1;
        std::string component;
        int spin = -1;
        int kpoint = -1;
        int frequency_index = -2;
        double frequency = 0.0;
        int matrix_row = -1;
        int matrix_column = -1;
        double real = 0.0;
        double imaginary = 0.0;
        values >> iteration >> channel >> component >> spin >> kpoint >>
            frequency_index >> frequency >> matrix_row >> matrix_column >>
            real >> imaginary;
        assert(values);
        assert(iteration == 1);
        assert(spin == 0);
        assert(kpoint == 2);
        if (rows < 4)
        {
            assert(channel == 0);
            assert(component == "exx");
            assert(frequency_index == -1);
            assert(frequency == 0.0);
        }
        else
        {
            assert(channel == 1);
            assert(component == "sigma_c_iw");
            assert(frequency_index == (rows - 4) / 4);
            assert(std::abs(frequency -
                            (frequency_index == 0 ? 0.25 : 1.50)) < 1.0e-14);
        }
        if (matrix_row == 0 && matrix_column == 1)
        {
            assert(real == 2.0);
            assert(imaginary == 3.0);
        }
        if (matrix_row == 0 && matrix_column == 0)
        {
            assert(real == roundtrip_value);
        }
        ++rows;
    }
    assert(rows == 12);
}

void test_invalid_matrix_component_trace_is_rejected()
{
    SpinKMatrixMap valid;
    valid[0][0] = librpa_int::Matz(1, 1, librpa_int::MAJOR::ROW);
    valid[0][0](0, 0) = 1.0;
    std::ostringstream output;
    assert_throws([&] {
        write_matrix_component_trace(
            output, -1, IterationChannel::Grid, "exx", valid);
    });
    assert_throws([&] {
        write_matrix_component_trace(
            output, 1, IterationChannel::Grid, "bad component", valid);
    });
    assert_throws([&] {
        write_matrix_component_trace(
            output, 1, IterationChannel::Grid, "empty", {});
    });
    valid[0][0](0, 0) = {
        std::numeric_limits<double>::quiet_NaN(), 0.0};
    assert_throws([&] {
        write_matrix_component_trace(
            output, 1, IterationChannel::Grid, "nonfinite", valid);
    });
}

void test_occupation_and_scalar_components_are_explicit()
{
    MeanField meanfield(1, 2, 2, 1, 1);
    meanfield.get_weight()[0](0, 0) = 0.5;
    meanfield.get_weight()[0](0, 1) = 0.0;
    meanfield.get_weight()[0](1, 0) = 0.25;
    meanfield.get_weight()[0](1, 1) = 0.25;

    std::ostringstream output;
    write_occupation_trace(
        output, 4, IterationChannel::Grid, meanfield);
    write_scalar_component_trace(
        output, 4, IterationChannel::Grid, "fermi_energy_ha", -0.125);
    write_scalar_component_trace(
        output, 4, IterationChannel::Grid, "electron_count", 2.0);
    write_scalar_component_trace(
        output, 4, IterationChannel::Grid, "gap_ha", 0.25);

    int occupation_rows = 0;
    int scalar_rows = 0;
    std::istringstream lines(output.str());
    std::string row;
    while (std::getline(lines, row))
    {
        if (row.empty()) continue;
        std::istringstream values(row);
        int iteration = -1;
        int channel = -1;
        std::string component;
        int spin = -1;
        int kpoint = -1;
        int frequency_index = -2;
        double frequency = -1.0;
        int matrix_row = -1;
        int matrix_column = -1;
        double real = 0.0;
        double imaginary = 0.0;
        values >> iteration >> channel >> component >> spin >> kpoint >>
            frequency_index >> frequency >> matrix_row >> matrix_column >>
            real >> imaginary;
        assert(values);
        assert(iteration == 4);
        assert(channel == 0);
        assert(frequency_index == -1);
        assert(frequency == 0.0);
        assert(imaginary == 0.0);
        if (component == "occupation")
        {
            assert(spin == 0);
            assert(kpoint == occupation_rows / 2);
            assert(matrix_row == 0);
            assert(matrix_column == occupation_rows % 2);
            ++occupation_rows;
        }
        else
        {
            assert(component == "fermi_energy_ha" ||
                   component == "electron_count" || component == "gap_ha");
            assert(spin == 0);
            assert(kpoint == 0);
            assert(matrix_row == 0);
            assert(matrix_column == 0);
            ++scalar_rows;
        }
    }
    assert(occupation_rows == 4);
    assert(scalar_rows == 3);

    assert_throws([&] {
        write_scalar_component_trace(
            output, 4, IterationChannel::Grid, "bad",
            std::numeric_limits<double>::quiet_NaN());
    });
}

void test_wavefunction_and_velocity_components_are_explicit()
{
    MeanField meanfield(1, 1, 2, 3, 1);
    auto& wfc = meanfield.get_eigenvectors()[0][0][0];
    wfc.create(2, 3);
    for (int row = 0; row < wfc.nr; ++row)
    {
        for (int column = 0; column < wfc.nc; ++column)
        {
            wfc(row, column) = {
                static_cast<double>(10 * row + column),
                static_cast<double>(row - column)};
        }
    }

    std::vector<std::vector<std::vector<librpa_int::ComplexMatrix>>> velocity(1);
    velocity[0].resize(1);
    velocity[0][0].resize(3);
    for (int direction = 0; direction < 3; ++direction)
    {
        velocity[0][0][direction].create(2, 2);
        for (int row = 0; row < 2; ++row)
        {
            for (int column = 0; column < 2; ++column)
            {
                velocity[0][0][direction](row, column) = {
                    static_cast<double>(direction + row + column),
                    static_cast<double>(row - column)};
            }
        }
    }

    std::ostringstream output;
    write_wavefunction_trace(
        output, 2, IterationChannel::Grid, meanfield);
    write_velocity_trace(
        output, 2, IterationChannel::Grid, velocity);

    int wfc_rows = 0;
    int velocity_rows = 0;
    std::istringstream lines(output.str());
    std::string row;
    while (std::getline(lines, row))
    {
        if (row.empty()) continue;
        std::istringstream values(row);
        int iteration = -1;
        int channel = -1;
        std::string component;
        values >> iteration >> channel >> component;
        assert(values);
        assert(iteration == 2);
        assert(channel == 0);
        if (component == "wfc_spinor0")
        {
            ++wfc_rows;
        }
        else
        {
            assert(component == "velocity_x" ||
                   component == "velocity_y" ||
                   component == "velocity_z");
            ++velocity_rows;
        }
    }
    assert(wfc_rows == 6);
    assert(velocity_rows == 12);

    std::ostringstream provenance_output;
    write_wavefunction_trace(
        provenance_output, 0, IterationChannel::Headwing, meanfield,
        "headwing_reader_wfc");
    write_velocity_trace(
        provenance_output, 0, IterationChannel::Headwing, velocity,
        "headwing_reader_velocity");
    assert(provenance_output.str().find(
               "headwing_reader_wfc_spinor0") != std::string::npos);
    assert(provenance_output.str().find(
               "headwing_reader_velocity_x") != std::string::npos);
    assert(provenance_output.str().find(
               "headwing_reader_velocity_y") != std::string::npos);
    assert(provenance_output.str().find(
               "headwing_reader_velocity_z") != std::string::npos);

    velocity[0][0].pop_back();
    assert_throws([&] {
        write_velocity_trace(
            output, 2, IterationChannel::Grid, velocity);
    });
}

} // namespace

int main()
{
    test_summary_is_numeric_through_the_convergence_column();
    test_iteration_zero_and_band_channel_eigenvalues_are_explicit();
    test_invalid_trace_inputs_are_rejected();
    test_matrix_components_include_full_complex_data_and_frequency_index();
    test_invalid_matrix_component_trace_is_rejected();
    test_occupation_and_scalar_components_are_explicit();
    test_wavefunction_and_velocity_components_are_explicit();
    return 0;
}
