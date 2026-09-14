#include "iteration_trace.h"

#include "../../src/utils/constants.h"

#include <cctype>
#include <cmath>
#include <iomanip>
#include <limits>
#include <ostream>
#include <stdexcept>

namespace librpa_int
{
namespace qsgw
{
namespace
{

constexpr int trace_precision = std::numeric_limits<double>::max_digits10;

void require_finite_nonnegative(const double value, const char* label)
{
    if (!std::isfinite(value) || value < 0.0)
    {
        throw std::invalid_argument(
            std::string("QSGW iteration trace ") + label +
            " must be finite and nonnegative");
    }
}

bool finite_complex(const cplxdb value)
{
    return std::isfinite(value.real()) && std::isfinite(value.imag());
}

void validate_matrix_component_label(const std::string& component)
{
    if (component.empty())
    {
        throw std::invalid_argument(
            "QSGW matrix trace component must not be empty");
    }
    for (const char character : component)
    {
        const unsigned char byte = static_cast<unsigned char>(character);
        if (!std::isalnum(byte) && character != '_' &&
            character != '-' && character != '.')
        {
            throw std::invalid_argument(
                "QSGW matrix trace component is not a machine token");
        }
    }
}

void write_matrix_rows(std::ostream& output,
                       const int iteration,
                       const IterationChannel channel,
                       const std::string& component,
                       const int spin,
                       const int kpoint,
                       const int frequency_index,
                       const double frequency,
                       const Matz& matrix)
{
    if (spin < 0 || kpoint < 0 || frequency_index < -1 ||
        !std::isfinite(frequency) || frequency < 0.0 ||
        matrix.nr() <= 0 || matrix.nc() <= 0)
    {
        throw std::invalid_argument(
            "QSGW matrix trace layout is invalid");
    }
    for (int row = 0; row < matrix.nr(); ++row)
    {
        for (int column = 0; column < matrix.nc(); ++column)
        {
            const cplxdb value = matrix(row, column);
            if (!finite_complex(value))
            {
                throw std::invalid_argument(
                    "QSGW matrix trace contains non-finite data");
            }
            output << iteration << " "
                   << static_cast<int>(channel) << " "
                   << component << " " << spin << " " << kpoint << " "
                   << frequency_index << " " << frequency << " "
                   << row << " " << column << " "
                   << value.real() << " " << value.imag() << "\n";
        }
    }
}

class StreamFormatGuard
{
public:
    explicit StreamFormatGuard(std::ostream& output)
        : output_(output), flags_(output.flags()), precision_(output.precision())
    {
    }

    ~StreamFormatGuard()
    {
        output_.flags(flags_);
        output_.precision(precision_);
    }

private:
    std::ostream& output_;
    std::ios::fmtflags flags_;
    std::streamsize precision_;
};

} // namespace

void write_iteration_summary_header(std::ostream& output)
{
    output << "# iter max_delta_eV residual_l2_Ha residual_max_Ha "
              "efermi_eV gap_eV electron_count beta converged\n";
}

void write_iteration_summary(std::ostream& output,
                             const IterationSummary& summary)
{
    if (summary.iteration < 0)
    {
        throw std::invalid_argument(
            "QSGW iteration trace index must be nonnegative");
    }
    require_finite_nonnegative(
        summary.maximum_eigenvalue_change_ev, "maximum eigenvalue change");
    require_finite_nonnegative(summary.residual_l2_ha, "residual L2 norm");
    require_finite_nonnegative(summary.residual_max_ha, "residual max norm");
    require_finite_nonnegative(summary.gap_ev, "gap");
    require_finite_nonnegative(summary.electron_count, "electron count");
    if (!std::isfinite(summary.fermi_energy_ev) ||
        !(summary.beta > 0.0 && summary.beta <= 1.0))
        throw std::invalid_argument(
            "QSGW iteration trace contains invalid scalar data");

    StreamFormatGuard guard(output);
    output << std::scientific << std::setprecision(trace_precision)
           << summary.iteration << " "
           << summary.maximum_eigenvalue_change_ev << " "
           << summary.residual_l2_ha << " "
           << summary.residual_max_ha << " "
           << summary.fermi_energy_ev << " "
           << summary.gap_ev << " "
           << summary.electron_count << " "
           << summary.beta << " "
           << (summary.converged ? 1 : 0) << "\n";
}

void write_eigenvalue_trace_header(std::ostream& output)
{
    output << "# iter channel spin kpoint kx ky kz band energy_eV\n";
}

void write_eigenvalue_trace(
    std::ostream& output,
    const int iteration,
    const IterationChannel channel,
    const MeanField& meanfield,
    const std::vector<Vector3_Order<double>>& kpoints)
{
    if (iteration < 0 || !meanfield.initialized() ||
        kpoints.size() !=
            static_cast<std::size_t>(meanfield.get_n_kpoints()))
    {
        throw std::invalid_argument(
            "QSGW eigenvalue trace layout is invalid");
    }

    StreamFormatGuard guard(output);
    output << std::scientific << std::setprecision(trace_precision);
    for (int spin = 0; spin < meanfield.get_n_spins(); ++spin)
    {
        for (int kpoint = 0; kpoint < meanfield.get_n_kpoints(); ++kpoint)
        {
            const auto& coordinate =
                kpoints[static_cast<std::size_t>(kpoint)];
            if (!std::isfinite(coordinate.x) ||
                !std::isfinite(coordinate.y) ||
                !std::isfinite(coordinate.z))
            {
                throw std::invalid_argument(
                    "QSGW eigenvalue trace k point is non-finite");
            }
            for (int band = 0; band < meanfield.get_n_bands(); ++band)
            {
                const double eigenvalue =
                    meanfield.get_eigenvals()[spin](kpoint, band);
                if (!std::isfinite(eigenvalue))
                {
                    throw std::invalid_argument(
                        "QSGW eigenvalue trace value is non-finite");
                }
                output << iteration << " "
                       << static_cast<int>(channel) << " "
                       << spin << " " << kpoint << " "
                       << coordinate.x << " " << coordinate.y << " "
                       << coordinate.z << " " << band << " "
                       << eigenvalue * HA2EV << "\n";
            }
        }
    }
}

void write_matrix_trace_header(std::ostream& output)
{
    output << "# iter channel component spin kpoint frequency_index "
              "frequency_Ha row column real_value imag_value\n";
}

void write_matrix_component_trace(
    std::ostream& output,
    const int iteration,
    const IterationChannel channel,
    const std::string& component,
    const SpinKMatrixMap& matrices)
{
    if (iteration < 0 || matrices.empty())
    {
        throw std::invalid_argument(
            "QSGW matrix component trace input is empty or invalid");
    }
    validate_matrix_component_label(component);
    StreamFormatGuard guard(output);
    output << std::scientific << std::setprecision(trace_precision);
    for (const auto& [spin, by_kpoint] : matrices)
    {
        if (by_kpoint.empty())
        {
            throw std::invalid_argument(
                "QSGW matrix component trace has an empty spin channel");
        }
        for (const auto& [kpoint, matrix] : by_kpoint)
        {
            write_matrix_rows(
                output, iteration, channel, component, spin, kpoint,
                -1, 0.0, matrix);
        }
    }
}

void write_frequency_matrix_component_trace(
    std::ostream& output,
    const int iteration,
    const IterationChannel channel,
    const std::string& component,
    const SpinKFrequencyMatrixMap& matrices)
{
    if (iteration < 0 || matrices.empty())
    {
        throw std::invalid_argument(
            "QSGW frequency matrix trace input is empty or invalid");
    }
    validate_matrix_component_label(component);
    StreamFormatGuard guard(output);
    output << std::scientific << std::setprecision(trace_precision);
    for (const auto& [spin, by_kpoint] : matrices)
    {
        if (by_kpoint.empty())
        {
            throw std::invalid_argument(
                "QSGW frequency matrix trace has an empty spin channel");
        }
        for (const auto& [kpoint, by_frequency] : by_kpoint)
        {
            if (by_frequency.empty())
            {
                throw std::invalid_argument(
                    "QSGW frequency matrix trace has an empty k point");
            }
            int frequency_index = 0;
            for (const auto& [frequency, matrix] : by_frequency)
            {
                write_matrix_rows(
                    output, iteration, channel, component, spin, kpoint,
                    frequency_index, frequency, matrix);
                ++frequency_index;
            }
        }
    }
}

void write_occupation_trace(
    std::ostream& output,
    const int iteration,
    const IterationChannel channel,
    const MeanField& meanfield)
{
    if (!meanfield.initialized())
    {
        throw std::invalid_argument(
            "QSGW occupation trace mean field is not initialized");
    }
    SpinKMatrixMap matrices;
    for (int spin = 0; spin < meanfield.get_n_spins(); ++spin)
    {
        for (int kpoint = 0; kpoint < meanfield.get_n_kpoints(); ++kpoint)
        {
            Matz occupations(1, meanfield.get_n_bands(), MAJOR::ROW);
            for (int band = 0; band < meanfield.get_n_bands(); ++band)
            {
                occupations(0, band) =
                    meanfield.get_weight()[spin](kpoint, band);
            }
            matrices[spin][kpoint] = std::move(occupations);
        }
    }
    write_matrix_component_trace(
        output, iteration, channel, "occupation", matrices);
}

void write_wavefunction_trace(
    std::ostream& output,
    const int iteration,
    const IterationChannel channel,
    const MeanField& meanfield,
    const std::string& component_prefix)
{
    if (component_prefix.empty())
    {
        throw std::invalid_argument(
            "QSGW wavefunction trace component prefix is empty");
    }
    if (!meanfield.initialized())
    {
        throw std::invalid_argument(
            "QSGW wavefunction trace mean field is not initialized");
    }
    for (int spinor = 0; spinor < meanfield.get_n_spinor(); ++spinor)
    {
        SpinKMatrixMap matrices;
        for (int spin = 0; spin < meanfield.get_n_spins(); ++spin)
        {
            for (int kpoint = 0; kpoint < meanfield.get_n_kpoints(); ++kpoint)
            {
                const ComplexMatrix* block =
                    meanfield.find_wfc(spin, spinor, kpoint);
                if (block == nullptr ||
                    block->nr != meanfield.get_n_bands() ||
                    block->nc != meanfield.get_n_aos())
                {
                    throw std::invalid_argument(
                        "QSGW wavefunction trace map is incomplete or has an invalid shape");
                }
                Matz matrix(block->nr, block->nc, MAJOR::ROW);
                for (int row = 0; row < block->nr; ++row)
                {
                    for (int column = 0; column < block->nc; ++column)
                    {
                        matrix(row, column) = (*block)(row, column);
                    }
                }
                matrices[spin][kpoint] = std::move(matrix);
            }
        }
        write_matrix_component_trace(
            output, iteration, channel,
            component_prefix + "_spinor" + std::to_string(spinor), matrices);
    }
}

} // namespace qsgw
} // namespace librpa_int
