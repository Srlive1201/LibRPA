#include "correlation_potential.h"

#include "../core/analycont.h"

#include <cmath>
#include <stdexcept>
#include <string>

namespace librpa_int
{
namespace qsgw
{
namespace
{

bool finite_complex(const cplxdb value)
{
    return std::isfinite(value.real()) && std::isfinite(value.imag());
}

void validate_parameter_count(const int count,
                              const std::size_t data_size,
                              const char* label)
{
    if (count < -1 || count == 0 || (count == 1 && data_size > 1))
    {
        throw std::invalid_argument(
            std::string(label) +
            " must be -1 or at least 2 for a multi-point continuation");
    }
}

void validate_contract(
    const MeanField& meanfield,
    const std::vector<double>& source_frequencies,
    const std::map<double, Matz>& sigma,
    const int spin,
    const int kpoint,
    const CorrelationPotentialSettings& settings)
{
    if (!meanfield.initialized() || spin < 0 ||
        spin >= meanfield.get_n_spins() || kpoint < 0 ||
        kpoint >= meanfield.get_n_kpoints())
    {
        throw std::invalid_argument(
            "QSGW correlation potential received an invalid mean-field index");
    }
    if (!std::isfinite(meanfield.get_efermi()))
    {
        throw std::invalid_argument(
            "QSGW correlation potential Fermi energy is non-finite");
    }
    for (int band = 0; band < meanfield.get_n_bands(); ++band)
    {
        if (!std::isfinite(meanfield.get_eigenvals()[spin](kpoint, band)))
        {
            throw std::invalid_argument(
                "QSGW correlation potential eigenvalue is non-finite");
        }
    }

    if (source_frequencies.size() < 2 ||
        sigma.size() != source_frequencies.size())
    {
        throw std::invalid_argument(
            "QSGW correlation potential requires matching multi-point frequency data");
    }
    double previous = -1.0;
    for (const double frequency : source_frequencies)
    {
        if (!std::isfinite(frequency) || frequency <= 0.0 ||
            frequency <= previous)
        {
            throw std::invalid_argument(
                "QSGW imaginary frequencies must be finite, positive, and strictly increasing");
        }
        previous = frequency;
        const auto sigma_it = sigma.find(frequency);
        if (sigma_it == sigma.end())
        {
            throw std::invalid_argument(
                "QSGW self-energy map is missing a requested frequency");
        }
        const Matz& matrix = sigma_it->second;
        if (matrix.nr() != meanfield.get_n_bands() ||
            matrix.nc() != meanfield.get_n_bands())
        {
            throw std::invalid_argument(
                "QSGW self-energy matrix has an invalid shape");
        }
        for (int row = 0; row < matrix.nr(); ++row)
        {
            for (int column = 0; column < matrix.nc(); ++column)
            {
                if (!finite_complex(matrix(row, column)))
                {
                    throw std::invalid_argument(
                        "QSGW self-energy matrix contains non-finite data");
                }
            }
        }
    }

    validate_parameter_count(settings.n_params_anacon,
                             settings.resample_imagfreqs.empty()
                                 ? source_frequencies.size()
                                 : settings.resample_imagfreqs.size(),
                             "n_params_anacon");
    if (!settings.resample_imagfreqs.empty())
    {
        validate_parameter_count(settings.n_params_anacon_resample,
                                 source_frequencies.size(),
                                 "n_params_anacon_resample");
        double previous_imaginary = -1.0;
        for (const cplxdb frequency : settings.resample_imagfreqs)
        {
            if (!finite_complex(frequency) || frequency.real() != 0.0 ||
                frequency.imag() <= 0.0 ||
                frequency.imag() <= previous_imaginary)
            {
                throw std::invalid_argument(
                    "QSGW resampling points must be strictly increasing positive imaginary frequencies");
            }
            previous_imaginary = frequency.imag();
        }
    }

    if (settings.mode != CorrelationPotentialMode::ModeA &&
        settings.mode != CorrelationPotentialMode::ModeB)
    {
        throw std::invalid_argument(
            "QSGW correlation potential mode is invalid");
    }
}

cplxdb continue_element(
    const std::vector<cplxdb>& source_points,
    const std::vector<cplxdb>& source_data,
    const cplxdb target,
    const CorrelationPotentialSettings& settings)
{
    cplxdb result;
    if (settings.resample_imagfreqs.empty())
    {
        result = AnalyContPade(settings.n_params_anacon,
                              source_points, source_data)
                     .get(target);
    }
    else
    {
        const AnalyContPade source_pade(
            settings.n_params_anacon_resample,
            source_points, source_data);
        std::vector<cplxdb> resampled_data;
        resampled_data.reserve(settings.resample_imagfreqs.size());
        for (const cplxdb frequency : settings.resample_imagfreqs)
        {
            const cplxdb value = source_pade.get(frequency);
            if (!finite_complex(value))
            {
                throw std::runtime_error(
                    "QSGW first-stage analytic continuation produced non-finite data");
            }
            resampled_data.push_back(value);
        }
        result = AnalyContPade(settings.n_params_anacon,
                              settings.resample_imagfreqs,
                              resampled_data)
                     .get(target);
    }
    if (!finite_complex(result))
    {
        throw std::runtime_error(
            "QSGW analytic continuation produced non-finite data");
    }
    return result;
}

Matz evaluate_hermitian_sigma(
    const std::vector<double>& source_frequencies,
    const std::map<double, Matz>& sigma,
    const double target_energy,
    const CorrelationPotentialSettings& settings)
{
    std::vector<cplxdb> source_points;
    source_points.reserve(source_frequencies.size());
    for (const double frequency : source_frequencies)
    {
        source_points.emplace_back(0.0, frequency);
    }

    const int dimension = sigma.begin()->second.nr();
    Matz continued(dimension, dimension, MAJOR::ROW);
    for (int row = 0; row < dimension; ++row)
    {
        for (int column = 0; column < dimension; ++column)
        {
            std::vector<cplxdb> source_data;
            source_data.reserve(source_frequencies.size());
            for (const double frequency : source_frequencies)
            {
                source_data.push_back(sigma.at(frequency)(row, column));
            }
            continued(row, column) = continue_element(
                source_points, source_data, cplxdb(target_energy, 0.0),
                settings);
        }
    }

    Matz hermitian(dimension, dimension, MAJOR::ROW);
    for (int row = 0; row < dimension; ++row)
    {
        for (int column = 0; column < dimension; ++column)
        {
            hermitian(row, column) =
                0.5 * (continued(row, column) +
                       std::conj(continued(column, row)));
        }
    }
    return hermitian;
}

} // namespace

Matz build_qsgw_correlation_potential(
    const MeanField& meanfield,
    const std::vector<double>& source_frequencies,
    const std::map<double, Matz>& sigma_imaginary_axis,
    const int spin,
    const int kpoint,
    const CorrelationPotentialSettings& settings)
{
    validate_contract(meanfield, source_frequencies, sigma_imaginary_axis,
                      spin, kpoint, settings);
    const int dimension = meanfield.get_n_bands();
    const double fermi_energy = meanfield.get_efermi();

    std::vector<Matz> sigma_at_state_energy;
    sigma_at_state_energy.reserve(static_cast<std::size_t>(dimension));
    for (int band = 0; band < dimension; ++band)
    {
        sigma_at_state_energy.push_back(evaluate_hermitian_sigma(
            source_frequencies, sigma_imaginary_axis,
            meanfield.get_eigenvals()[spin](kpoint, band) - fermi_energy,
            settings));
    }

    Matz result(dimension, dimension, MAJOR::ROW);
    if (settings.mode == CorrelationPotentialMode::ModeB)
    {
        const Matz sigma_at_fermi = evaluate_hermitian_sigma(
            source_frequencies, sigma_imaginary_axis, 0.0, settings);
        for (int row = 0; row < dimension; ++row)
        {
            for (int column = 0; column < dimension; ++column)
            {
                result(row, column) =
                    row == column
                        ? sigma_at_state_energy[row](row, row)
                        : sigma_at_fermi(row, column);
            }
        }
        return result;
    }

    for (int row = 0; row < dimension; ++row)
    {
        for (int column = 0; column < dimension; ++column)
        {
            result(row, column) =
                0.5 * (sigma_at_state_energy[row](row, column) +
                       sigma_at_state_energy[column](row, column));
        }
    }
    return result;
}

} // namespace qsgw
} // namespace librpa_int
