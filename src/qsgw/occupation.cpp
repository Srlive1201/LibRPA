#include "occupation.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace librpa_int
{
namespace qsgw
{

namespace
{

struct State
{
    double energy;
    double capacity;
    int spin;
    int kpoint;
    int band;
};

void validate_same_layout(const MeanField& live, const MeanField& reference)
{
    if (!live.initialized() || !reference.initialized())
    {
        throw std::invalid_argument("QSGW occupation mean fields must be initialized");
    }
    if (live.get_n_spins() != reference.get_n_spins() ||
        live.get_n_kpoints() != reference.get_n_kpoints() ||
        live.get_n_bands() != reference.get_n_bands() ||
        live.get_n_spinor() != reference.get_n_spinor())
    {
        throw std::invalid_argument("QSGW live and reference mean-field layouts differ");
    }
}

double state_capacity_factor(const MeanField& meanfield)
{
    return 2.0 /
           static_cast<double>(meanfield.get_n_spins() * meanfield.get_n_spinor());
}

void validate_kpoint_weights(const MeanField& meanfield,
                             const std::vector<double>& kpoint_weights,
                             const double tolerance)
{
    if (!meanfield.initialized())
    {
        throw std::invalid_argument("QSGW occupation mean field must be initialized");
    }
    if (!std::isfinite(tolerance) || !(tolerance > 0.0))
    {
        throw std::invalid_argument("QSGW occupation tolerance must be finite and positive");
    }
    if (kpoint_weights.size() !=
        static_cast<std::size_t>(meanfield.get_n_kpoints()))
    {
        throw std::invalid_argument("QSGW k-point weight count does not match mean field");
    }

    double weight_sum = 0.0;
    for (const double weight: kpoint_weights)
    {
        if (!std::isfinite(weight) || weight < 0.0)
        {
            throw std::invalid_argument(
                "QSGW k-point weights must be finite and nonnegative");
        }
        weight_sum += weight;
    }
    if (!std::isfinite(weight_sum) ||
        std::abs(weight_sum - 1.0) > tolerance)
    {
        throw std::invalid_argument("QSGW k-point weights must sum to one");
    }
}

} // namespace

double physical_electron_count(
    const MeanField& meanfield,
    const std::vector<double>& kpoint_weights,
    const double tolerance)
{
    validate_kpoint_weights(meanfield, kpoint_weights, tolerance);
    const double capacity_factor = state_capacity_factor(meanfield);

    double result = 0.0;
    for (int spin = 0; spin < meanfield.get_n_spins(); ++spin)
    {
        for (int kpoint = 0; kpoint < meanfield.get_n_kpoints(); ++kpoint)
        {
            const double state_capacity =
                capacity_factor * kpoint_weights[kpoint];
            for (int band = 0; band < meanfield.get_n_bands(); ++band)
            {
                const double stored_weight =
                    meanfield.get_weight()[spin](kpoint, band);
                if (!std::isfinite(stored_weight) ||
                    stored_weight < -tolerance ||
                    stored_weight > state_capacity + tolerance)
                {
                    throw std::invalid_argument(
                        "QSGW occupation exceeds the MeanField k-weighted state capacity");
                }
                result += stored_weight;
            }
        }
    }
    if (!std::isfinite(result))
    {
        throw std::invalid_argument("QSGW physical electron count is not finite");
    }
    return result;
}

OccupationResult analyze_qsgw_occupations(
    const MeanField& meanfield,
    const std::vector<double>& kpoint_weights,
    const double total_electrons,
    const OccupationSettings& settings)
{
    MeanField probe = meanfield;
    return update_qsgw_occupations(
        probe, meanfield, kpoint_weights, total_electrons, settings);
}

OccupationResult update_qsgw_occupations(
    MeanField& live_meanfield,
    const MeanField& reference_meanfield,
    const std::vector<double>& kpoint_weights,
    const double total_electrons,
    const OccupationSettings& settings)
{
    validate_same_layout(live_meanfield, reference_meanfield);
    if (&live_meanfield == &reference_meanfield)
    {
        throw std::invalid_argument(
            "QSGW live and immutable reference mean fields must not alias");
    }
    if (settings.temperature_kelvin != 0.0)
    {
        throw std::invalid_argument("Finite-temperature QSGW occupations are not implemented");
    }
    if (!std::isfinite(settings.degeneracy_tolerance_ha) ||
        !std::isfinite(settings.electron_tolerance) ||
        !(settings.degeneracy_tolerance_ha > 0.0) ||
        !(settings.electron_tolerance > 0.0))
    {
        throw std::invalid_argument("QSGW occupation tolerances must be finite and positive");
    }
    validate_kpoint_weights(
        live_meanfield, kpoint_weights, settings.electron_tolerance);

    const double capacity_factor = state_capacity_factor(reference_meanfield);
    const double reference_electrons = physical_electron_count(
        reference_meanfield, kpoint_weights, settings.electron_tolerance);
    double total_capacity = 0.0;
    std::vector<State> states;
    states.reserve(static_cast<std::size_t>(live_meanfield.get_n_spins()) *
                   live_meanfield.get_n_kpoints() * live_meanfield.get_n_bands());

    for (int spin = 0; spin < live_meanfield.get_n_spins(); ++spin)
    {
        for (int kpoint = 0; kpoint < live_meanfield.get_n_kpoints(); ++kpoint)
        {
            const double state_capacity =
                capacity_factor * kpoint_weights[kpoint];
            for (int band = 0; band < live_meanfield.get_n_bands(); ++band)
            {
                const double energy =
                    live_meanfield.get_eigenvals()[spin](kpoint, band);
                if (!std::isfinite(energy))
                {
                    throw std::invalid_argument(
                        "QSGW eigenvalues must be finite before occupation filling");
                }
                total_capacity += state_capacity;
                states.push_back(
                    {energy, state_capacity, spin, kpoint, band});
            }
        }
    }

    if (!std::isfinite(total_electrons) || total_electrons < 0.0 ||
        total_electrons > total_capacity + settings.electron_tolerance)
    {
        throw std::invalid_argument("QSGW target electron count is outside state capacity");
    }
    if (!std::isfinite(reference_electrons) || !std::isfinite(total_capacity))
    {
        throw std::invalid_argument(
            "QSGW reference electron count or state capacity is not finite");
    }
    if (std::abs(reference_electrons - total_electrons) > settings.electron_tolerance)
    {
        throw std::invalid_argument(
            "QSGW target electron count differs from immutable reference occupations");
    }

    std::sort(states.begin(), states.end(), [](const State& lhs, const State& rhs) {
        if (lhs.energy != rhs.energy)
        {
            return lhs.energy < rhs.energy;
        }
        if (lhs.spin != rhs.spin)
        {
            return lhs.spin < rhs.spin;
        }
        if (lhs.kpoint != rhs.kpoint)
        {
            return lhs.kpoint < rhs.kpoint;
        }
        return lhs.band < rhs.band;
    });

    auto updated_weights = live_meanfield.get_weight();
    for (auto& spin_weights: updated_weights)
    {
        spin_weights.zero_out();
    }
    double remaining = total_electrons;
    for (std::size_t begin = 0; begin < states.size();)
    {
        std::size_t end = begin + 1;
        while (end < states.size() &&
               std::abs(states[end].energy - states[begin].energy) <=
                   settings.degeneracy_tolerance_ha)
        {
            ++end;
        }

        double group_capacity = 0.0;
        for (std::size_t index = begin; index < end; ++index)
        {
            group_capacity += states[index].capacity;
        }
        const double fraction = group_capacity > 0.0
                                    ? std::clamp(remaining / group_capacity, 0.0, 1.0)
                                    : 0.0;
        for (std::size_t index = begin; index < end; ++index)
        {
            const State& state = states[index];
            updated_weights[state.spin](state.kpoint, state.band) =
                fraction * state.capacity;
        }
        remaining -= fraction * group_capacity;
        if (remaining < settings.electron_tolerance)
        {
            remaining = 0.0;
        }
        begin = end;
    }
    if (remaining > settings.electron_tolerance)
    {
        throw std::runtime_error("QSGW global occupation filling did not conserve charge");
    }

    OccupationResult result;
    double filling_vbm = -std::numeric_limits<double>::infinity();
    double filling_cbm = std::numeric_limits<double>::infinity();
    double lowest_active_energy = std::numeric_limits<double>::infinity();
    double highest_active_energy = -std::numeric_limits<double>::infinity();
    for (const State& state: states)
    {
        const double occupation =
            updated_weights[state.spin](state.kpoint, state.band);
        result.electron_count += occupation;
        if (state.capacity <= settings.electron_tolerance)
        {
            continue;
        }
        lowest_active_energy = std::min(lowest_active_energy, state.energy);
        highest_active_energy = std::max(highest_active_energy, state.energy);
        if (occupation > settings.electron_tolerance)
        {
            filling_vbm = std::max(filling_vbm, state.energy);
        }
        if (occupation < state.capacity - settings.electron_tolerance)
        {
            filling_cbm = std::min(filling_cbm, state.energy);
        }
        if (occupation > settings.electron_tolerance &&
            occupation < state.capacity - settings.electron_tolerance)
        {
            result.metallic = true;
        }
    }

    if (!std::isfinite(filling_vbm))
    {
        filling_vbm = std::isfinite(lowest_active_energy)
                           ? lowest_active_energy
                           : states.front().energy;
    }
    if (!std::isfinite(filling_cbm))
    {
        filling_cbm = std::isfinite(highest_active_energy)
                           ? highest_active_energy
                           : states.back().energy;
    }
    result.chemical_potential =
        0.5 * (filling_vbm + filling_cbm);

    result.vbm = -std::numeric_limits<double>::infinity();
    result.cbm = std::numeric_limits<double>::infinity();
    for (const State& state: states)
    {
        const double reference_occupation =
            reference_meanfield.get_weight()[state.spin](
                state.kpoint, state.band);
        if (reference_occupation > settings.electron_tolerance)
        {
            result.vbm = std::max(result.vbm, state.energy);
        }
        if (reference_occupation <
            state.capacity - settings.electron_tolerance)
        {
            result.cbm = std::min(result.cbm, state.energy);
        }
        if (reference_occupation > settings.electron_tolerance &&
            reference_occupation <
                state.capacity - settings.electron_tolerance)
        {
            result.metallic = true;
        }
    }
    if (!std::isfinite(result.vbm)) result.vbm = filling_vbm;
    if (!std::isfinite(result.cbm)) result.cbm = filling_cbm;
    const double manifold_gap = result.cbm - result.vbm;
    if (manifold_gap <= settings.degeneracy_tolerance_ha)
    {
        result.metallic = true;
    }
    result.gap = result.metallic ? 0.0 : manifold_gap;

    if (std::abs(result.electron_count - total_electrons) >
        settings.electron_tolerance)
    {
        throw std::runtime_error("QSGW updated occupations do not conserve charge");
    }
    live_meanfield.get_weight() = std::move(updated_weights);
    live_meanfield.get_efermi() = result.chemical_potential;
    return result;
}

} // namespace qsgw
} // namespace librpa_int
