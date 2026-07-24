#include "hamiltonian_cut.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>

namespace librpa_int
{
namespace qsgw
{
namespace
{

void validate_options(const HamiltonianCutOptions& options)
{
    if (options.unoccupied_keep < 0)
        throw std::invalid_argument(
            "QSGW Hamiltonian cut requires non-negative unoccupied_keep");
    if (!std::isfinite(options.shift_ha))
        throw std::invalid_argument(
            "QSGW Hamiltonian cut shift must be finite");
    switch (options.mode)
    {
        case HamiltonianCutMode::Uncut:
        case HamiltonianCutMode::ReferenceDiagonal:
        case HamiltonianCutMode::ShiftedReferenceDiagonal:
            return;
    }
    throw std::invalid_argument("Unsupported QSGW Hamiltonian cut mode");
}

void validate_matrix(const Matz& matrix, const int dimension,
                     const std::string& label)
{
    if (matrix.nr() != dimension || matrix.nc() != dimension)
        throw std::invalid_argument(
            "QSGW Hamiltonian cut " + label + " dimension differs");
    for (int row = 0; row < dimension; ++row)
    {
        for (int column = 0; column < dimension; ++column)
        {
            const cplxdb value = matrix(row, column);
            if (!std::isfinite(value.real()) || !std::isfinite(value.imag()))
                throw std::invalid_argument(
                    "QSGW Hamiltonian cut " + label +
                    " contains non-finite data");
        }
    }
}

} // namespace

HamiltonianCutMode hamiltonian_cut_mode_from_int(const int mode)
{
    switch (mode)
    {
        case 0: return HamiltonianCutMode::Uncut;
        case 1: return HamiltonianCutMode::ReferenceDiagonal;
        case 2: return HamiltonianCutMode::ShiftedReferenceDiagonal;
        default:
            throw std::invalid_argument(
                "QSGW Hamiltonian cut mode must be 0, 1, or 2");
    }
}

SpinKMatrixMap apply_hamiltonian_cut(
    const SpinKMatrixMap& hamiltonian,
    const SpinKMatrixMap& reference_hamiltonian,
    const MeanField& live_meanfield,
    const HamiltonianCutOptions& options)
{
    validate_options(options);
    if (!live_meanfield.initialized())
        throw std::invalid_argument(
            "QSGW Hamiltonian cut mean field is not initialized");
    if (options.mode != HamiltonianCutMode::Uncut &&
        !std::isfinite(live_meanfield.get_efermi()))
        throw std::invalid_argument(
            "QSGW Hamiltonian cut Fermi energy is non-finite");
    if (hamiltonian.size() != reference_hamiltonian.size() ||
        static_cast<int>(hamiltonian.size()) != live_meanfield.get_n_spins())
        throw std::invalid_argument(
            "QSGW Hamiltonian cut spin layout differs");

    const int dimension = live_meanfield.get_n_bands();
    SpinKMatrixMap result;
    for (const auto& [spin, by_kpoint] : hamiltonian)
    {
        const auto reference_spin = reference_hamiltonian.find(spin);
        if (spin < 0 || spin >= live_meanfield.get_n_spins() ||
            reference_spin == reference_hamiltonian.end() ||
            by_kpoint.size() != reference_spin->second.size() ||
            static_cast<int>(by_kpoint.size()) !=
                live_meanfield.get_n_kpoints())
            throw std::invalid_argument(
                "QSGW Hamiltonian cut k-point layout differs");

        for (const auto& [kpoint, matrix] : by_kpoint)
        {
            const auto reference_kpoint =
                reference_spin->second.find(kpoint);
            if (kpoint < 0 || kpoint >= live_meanfield.get_n_kpoints() ||
                reference_kpoint == reference_spin->second.end())
                throw std::invalid_argument(
                    "QSGW Hamiltonian cut k-point keys differ");
            validate_matrix(matrix, dimension, "input");
            validate_matrix(reference_kpoint->second, dimension,
                            "reference");

            Matz cut = matrix.copy();
            if (options.mode != HamiltonianCutMode::Uncut)
            {
                int occupied = 0;
                for (int band = 0; band < dimension; ++band)
                {
                    const double energy =
                        live_meanfield.get_eigenvals()[spin](kpoint, band);
                    if (!std::isfinite(energy))
                        throw std::invalid_argument(
                            "QSGW Hamiltonian cut eigenvalue is non-finite");
                    if (energy < live_meanfield.get_efermi()) ++occupied;
                }
                const int active_limit = std::min(
                    dimension, occupied + options.unoccupied_keep);
                const double shift =
                    options.mode ==
                            HamiltonianCutMode::ShiftedReferenceDiagonal
                        ? options.shift_ha
                        : 0.0;
                for (int row = 0; row < dimension; ++row)
                {
                    for (int column = 0; column < dimension; ++column)
                    {
                        if (row < active_limit && column < active_limit)
                            continue;
                        cut(row, column) = row == column
                            ? reference_kpoint->second(row, row) + shift
                            : cplxdb(0.0, 0.0);
                    }
                }
            }
            result[spin][kpoint] = std::move(cut);
        }
    }
    return result;
}

} // namespace qsgw
} // namespace librpa_int
