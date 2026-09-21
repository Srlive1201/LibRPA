#include "projection_target.h"

#include <cmath>
#include <stdexcept>

namespace librpa_int
{
namespace qsgw
{

ProjectionTargetShape validate_projection_target(
    const MeanField& reference,
    const std::vector<Vector3_Order<double>>& kpoints,
    const int expected_n_spins,
    const int expected_n_spinors,
    const int expected_n_aos,
    const std::string& label)
{
    if (!reference.initialized())
    {
        throw std::invalid_argument(
            "QSGW " + label + " projection mean field is not initialized");
    }
    if (reference.get_n_spins() != expected_n_spins ||
        reference.get_n_spinor() != expected_n_spinors ||
        reference.get_n_aos() != expected_n_aos)
    {
        throw std::invalid_argument(
            "QSGW " + label +
            " projection target is incompatible with the source basis");
    }
    if (reference.get_n_kpoints() != static_cast<int>(kpoints.size()))
    {
        throw std::invalid_argument(
            "QSGW " + label +
            " projection k-point coordinates do not match the mean field");
    }

    for (const auto& kpoint : kpoints)
    {
        if (!std::isfinite(kpoint.x) || !std::isfinite(kpoint.y) ||
            !std::isfinite(kpoint.z))
        {
            throw std::invalid_argument(
                "QSGW " + label +
                " projection contains a non-finite k-point");
        }
    }

    for (int spin = 0; spin < reference.get_n_spins(); ++spin)
    {
        for (int spinor = 0; spinor < reference.get_n_spinor(); ++spinor)
        {
            for (int kpoint = 0; kpoint < reference.get_n_kpoints(); ++kpoint)
            {
                const ComplexMatrix* wavefunction =
                    reference.find_wfc(spin, spinor, kpoint);
                if (wavefunction == nullptr ||
                    wavefunction->nr != reference.get_n_bands() ||
                    wavefunction->nc != reference.get_n_aos())
                {
                    throw std::invalid_argument(
                        "QSGW " + label +
                        " projection wavefunction map is incomplete or malformed");
                }
                for (int index = 0; index < wavefunction->size; ++index)
                {
                    const cplxdb value = wavefunction->c[index];
                    if (!std::isfinite(value.real()) ||
                        !std::isfinite(value.imag()))
                    {
                        throw std::invalid_argument(
                            "QSGW " + label +
                            " projection wavefunction contains non-finite data");
                    }
                }
            }
        }
    }

    return {
        reference.get_n_spins(), reference.get_n_spinor(),
        reference.get_n_kpoints(), reference.get_n_bands(),
        reference.get_n_aos()};
}

} // namespace qsgw
} // namespace librpa_int
