#include "band_output.h"

#include "../utils/constants.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <ostream>
#include <stdexcept>
#include <string>

namespace librpa_int
{
namespace qsgw
{
namespace
{

void validate_layout(
    const MeanField& live,
    const MeanField& reference,
    const std::vector<Vector3_Order<double>>& kpoints,
    const int spin,
    const double chemical_potential)
{
    if (!live.initialized() || !reference.initialized())
        throw std::invalid_argument(
            "QSGW band output mean fields must be initialized");
    if (live.get_n_spins() != reference.get_n_spins() ||
        live.get_n_kpoints() != reference.get_n_kpoints() ||
        live.get_n_bands() != reference.get_n_bands() ||
        live.get_n_spinor() != reference.get_n_spinor())
        throw std::invalid_argument(
            "QSGW band output mean-field layouts differ");
    if (spin < 0 || spin >= live.get_n_spins())
        throw std::invalid_argument(
            "QSGW band output spin index is out of range");
    if (kpoints.size() !=
        static_cast<std::size_t>(live.get_n_kpoints()))
        throw std::invalid_argument(
            "QSGW band output k-point count differs");
    if (!std::isfinite(chemical_potential))
        throw std::invalid_argument(
            "QSGW band output chemical potential is non-finite");
}

const Matz& matrix_at(const SpinKMatrixMap& matrices,
                      const int spin,
                      const int kpoint,
                      const int dimension,
                      const char* label)
{
    const auto spin_it = matrices.find(spin);
    if (spin_it == matrices.end())
        throw std::invalid_argument(
            std::string("QSGW band output is missing ") + label +
            " spin data");
    const auto kpoint_it = spin_it->second.find(kpoint);
    if (kpoint_it == spin_it->second.end() ||
        kpoint_it->second.nr() != dimension ||
        kpoint_it->second.nc() != dimension)
        throw std::invalid_argument(
            std::string("QSGW band output has malformed ") + label +
            " k-point data");
    return kpoint_it->second;
}

double checked_real(const cplxdb value, const char* label)
{
    if (!std::isfinite(value.real()) || !std::isfinite(value.imag()))
        throw std::invalid_argument(
            std::string("QSGW band output has invalid ") + label +
            " diagonal data");
    return value.real();
}

void write_kpoint_prefix(std::ostream& output,
                         const int kpoint,
                         const Vector3_Order<double>& coordinate)
{
    output << std::setw(5) << kpoint + 1
           << std::setw(15) << std::setprecision(7) << coordinate.x
           << std::setw(15) << std::setprecision(7) << coordinate.y
           << std::setw(15) << std::setprecision(7) << coordinate.z;
}

} // namespace

void write_qsgw_band_spin_tables(
    std::ostream& ks_output,
    std::ostream& exx_output,
    std::ostream& qsgw_output,
    const MeanField& live_band,
    const MeanField& reference_band,
    const std::vector<Vector3_Order<double>>& kpoints,
    const SpinKMatrixMap& reference_hamiltonian,
    const SpinKMatrixMap& dft_vxc,
    const SpinKMatrixMap& exchange,
    const int spin,
    const double chemical_potential_ha)
{
    validate_layout(
        live_band, reference_band, kpoints, spin,
        chemical_potential_ha);
    if (!ks_output.good() || !exx_output.good() || !qsgw_output.good())
        throw std::invalid_argument(
            "QSGW band output stream is not writable");

    ks_output << std::fixed;
    exx_output << std::fixed;
    qsgw_output << std::fixed;
    const double occupation_scale = static_cast<double>(
        live_band.get_n_kpoints() * live_band.get_n_spins());
    for (int kpoint = 0; kpoint < live_band.get_n_kpoints(); ++kpoint)
    {
        const Matz& h_ks = matrix_at(
            reference_hamiltonian, spin, kpoint,
            live_band.get_n_bands(), "KS Hamiltonian");
        const Matz& vxc = matrix_at(
            dft_vxc, spin, kpoint, live_band.get_n_bands(), "Vxc");
        const Matz& exx = matrix_at(
            exchange, spin, kpoint, live_band.get_n_bands(), "EXX");
        write_kpoint_prefix(ks_output, kpoint, kpoints[kpoint]);
        write_kpoint_prefix(exx_output, kpoint, kpoints[kpoint]);
        write_kpoint_prefix(qsgw_output, kpoint, kpoints[kpoint]);

        for (int band = 0; band < live_band.get_n_bands(); ++band)
        {
            const double reference_occupation =
                reference_band.get_weight()[spin](kpoint, band) *
                occupation_scale;
            if (!std::isfinite(reference_occupation) ||
                reference_occupation < 0.0)
            {
                throw std::invalid_argument(
                    "QSGW band output reference occupation is invalid");
            }
            const double qsgw_energy =
                live_band.get_eigenvals()[spin](kpoint, band);
            if (!std::isfinite(qsgw_energy))
                throw std::invalid_argument(
                    "QSGW band output eigenvalue is non-finite");
            const double ks_energy = checked_real(
                h_ks(band, band), "KS Hamiltonian");
            const double exx_energy = checked_real(
                h_ks(band, band) - vxc(band, band) +
                    exx(band, band),
                "EXX Hamiltonian");

            ks_output << std::setw(15) << std::setprecision(5)
                      << reference_occupation
                      << std::setw(15) << std::setprecision(5)
                      << ks_energy * HA2EV;
            exx_output << std::setw(15) << std::setprecision(5)
                       << reference_occupation
                       << std::setw(15) << std::setprecision(5)
                       << exx_energy * HA2EV;
            qsgw_output << std::setw(15) << std::setprecision(5)
                        << reference_occupation
                        << std::setw(15) << std::setprecision(5)
                        << qsgw_energy * HA2EV;
        }
        ks_output << '\n';
        exx_output << '\n';
        qsgw_output << '\n';
    }
    if (!ks_output.good() || !exx_output.good() || !qsgw_output.good())
        throw std::runtime_error("Failed to write QSGW band output tables");
}

} // namespace qsgw
} // namespace librpa_int
