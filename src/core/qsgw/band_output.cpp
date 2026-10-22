#include "band_output.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <ostream>
#include <stdexcept>
#include <string>

#include "../../utils/constants.h"

namespace librpa_int
{
namespace qsgw
{
namespace
{

void validate_layout(const MeanField& live, const MeanField& reference,
                     const std::vector<Vector3_Order<double>>& kpoints, const int spin,
                     const double chemical_potential)
{
    if (!live.initialized() || !reference.initialized())
        throw std::invalid_argument("QSGW band output mean fields must be initialized");
    if (live.get_n_spins() != reference.get_n_spins() ||
        live.get_n_kpoints() != reference.get_n_kpoints() ||
        live.get_n_bands() != reference.get_n_bands() ||
        live.get_n_spinor() != reference.get_n_spinor())
        throw std::invalid_argument("QSGW band output mean-field layouts differ");
    if (spin < 0 || spin >= live.get_n_spins())
        throw std::invalid_argument("QSGW band output spin index is out of range");
    if (kpoints.size() != static_cast<std::size_t>(live.get_n_kpoints()))
        throw std::invalid_argument("QSGW band output k-point count differs");
    if (!std::isfinite(chemical_potential))
        throw std::invalid_argument("QSGW band output chemical potential is non-finite");
}

const Matz& matrix_at(const SpinKMatrixMap& matrices, const int spin, const int kpoint,
                      const int dimension, const char* label)
{
    const auto spin_it = matrices.find(spin);
    if (spin_it == matrices.end())
        throw std::invalid_argument(std::string("QSGW band output is missing ") + label +
                                    " spin data");
    const auto kpoint_it = spin_it->second.find(kpoint);
    if (kpoint_it == spin_it->second.end() || kpoint_it->second.nr() != dimension ||
        kpoint_it->second.nc() != dimension)
        throw std::invalid_argument(std::string("QSGW band output has malformed ") + label +
                                    " k-point data");
    return kpoint_it->second;
}

double checked_real(const cplxdb value, const char* label)
{
    if (!std::isfinite(value.real()) || !std::isfinite(value.imag()))
        throw std::invalid_argument(std::string("QSGW band output has invalid ") + label +
                                    " diagonal data");
    return value.real();
}

void write_kpoint_prefix(std::ostream& output, const int kpoint,
                         const Vector3_Order<double>& coordinate)
{
    output << std::setw(5) << kpoint + 1 << std::setw(15) << std::setprecision(7) << coordinate.x
           << std::setw(15) << std::setprecision(7) << coordinate.y << std::setw(15)
           << std::setprecision(7) << coordinate.z;
}

}  // namespace

void write_qsgw_band_spin_tables(std::ostream& ks_output, std::ostream& exx_output,
                                 std::ostream& qsgw_output, const MeanField& live_band,
                                 const MeanField& reference_band,
                                 const std::vector<Vector3_Order<double>>& kpoints,
                                 const SpinKMatrixMap& reference_hamiltonian,
                                 const SpinKMatrixMap& dft_vxc, const SpinKMatrixMap& exchange,
                                 const int spin, const double chemical_potential_ha)
{
    validate_layout(live_band, reference_band, kpoints, spin, chemical_potential_ha);
    if (!ks_output.good() || !exx_output.good() || !qsgw_output.good())
        throw std::invalid_argument("QSGW band output stream is not writable");

    ks_output << std::fixed;
    exx_output << std::fixed;
    qsgw_output << std::fixed;
    const double occupation_scale = static_cast<double>(live_band.get_n_kpoints());
    for (int kpoint = 0; kpoint < live_band.get_n_kpoints(); ++kpoint)
    {
        const Matz& h_ks = matrix_at(reference_hamiltonian, spin, kpoint, live_band.get_n_bands(),
                                     "KS Hamiltonian");
        const Matz& vxc = matrix_at(dft_vxc, spin, kpoint, live_band.get_n_bands(), "Vxc");
        const Matz& exx = matrix_at(exchange, spin, kpoint, live_band.get_n_bands(), "EXX");
        write_kpoint_prefix(ks_output, kpoint, kpoints[kpoint]);
        write_kpoint_prefix(exx_output, kpoint, kpoints[kpoint]);
        write_kpoint_prefix(qsgw_output, kpoint, kpoints[kpoint]);

        for (int band = 0; band < live_band.get_n_bands(); ++band)
        {
            const double reference_occupation =
                reference_band.get_weight()[spin](kpoint, band) * occupation_scale;
            if (!std::isfinite(reference_occupation) || reference_occupation < 0.0)
            {
                throw std::invalid_argument("QSGW band output reference occupation is invalid");
            }
            const double qsgw_energy = live_band.get_eigenvals()[spin](kpoint, band);
            if (!std::isfinite(qsgw_energy))
                throw std::invalid_argument("QSGW band output eigenvalue is non-finite");
            const double ks_energy = checked_real(h_ks(band, band), "KS Hamiltonian");
            const double exx_energy = checked_real(
                h_ks(band, band) - vxc(band, band) + exx(band, band), "EXX Hamiltonian");

            ks_output << std::setw(15) << std::setprecision(5) << reference_occupation
                      << std::setw(15) << std::setprecision(5) << ks_energy * HA2EV;
            exx_output << std::setw(15) << std::setprecision(5) << reference_occupation
                       << std::setw(15) << std::setprecision(5) << exx_energy * HA2EV;
            qsgw_output << std::setw(15) << std::setprecision(5) << reference_occupation
                        << std::setw(15) << std::setprecision(5) << qsgw_energy * HA2EV;
        }
        ks_output << '\n';
        exx_output << '\n';
        qsgw_output << '\n';
    }
    if (!ks_output.good() || !exx_output.good() || !qsgw_output.good())
        throw std::runtime_error("Failed to write QSGW band output tables");
}

}  // namespace qsgw
}  // namespace librpa_int

namespace librpa_int::qsgw
{
double occupation_number_from_kpoint_weight(const double weighted_occupation,
                                            const double kpoint_weight)
{
    if (!std::isfinite(kpoint_weight) || kpoint_weight <= 0.0)
        throw std::invalid_argument("k-point weight must be finite and positive");
    // SCF occupations include the weight of their loaded (possibly IBZ) k point.
    return weighted_occupation / kpoint_weight;
}

void write_energy_qp(std::ostream &ofs, const librpa_int::MeanField &mf,
                     const std::vector<librpa_int::Vector3_Order<double>> &kfrac_output,
                     const std::vector<int> &output_to_input_kpoint,
                     const std::vector<librpa_int::matrix> &vxc, const std::vector<double> &vexx,
                     const std::vector<librpa_int::cplxdb> &sigc, const int n_kpoints_data,
                     const int i_state_low, const int n_states_calc,
                     const std::vector<double> &kpoint_weights)
{
    if (kpoint_weights.size() != static_cast<std::size_t>(mf.get_n_kpoints()))
        throw std::invalid_argument("k-point weight count does not match mean-field data");
    const std::string sep =
        "-----------------------------------------"
        "----------------------------------------------------------";
    ofs << "  state     occ_num        e_gs(Ha)        e_qp(Ha)" << std::endl;
    ofs << sep << std::endl;
    for (int i_kpoint = 0; i_kpoint < static_cast<int>(kfrac_output.size()); i_kpoint++)
    {
        const int i_kpoint_input = output_to_input_kpoint.empty()
                                       ? i_kpoint
                                       : output_to_input_kpoint[static_cast<size_t>(i_kpoint)];
        const auto &k = kfrac_output[static_cast<size_t>(i_kpoint)];
        for (int i_spin = 0; i_spin < mf.get_n_spins(); i_spin++)
        {
            if (mf.get_n_spins() == 2)
            {
                ofs << std::setw(35) << "" << (i_spin == 0 ? "Spin Up" : "Spin Down") << std::endl;
            }
            const size_t start_k = (i_spin * n_kpoints_data + i_kpoint_input) * n_states_calc;
            ofs << "  K_point " << std::setw(4) << i_kpoint + 1 << " : " << std::fixed
                << std::setprecision(4) << std::setw(16) << k.x << std::setw(16) << k.y
                << std::setw(16) << k.z << std::endl;
            ofs << sep << std::endl;
            for (int i = 0; i < n_states_calc; i++)
            {
                const int i_state = i + i_state_low;
                const auto occ_state = occupation_number_from_kpoint_weight(
                    mf.get_weight()[i_spin](i_kpoint_input, i_state),
                    kpoint_weights.at(static_cast<std::size_t>(i_kpoint_input)));
                const auto eks_state = mf.get_eigenvals()[i_spin](i_kpoint_input, i_state);
                const auto eqp = eks_state - vxc[i_spin](i_kpoint_input, i_state) +
                                 vexx[start_k + i] + sigc[start_k + i].real();
                ofs << "  " << std::setw(6) << i_state + 1 << "  " << std::fixed
                    << std::setprecision(4) << std::setw(8) << occ_state << std::scientific
                    << std::uppercase << std::setprecision(10) << std::setw(20) << eks_state
                    << std::setw(20) << eqp << std::endl;
            }
            if (mf.get_n_spins() == 2 && i_spin == 0)
            {
                ofs << sep << std::endl;
            }
        }
        ofs << sep << std::endl;
        ofs << std::endl;
    }
    if (!ofs) throw std::runtime_error("Cannot write quasiparticle energy output");
}

void write_gw_band_spin_table(std::ostream &output, const librpa_int::MeanField &mf,
                              const std::vector<librpa_int::Vector3_Order<double>> &kfrac,
                              const std::vector<librpa_int::matrix> &vxc,
                              const std::vector<double> &vexx,
                              const std::vector<librpa_int::cplxdb> &sigc, const int spin,
                              const int i_state_low, const int n_states_calc)
{
    using namespace librpa_int;
    const int nk = mf.get_n_kpoints();
    if (spin < 0 || spin >= mf.get_n_spins() || kfrac.size() != static_cast<std::size_t>(nk) ||
        i_state_low < 0 || n_states_calc <= 0 || i_state_low + n_states_calc > mf.get_n_states())
        throw std::invalid_argument("Invalid GW band table layout");
    const auto count = static_cast<std::size_t>(mf.get_n_spins()) * nk * n_states_calc;
    if (vxc.size() != static_cast<std::size_t>(mf.get_n_spins()) || vxc[spin].nr != nk ||
        vxc[spin].nc != mf.get_n_states() || vexx.size() != count || sigc.size() != count)
        throw std::invalid_argument("Incomplete GW band table energies");

    output << std::fixed;
    for (int ik = 0; ik < nk; ++ik)
    {
        const auto &k = kfrac[ik];
        const std::size_t start = (spin * nk + ik) * n_states_calc;
        output << std::setw(5) << ik + 1 << std::setw(15) << std::setprecision(7) << k.x
               << std::setw(15) << k.y << std::setw(15) << k.z;
        for (int i = 0; i < n_states_calc; ++i)
        {
            const int state = i_state_low + i;
            const double occupation = mf.get_weight()[spin](ik, state) * nk;
            const double eigenvalue = mf.get_eigenvals()[spin](ik, state) * HA2EV;
            const double potential = vxc[spin](ik, state) * HA2EV;
            const double exchange = vexx[start + i] * HA2EV;
            const double correlation = sigc[start + i].real() * HA2EV;
            const double eqp = eigenvalue - potential + exchange + correlation;
            output << std::setw(15) << std::setprecision(5) << occupation << std::setw(15) << eqp;
        }
        output << '\n';
    }
    output.flush();
    if (!output) throw std::runtime_error("Cannot write GW band table");
}
}  // namespace librpa_int::qsgw
