#include "task_helper.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <string>

#include "../../src/io/global_io.h"
#include "../../src/utils/constants.h"

void print_band_analysis(const librpa_int::MeanField &mf,
                         const std::vector<librpa_int::Vector3_Order<double>> &kfrac_output,
                         const std::vector<int> &output_to_input_kpoint,
                         const std::vector<librpa_int::matrix> &vxc,
                         const std::vector<double> &vexx,
                         const std::vector<librpa_int::cplxdb> &sigc,
                         const int i_state_low, const int n_states_calc)
{
    using namespace librpa_int;
    using global::lib_printf;
    constexpr auto output_level = LIBRPA_VERBOSE_CRITICAL;
    const int nk = mf.get_n_kpoints();
    if (nk == 0 || kfrac_output.empty()) return;

    // Retain the original occupied/empty band partition, separately for each spin.
    for (int isp = 0; isp < mf.get_n_spins(); ++isp)
    {
        lib_printf(output_level, "Band analysis, spin %2d\n", isp + 1);
        const int nocc = mf.find_highest_occupied_state(isp, 0).second + 1;
        bool fixed_occupation = true;
        for (int ik = 1; ik < nk; ++ik)
            if (mf.find_highest_occupied_state(isp, ik).second + 1 != nocc)
                fixed_occupation = false;
        if (!fixed_occupation)
        {
            lib_printf(output_level,
                       "Band analysis unavailable: occupied-band count varies across k-points.\n");
            continue;
        }
        lib_printf(output_level, "Bands of occupation: %4d\n", nocc);

        // Higher empty bands may acquire spuriously low QP energies from RI errors.
        // Check them without including them in the CBM search or changing occupations.
        size_t n_anomalous = 0, ik_lowest = 0;
        int state_lowest = -1;
        double lowest_qp = mf.get_efermi();
        for (size_t ik_out = 0; ik_out < kfrac_output.size(); ++ik_out)
        {
            const int ik = output_to_input_kpoint.empty()
                ? static_cast<int>(ik_out) : output_to_input_kpoint[ik_out];
            const size_t start_k = (isp * nk + ik) * n_states_calc;
            for (int state = std::max(i_state_low, nocc + 1);
                 state < i_state_low + n_states_calc; ++state)
            {
                const size_t index = start_k + state - i_state_low;
                const double eqp = mf.get_eigenvals()[isp](ik, state)
                    - vxc[isp](ik, state) + vexx[index] + sigc[index].real();
                if (!std::isfinite(eqp) || eqp >= mf.get_efermi()) continue;
                ++n_anomalous;
                if (eqp < lowest_qp)
                {
                    lowest_qp = eqp;
                    ik_lowest = ik_out;
                    state_lowest = state;
                }
            }
        }
        if (n_anomalous > 0)
        {
            const auto &k = kfrac_output[ik_lowest];
            lib_printf(LIBRPA_VERBOSE_WARN,
                       "Warning! Spin %d: %zu higher empty GW states lie below the mean-field Fermi level"
                       " (%.7f eV); possible RI error. These states are excluded from the CBM search.\n",
                       isp + 1, n_anomalous, mf.get_efermi() * HA2EV);
            lib_printf(LIBRPA_VERBOSE_WARN,
                       "Lowest offending state: band %d, k-point %zu: (%.5f, %.5f, %.5f), e_qp = %.7f eV\n",
                       state_lowest + 1, ik_lowest + 1, k.x, k.y, k.z, lowest_qp * HA2EV);
        }

        if (nocc == 0 || nocc >= mf.get_n_states() || nocc - 1 < i_state_low
            || nocc >= i_state_low + n_states_calc)
        {
            lib_printf(output_level,
                       "Band analysis unavailable: calculated state window does not contain both band edges.\n");
            continue;
        }

        const char *labels[] = {"GW", "EXX", "DFT"};
        for (int method = 0; method < 3; ++method)
        {
            if (method > 0) lib_printf(output_level, "\n");
            double valence = -std::numeric_limits<double>::infinity();
            double conduct = std::numeric_limits<double>::infinity();
            int ik_val = -1, ik_cond = -1;
            bool finite = true;
            for (size_t ik_out = 0; ik_out < kfrac_output.size(); ++ik_out)
            {
                const int ik = output_to_input_kpoint.empty()
                    ? static_cast<int>(ik_out) : output_to_input_kpoint[ik_out];
                const size_t start_k = (isp * nk + ik) * n_states_calc;
                // Search all calculated valence bands, but retain the first empty band as CBM.
                for (int state = i_state_low; state <= nocc; ++state)
                {
                    double energy = mf.get_eigenvals()[isp](ik, state);
                    const size_t index = start_k + state - i_state_low;
                    if (method < 2) energy += vexx[index] - vxc[isp](ik, state);
                    if (method == 0) energy += sigc[index].real();
                    energy *= HA2EV;
                    if (!std::isfinite(energy)) finite = false;
                    if (state < nocc && energy > valence)
                    {
                        valence = energy;
                        ik_val = static_cast<int>(ik_out);
                    }
                    if (state == nocc && energy < conduct)
                    {
                        conduct = energy;
                        ik_cond = static_cast<int>(ik_out);
                    }
                }
            }
            if (!finite || ik_val < 0 || ik_cond < 0)
            {
                lib_printf(output_level, "%s band analysis unavailable: non-finite energy in VBM/CBM search window.\n",
                           labels[method]);
                continue;
            }
            const auto &kv = kfrac_output[ik_val];
            const auto &kc = kfrac_output[ik_cond];
            lib_printf(output_level, "%s VBM: k-point %4d: (%.5f, %.5f, %.5f)\n",
                       labels[method], ik_val + 1, kv.x, kv.y, kv.z);
            lib_printf(output_level, "%s CBM: k-point %4d: (%.5f, %.5f, %.5f)\n",
                       labels[method], ik_cond + 1, kc.x, kc.y, kc.z);
            lib_printf(output_level, "%s bandgap(eV): %12.7f\n", labels[method], conduct - valence);
        }
    }
}

void write_energy_qp(const librpa_int::MeanField &mf,
                     const std::vector<librpa_int::Vector3_Order<double>> &kfrac_output,
                     const std::vector<int> &output_to_input_kpoint,
                     const std::vector<librpa_int::matrix> &vxc, const std::vector<double> &vexx,
                     const std::vector<librpa_int::cplxdb> &sigc, const int n_kpoints_data,
                     const int i_state_low, const int n_states_calc, const double occupation_scale)
{
    const std::string sep =
        "-----------------------------------------"
        "----------------------------------------------------------";
    std::ofstream ofs("energy_qp");
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
                ofs << std::setw(35) << "" << (i_spin == 0 ? "Spin Up" : "Spin Down")
                    << std::endl;
            }
            const size_t start_k = (i_spin * n_kpoints_data + i_kpoint_input) * n_states_calc;
            ofs << "  K_point " << std::setw(4) << i_kpoint + 1 << " : " << std::fixed
                << std::setprecision(4) << std::setw(16) << k.x << std::setw(16) << k.y
                << std::setw(16) << k.z << std::endl;
            ofs << sep << std::endl;
            for (int i = 0; i < n_states_calc; i++)
            {
                const int i_state = i + i_state_low;
                const auto occ_state =
                    mf.get_weight()[i_spin](i_kpoint_input, i_state) * occupation_scale;
                const auto eks_state = mf.get_eigenvals()[i_spin](i_kpoint_input, i_state);
                const auto eqp = eks_state - vxc[i_spin](i_kpoint_input, i_state) +
                                 vexx[start_k+i] + sigc[start_k+i].real();
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
}
