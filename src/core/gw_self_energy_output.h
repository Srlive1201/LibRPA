#pragma once

#include <algorithm>
#include <numeric>
#include <string>
#include <vector>

#include "../io/output_gw.h"
#include "../utils/error.h"
#include "gw.h"

namespace librpa_int
{

inline void write_self_energy_omega(const char *fn, const G0W0 &s_g0w0, const std::vector<int> &iks,
                                    const int n_bands)
{
    const auto &comm_h = s_g0w0.comm_h;
    const int n_spins = s_g0w0.mf.get_n_spins();
    const auto &freqs = s_g0w0.tfg.get_freq_nodes();
    const int n_kpts = static_cast<int>(iks.size());

    SelfEnergyOmegaWriter writer;
    int file_ok = 1;
    if (comm_h.is_root())
    {
        file_ok = writer.open(fn) ? 1 : 0;
    }
    comm_h.bcast(&file_ok, 1, 0);
    if (!file_ok)
        throw LIBRPA_RUNTIME_ERROR("failed to open self-energy output file: " + std::string(fn));

    if (comm_h.is_root()) writer.write_header(freqs, n_spins, n_kpts, n_bands);

    const int nfreq = static_cast<int>(freqs.size());
    const int n_value = nfreq * n_bands;
    std::vector<int> owners(nfreq, 0);
    std::vector<int> owners_sum(nfreq, 0);
    std::vector<cplxdb> values(n_value, cplxdb{0.0, 0.0});
    std::vector<cplxdb> values_sum(n_value, cplxdb{0.0, 0.0});

    const auto diag_index = [n_bands](int ifreq, int ib) { return ifreq * n_bands + ib; };

    for (int ispin = 0; ispin != n_spins; ispin++)
    {
        for (const int ik : iks)
        {
            std::fill(owners.begin(), owners.end(), 0);
            std::fill(owners_sum.begin(), owners_sum.end(), 0);
            std::fill(values.begin(), values.end(), cplxdb{0.0, 0.0});
            std::fill(values_sum.begin(), values_sum.end(), cplxdb{0.0, 0.0});

            const auto it_sp = s_g0w0.sigc_diag_is_ik_f_KS.find(ispin);
            if (it_sp != s_g0w0.sigc_diag_is_ik_f_KS.cend())
            {
                const auto it_k = it_sp->second.find(ik);
                if (it_k != it_sp->second.cend())
                {
                    for (int ifreq = 0; ifreq != nfreq; ++ifreq)
                    {
                        const double freq = freqs[ifreq];
                        const auto it_freq = it_k->second.find(freq);
                        if (it_freq == it_k->second.cend()) continue;

                        const auto &sigc_diag = it_freq->second;
                        if (static_cast<int>(sigc_diag.size()) != n_bands) continue;

                        owners[ifreq] = 1;
                        for (int ib = 0; ib != n_bands; ib++)
                        {
                            values[diag_index(ifreq, ib)] = sigc_diag[ib];
                        }
                    }
                }
            }

            comm_h.allreduce(owners.data(), owners_sum.data(), nfreq, MPI_SUM);
            comm_h.allreduce(values.data(), values_sum.data(), n_value, MPI_SUM);

            for (int ifreq = 0; ifreq != nfreq; ++ifreq)
            {
                if (owners_sum[ifreq] != 1)
                {
                    throw LIBRPA_RUNTIME_ERROR(
                        "failed to locate a unique SigC matrix owner for spin = " +
                        std::to_string(ispin) + " ik = " + std::to_string(ik) +
                        " freq = " + std::to_string(freqs[ifreq]) +
                        ", owner count = " + std::to_string(owners_sum[ifreq]));
                }
            }

            if (!comm_h.is_root()) continue;

            writer.write_channel(values_sum, n_bands);
        }
    }
    if (comm_h.is_root()) writer.close();
}

inline void write_self_energy_omega(const char *fn, const G0W0 &s_g0w0, const int n_kpts,
                                    const int n_bands)
{
    std::vector<int> iks(n_kpts);
    std::iota(iks.begin(), iks.end(), 0);
    write_self_energy_omega(fn, s_g0w0, iks, n_bands);
}

}  // namespace librpa_int
