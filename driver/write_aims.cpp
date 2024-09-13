#include "write_aims.h"

#include "gw.h"

#include <fstream>
#include <iomanip>

void write_self_energy_omega(const char *fn, const LIBRPA::G0W0 &s_g0w0)
{
    const auto &n_spins = s_g0w0.mf.get_n_spins();
    const auto &n_kpts = s_g0w0.mf.get_n_kpoints();
    const auto &n_bands = s_g0w0.mf.get_n_bands();
    const auto &freqs = s_g0w0.tfg.get_freq_nodes();
    std::ofstream ofs;
    ofs.open(fn);
    ofs << freqs.size() << " "
        << n_spins << " "
        << n_kpts  << " "
        << n_bands << "\n";
    for (const auto& freq: freqs)
    {
        ofs << fixed << setprecision(12) << freq << "\n";
    }

    for (int ispin = 0; ispin != n_spins; ispin++)
    {
        for (int ik = 0; ik != n_kpts; ik++)
        {
            for (int ib = 0; ib != n_bands; ib++)
            {
                for (const auto& freq: freqs)
                {
                    const auto &sigc_mat = s_g0w0.sigc_is_ik_f_KS.at(ispin).at(ik).at(freq);
                    const auto sigc = sigc_mat(ib, ib);
                    // NOTE: sprintf in GCC 14 somehow may lead to segfault. Use ofstream instead.
                    ofs << fixed << setw(23) << setprecision(16) << sigc.real() << setw(23) << setprecision(16) << sigc.imag() << "\n";
                }
            }
        }
    }
    ofs.close();
}
