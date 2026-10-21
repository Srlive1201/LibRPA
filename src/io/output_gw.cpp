#include "output_gw.h"

#include <iomanip>
#include <string>

#include "../utils/error.h"

namespace librpa_int
{

bool SelfEnergyOmegaWriter::open(const char *filename)
{
    output_.open(filename);
    return output_.good();
}

void SelfEnergyOmegaWriter::write_header(const std::vector<double> &frequencies, const int n_spins,
                                         const int n_kpoints, const int n_bands)
{
    n_frequencies_ = static_cast<int>(frequencies.size());
    output_ << frequencies.size() << " " << n_spins << " " << n_kpoints << " " << n_bands << "\n";
    for (const auto frequency : frequencies)
        output_ << std::scientific << std::setprecision(16) << frequency << "\n";
}

void SelfEnergyOmegaWriter::write_channel(const std::vector<cplxdb> &values, const int n_bands)
{
    for (int band = 0; band != n_bands; ++band)
        for (int frequency = 0; frequency != n_frequencies_; ++frequency)
        {
            const auto value = values[frequency * n_bands + band];
            output_ << std::scientific << std::setw(23) << std::setprecision(16) << value.real()
                    << " " << std::setw(23) << std::setprecision(16) << value.imag() << "\n";
        }
}

void SelfEnergyOmegaWriter::close() { output_.close(); }

void write_self_energy_omega_kpoints(const char *fn, const MpiCommHandler &comm_h,
                                     const std::vector<int> &iks)
{
    std::ofstream ofs;
    int file_ok = 1;
    if (comm_h.is_root())
    {
        ofs.open(fn);
        file_ok = ofs.good() ? 1 : 0;
    }
    comm_h.bcast(&file_ok, 1, 0);
    if (!file_ok)
        throw LIBRPA_RUNTIME_ERROR("failed to open self-energy k-index output file: " +
                                   std::string(fn));

    if (comm_h.is_root())
    {
        ofs << iks.size() << "\n";
        for (const int ik : iks)
        {
            ofs << ik << "\n";
        }
        ofs.close();
    }
}

}  // namespace librpa_int
