#pragma once

#include <fstream>
#include <vector>

#include "../math/matrix_m.h"
#include "../mpi/base_mpi.h"

namespace librpa_int
{

//! Serialize a gathered spin/k channel of KS-diagonal SigC(iw) values.
//! Values are frequency-major; the text layout is band-major.
class SelfEnergyOmegaWriter
{
public:
    bool open(const char *filename);
    void write_header(const std::vector<double> &frequencies, int n_spins, int n_kpoints,
                      int n_bands);
    void write_channel(const std::vector<cplxdb> &values, int n_bands);
    void close();

private:
    std::ofstream output_;
    int n_frequencies_ = 0;
};

//! Collective over comm_h. Rank 0 writes the supplied output k-point indices.
void write_self_energy_omega_kpoints(const char *filename, const MpiCommHandler &comm_h,
                                     const std::vector<int> &iks);

}  // namespace librpa_int
