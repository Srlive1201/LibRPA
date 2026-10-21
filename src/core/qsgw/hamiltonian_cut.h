#pragma once

#include "../meanfield.h"
#include "matrix_map.h"

namespace librpa_int
{
namespace qsgw
{

enum class HamiltonianCutMode
{
    Uncut = 0,
    ShiftedReferenceDiagonal = 1,
};

struct HamiltonianCutOptions
{
    int unoccupied_keep = 10;
    HamiltonianCutMode mode =
        HamiltonianCutMode::ShiftedReferenceDiagonal;
    double shift_ha = 0.00367493;
};

HamiltonianCutMode hamiltonian_cut_mode_from_int(int mode);

SpinKMatrixMap apply_hamiltonian_cut(
    const SpinKMatrixMap& hamiltonian,
    const SpinKMatrixMap& reference_hamiltonian,
    const MeanField& live_meanfield,
    const HamiltonianCutOptions& options);

} // namespace qsgw
} // namespace librpa_int
