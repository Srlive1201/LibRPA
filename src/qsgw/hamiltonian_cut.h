#pragma once

#include "matrix_map.h"

#include "../core/meanfield.h"

namespace librpa_int
{
namespace qsgw
{

enum class HamiltonianCutMode
{
    Uncut = 0,
    ReferenceDiagonal = 1,
    ShiftedReferenceDiagonal = 2,
};

struct HamiltonianCutOptions
{
    int unoccupied_keep = 10;
    HamiltonianCutMode mode =
        HamiltonianCutMode::ShiftedReferenceDiagonal;
    double shift_ha = 20.0;
};

HamiltonianCutMode hamiltonian_cut_mode_from_int(int mode);

SpinKMatrixMap apply_hamiltonian_cut(
    const SpinKMatrixMap& hamiltonian,
    const SpinKMatrixMap& reference_hamiltonian,
    const MeanField& live_meanfield,
    const HamiltonianCutOptions& options);

} // namespace qsgw
} // namespace librpa_int
