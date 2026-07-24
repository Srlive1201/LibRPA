#pragma once

#include "matrix_map.h"

#include "../core/meanfield.h"

namespace librpa_int
{
namespace qsgw
{

SpinKMatrixMap build_reference_hamiltonian(
    const MeanField& reference);

SpinKMatrixMap assemble_effective_hamiltonian(
    const SpinKMatrixMap& reference_hamiltonian,
    const SpinKMatrixMap& dft_vxc,
    const SpinKMatrixMap& exchange,
    const SpinKMatrixMap& correlation);

} // namespace qsgw
} // namespace librpa_int
