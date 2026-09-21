#pragma once

#include "matrix_map.h"

namespace librpa_int
{
namespace qsgw
{

struct SpinKHamiltonianResidual
{
    double l2 = 0.0;
    double maximum = 0.0;
};

// Norms of H_out - H_in on the SCF grid, before mixing.
SpinKHamiltonianResidual measure_spin_k_hamiltonian_residual(
    const SpinKMatrixMap& output,
    const SpinKMatrixMap& input);

// H_next = H_in + beta * (H_out - H_in); beta=1 accepts H_out.
// The driver owns the current Hamiltonians and uses the same beta on grid/path.
SpinKMatrixMap mix_spin_k_hamiltonian(
    const SpinKMatrixMap& input,
    const SpinKMatrixMap& output,
    double beta);

} // namespace qsgw
} // namespace librpa_int
