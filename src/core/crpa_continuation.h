#pragma once

#include "crpa.h"

namespace librpa_int
{
struct CrpaContinuedTensor
{
    int atom_index, spin_left, spin_right, n_orbitals;
    std::vector<ComplexMatrix> values;
    std::vector<std::vector<int>> status;
    std::vector<std::vector<double>> spread;
    double max_source_residual = 0, max_holdout_residual = 0;
};

//! Native Thiele Padé on Delta U(z)=U(z)-v, z=i*nu. No fermionic shift.
//! Status bits: 1 primary failure, 2 subset failure, 4 subset spread,
//! 8 source/holdout residual, 16 insufficient nodes, 32 nonfinite evaluation.
//! These are numerical diagnostics, not passivity or material convergence certificates.
std::vector<CrpaContinuedTensor> continue_crpa_interactions(
    const CrpaResult& result, const std::vector<std::complex<double>>& queries, int n_parameters,
    double tolerance);
}  // namespace librpa_int
