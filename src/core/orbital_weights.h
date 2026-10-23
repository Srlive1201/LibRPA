#pragma once

#include "atomic_basis.h"
#include "meanfield.h"

namespace librpa_int
{
struct OrbitalWeightGroup
{
    int atom;
    int type;
    int l;
    std::string channel;
};

struct OrbitalWeightState
{
    double energy;
    double occupation;
    double kweight;
    double norm;
    std::vector<double> mulliken_ao, lowdin_ao;
    std::vector<double> mulliken, lowdin;
};

struct OrbitalWeightResult
{
    int n_spins, n_spinor, n_kpoints, n_bands, n_aos;
    double fermi;
    double highest_occupied;
    double min_overlap_eigenvalue;
    double max_norm_error = 0;
    double max_sum_error = 0;
    double max_sqrt_residual = 0;
    std::vector<OrbitalWeightGroup> groups;
    std::vector<OrbitalWeightState> states;  // [spin][k][band]
};

// Uses existing scalar/spinor coefficients. Only S is diagonalized; no H/KS solve.
// Spinor components are summed as charge weights within one stored spin channel.
// All radial functions contribute to atom/l groups; d splits into global eg/t2g.
OrbitalWeightResult analyze_orbital_weights(
    const AtomicBasis& basis, const BasisConvention& convention, const std::vector<int>& atom_types,
    const MeanField& meanfield, const std::vector<ComplexMatrix>& overlap,
    const std::vector<double>& kweights, double residual_tol);

// Grid, poles and sigma have the same energy unit. Density uses its inverse.
// Stored spin channels remain separate; no occupation or degeneracy factor.
void gaussian_orbital_pdos(const OrbitalWeightResult& result, const std::vector<double>& grid,
                           double sigma, std::vector<double>& mulliken,
                           std::vector<double>& lowdin);
}  // namespace librpa_int
