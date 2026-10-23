#pragma once

#include "../headwing_velocity.h"
#include "fixed_basis.h"

namespace librpa_int::qsgw
{
using ::librpa_int::compute_localized_velocity_in_frame;
using ::librpa_int::HeadwingAoOperators;
using ::librpa_int::LocalizedMatrix;
using ::librpa_int::LocalizedVelocityResult;

// Model-control eigensolve. H/S/r conventions and the resulting velocity
// projection are shared with the native velocity kernel in core/headwing_velocity.h.
// Production QSGW refreshes retain their live entering E/C frame instead.
LocalizedVelocityResult solve_localized_velocity(const LocalizedMatrix& hamiltonian,
                                                 const HeadwingAoOperators& operators,
                                                 const Vector3_Order<double>& kpoint_fractional,
                                                 const Matrix3& lattice_bohr);

struct HamiltonianReconstructionResult
{
    std::vector<LocalizedMatrix> hamiltonian;
    double reference_metric_residual = 0.0;
    double fourier_roundtrip_residual_ha = 0.0;
    double native_dft_reference_residual_ha = 0.0;
    double native_dft_hermiticity_residual_ha = 0.0;
};

// Embed the complete fixed-reference operator as D h D^dagger, where
// D = S U0 (U0^dagger S U0)^-1 is the dual reference frame, then
// reconstruct its complex H(R). With native DFT baselines, the result is
// H_DFT(R) + IFT[H_current(k) - FT(H_DFT)(k)]; otherwise reconstruct whole H.
// The finite-BvK part uses equally weighted even-direction boundary images.
// Requires a complete unshifted mesh and the immutable square reference;
// neither live eigenpairs nor r(R) are needed. A supplied baseline must match
// the original DFT reference within the existing 1e-4 eV eigenframe gate.
HamiltonianReconstructionResult reconstruct_qsgw_hamiltonian(
    const MeanField& reference, const SpinKMatrixMap& current_hamiltonian,
    const HeadwingAoOperators& operators, const std::vector<Vector3_Order<double>>& kpoints,
    const Vector3_Order<int>& period, const Matrix3& lattice_bohr);

struct LiveVelocityResult
{
    VelocityMatrix velocity;
    double reference_metric_residual = 0.0;
    double live_metric_residual = 0.0;
    double hamiltonian_frame_residual_ha = 0.0;
    double fourier_roundtrip_residual_ha = 0.0;
    double native_dft_reference_residual_ha = 0.0;
    double native_dft_hermiticity_residual_ha = 0.0;
    double velocity_upper_diagonal_imaginary_ha_bohr = 0.0;
    // Complete reconstructed Hamiltonian, preserving an optional native DFT
    // baseline. The correction uses equally weighted BvK boundary images.
    // The velocity differentiates this same final H(R).
    std::vector<LocalizedMatrix> hamiltonian;
};

LiveVelocityResult compute_live_headwing_velocity(const MeanField& reference, const MeanField& live,
                                                  const SpinKMatrixMap& current_hamiltonian,
                                                  const HeadwingAoOperators& operators,
                                                  const std::vector<Vector3_Order<double>>& kpoints,
                                                  const Vector3_Order<int>& period,
                                                  const Matrix3& lattice_bohr);

}  // namespace librpa_int::qsgw
