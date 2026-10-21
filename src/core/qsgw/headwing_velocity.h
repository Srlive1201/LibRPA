#pragma once

#include <array>
#include <map>
#include <vector>

#include "../../math/matrix3.h"
#include "../../math/vector3_order.h"
#include "fixed_basis.h"

namespace librpa_int::qsgw
{

// Localized matrices with independent real-space supports. H is full complex;
// S and r use the upper-triangle Fourier completion specified below. Matrices
// may be row-major or column-major. No filenames or reader units enter here.
using LocalizedMatrix = std::map<Vector3_Order<int>, Matz>;

struct HeadwingAoOperators
{
    int dimension = 0;
    LocalizedMatrix overlap;
    std::array<LocalizedMatrix, 3> position;
};

struct LocalizedVelocityResult
{
    std::vector<double> eigenvalues_ha;
    // AO rows, band columns: H C = S C E and C^dagger S C = I.
    Matz ket_coefficients;
    Matz overlap;
    std::array<Matz, 3> velocity_ha_bohr;
    double metric_residual = 0.0;
    double hamiltonian_frame_residual_ha = 0.0;
    double hamiltonian_hermiticity_residual_ha = 0.0;
    double velocity_upper_diagonal_imaginary_ha_bohr = 0.0;
};

// Model-only native interface. H is in Ha, S is dimensionless, r and the row
// lattice vectors are in Bohr. k is fractional and the Fourier phase is
// exp(+i 2*pi*k.R); derivatives are with respect to Cartesian k in Bohr^-1.
// At each k, native S and dS completion takes their upper triangles with a real
// diagonal and conjugate lower triangle. Position is completed as r_ij =
// conj(raw_r_ji + i*dS_ji) for i >= j. Band velocity uses the Hermitian upper
// triangle and an explicitly real diagonal. These inherited real-diagonal
// conventions agree with executed PyATB on physically Hermitian models;
// PyATB's upper selfadjoint view preserves invalid complex diagonal entries,
// which are outside the claimed model consistency scope.
// The generalized eigensolve supplies matching E/C/v in a single frame.
LocalizedVelocityResult solve_localized_velocity(const LocalizedMatrix& hamiltonian,
                                                 const HeadwingAoOperators& operators,
                                                 const Vector3_Order<double>& kpoint_fractional,
                                                 const Matrix3& lattice_bohr);

// Evaluate the same model in an already selected complete eigenframe. This
// preserves the live QSGW C/E frame instead of introducing a new solve and
// arbitrary phases or rotations inside a degenerate subspace.
LocalizedVelocityResult compute_localized_velocity_in_frame(
    const LocalizedMatrix& hamiltonian, const HeadwingAoOperators& operators,
    const Vector3_Order<double>& kpoint_fractional, const Matrix3& lattice_bohr,
    const std::vector<double>& eigenvalues_ha, const Matz& ket_coefficients);

}  // namespace librpa_int::qsgw
