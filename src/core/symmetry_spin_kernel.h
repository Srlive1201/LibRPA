/*!
 * @file symmetry_spin_kernel.h
 * @brief Generic four-block spinor bilinear kernel (Phase 3, library-only).
 *
 * This kernel implements the spinor extension of the AO bilinear transform
 * for one spin-space operation (g, U_s, eta):
 *
 *   unitary:      X' = (D orb U) X (D orb U)^dagger
 *   antiunitary:  X' = J_AO [(D orb U) X (D orb U)^dagger]* J_AO^dagger,
 *                 J_AO = I orb i sigma_y
 *
 * Convention notes (locked in Phase 0):
 * - C4: the LibRPA orbital kernel is source-to-target, X' = D X D^dagger
 *   (dense A X A^dagger, or block form M_I^T X conj(M_J)). The SU(2) mixing
 *   below follows the same convention: Y^{ab} = sum_cd U_ac conj(U_bd) G^{cd}
 *   with U row-major [u00 u01; u10 u11]. Do not rewrite as D^dagger X D.
 * - C8 layout: LibRPA keeps spin channels "channel-outermost" — a spinor
 *   operator is four independent orbital matrices X^{ab}, where a is the bra
 *   (row) spin and b the ket (column) spin. Merged into one dense 2N x 2N
 *   matrix, element (i, a; j, b) sits at row a*N+i, col b*N+j (orbital index
 *   fastest). ABACUS uses the opposite spin-fast interleaved layout
 *   I = 2*iw + s; shuffling the merged matrix with I = 2*iw + s reproduces
 *   the ABACUS layout. All formulas here use the channel-outermost order.
 * - Time reversal is applied after the spatial unitary action (Theta M
 *   Theta^-1 = M*); the order does not commute. The remap
 *   {Y11*, -Y10*, -Y01*, Y00*} is an involution (sigma_y* = -sigma_y).
 *
 * Production use (Phase 4): the spinor Green's-function k-star restore in
 * meanfield.cpp / meanfield_mpi.cpp drives this kernel through the call
 * sites that are enabled only for n_spinor == 2 input. The scalar
 * (n_spinor == 1) restore path does not include this header.
 */
#pragma once

#include <array>
#include <cassert>
#include <complex>

#include "symmetry_types.h"

namespace librpa_int
{

//! Four spin blocks of one AO bilinear operator, channel-outermost (C8).
template <class Matrix>
struct SpinorBlocks4
{
    Matrix b00, b01, b10, b11;
};

//! Bilinear transform convention. Only the confirmed C4 form exists so far.
enum class BilinearConvention
{
    //! X' = D X D^dagger, source-to-target (LibRPA convention C4).
    SourceToTarget_DXDdag
};

/*!
 * @brief Apply one spin-space operation to a four-block AO bilinear operator.
 *
 * @param op     spin-space operation (spatial_id, spin_u, antiunitary)
 * @param in     source four blocks; missing blocks must be zero-filled by the caller
 * @param orbit  orbital rotation callback reusing the existing scalar kernel,
 *               orbit(spatial_id, X) -> A X A^dagger (dense) or the block form
 *               M_I^T X conj(M_J)
 * @param conv   bilinear convention; only SourceToTarget_DXDdag is supported
 *
 * Steps: (1) per-block spatial rotation G^{cd} = orbit(spatial_id, X^{cd});
 * (2) SU(2) mixing Y^{ab} = sum_cd U_ac conj(U_bd) G^{cd}; (3) when
 * antiunitary, the Theta remap {conj(Y11), -conj(Y10), -conj(Y01), conj(Y00)}.
 * Fast path: Identity spin source with U = I and eta = 0 skips steps 2-3 and
 * returns the step-1 result unchanged, preserving the old scalar-kernel cost.
 */
template <class Matrix, class OrbitalRotate>
SpinorBlocks4<Matrix> transform_spinor_bilinear(
    const SymmetrySpinOperation &op,
    const SpinorBlocks4<Matrix> &in,
    OrbitalRotate &&orbit,
    BilinearConvention conv)
{
    assert(conv == BilinearConvention::SourceToTarget_DXDdag);
    (void)conv;

    // step 1: per-block spatial rotation (existing scalar kernel)
    SpinorBlocks4<Matrix> g{orbit(op.spatial_id, in.b00), orbit(op.spatial_id, in.b01),
                            orbit(op.spatial_id, in.b10), orbit(op.spatial_id, in.b11)};

    // A supplied matrix is authoritative even with the Identity source tag.
    // Only an actually trivial spin action can skip the four-block mixing.
    if (op.spin_source == SymmetrySpinActionSource::Identity && !op.antiunitary &&
        op.spin_u == SymmetrySpinOperation{}.spin_u)
    {
        return g;
    }

    // step 2: SU(2) mixing, Y^{ab} = sum_cd U[a][c] conj(U[b][d]) G^{cd}
    const auto &U = op.spin_u; // row-major [u00 u01; u10 u11]
    const auto u = [&U](int r, int c) -> const std::complex<double> & {
        return U[2 * r + c];
    };
    const auto mix = [&u](int a, int b, const SpinorBlocks4<Matrix> &blocks) -> Matrix {
        Matrix out = u(a, 0) * std::conj(u(b, 0)) * blocks.b00;
        out += u(a, 0) * std::conj(u(b, 1)) * blocks.b01;
        out += u(a, 1) * std::conj(u(b, 0)) * blocks.b10;
        out += u(a, 1) * std::conj(u(b, 1)) * blocks.b11;
        return out;
    };
    SpinorBlocks4<Matrix> y{mix(0, 0, g), mix(0, 1, g), mix(1, 0, g), mix(1, 1, g)};

    if (!op.antiunitary)
    {
        return y;
    }

    // step 3: Theta remap sigma_y Y* sigma_y (involution; spatial first, Theta
    // last, the order does not commute). X' = J_AO Y* J_AO^dagger gives
    // X'00 = Y11*, X'01 = -Y10*, X'10 = -Y01*, X'11 = Y00*.
    const std::complex<double> minus_one(-1.0, 0.0);
    return {conj(y.b11), minus_one * conj(y.b10),
            minus_one * conj(y.b01), conj(y.b00)};
}

//! How a spin action acts on a collinear two-channel storage.
enum class CollinearChannelAction
{
    Keep,         //!< U diagonal (up to phase): channels evolve independently
    Swap,         //!< U off-diagonal: the two channels are exchanged
    Incompatible  //!< genuine spinor action; collinear storage cannot hold it
};

/*!
 * @brief Classify the collinear action of one spin operation.
 *
 * Builds U P_up U^dagger with P_up = diag(1, 0) and compares against P_up and
 * P_dn. The projector form removes the SU(2) double-cover sign (+-U give the
 * same classification). A diagonal result means Keep, an off-diagonal (fully
 * swapped) result means Swap, anything else is Incompatible and the caller
 * must upgrade to Spinor2 storage or reject the operation list.
 *
 * `antiunitary` is part of the calling convention only: the complex
 * conjugation of an antiunitary operation is applied by the Theta remap in
 * transform_spinor_bilinear, while the channel classification itself depends
 * on U alone.
 */
CollinearChannelAction classify_collinear_action(
    const std::array<std::complex<double>, 4> &U,
    bool antiunitary,
    double tol);

/*!
 * @brief Effective channel action on collinear two-channel storage.
 *
 * The antiunitary Theta remap exchanges the up/down channels on its own
 * (X'00 = conj(X11), X'11 = conj(X00)), so the effective permutation is the
 * XOR of the Theta swap and the U_s swap: an antiunitary operation with
 * diagonal U_s swaps the channels, while an antiunitary operation with
 * off-diagonal U_s keeps them. Incompatible stays Incompatible.
 */
CollinearChannelAction classify_collinear_action_effective(
    const std::array<std::complex<double>, 4> &U,
    bool antiunitary,
    double tol);

} // namespace librpa_int
