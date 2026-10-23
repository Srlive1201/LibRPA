/*!
 * @file test_symmetry_spin_kernel.cpp
 * @brief Tests for the Phase-3 four-block spinor kernel and SO(3)->SU(2) lift.
 *
 * Coverage (planning report section 13, Phase 3 exit criteria):
 * - explicit Kronecker reference in the channel-outermost (C8) layout and in
 *   the ABACUS spin-fast interleaved layout I = 2*iw + s
 * - +-U double-cover invariance
 * - time-reversal remap involution
 * - group multiplication closure, including antiunitary factors
 * - improper rotations going through the axial (proper) part
 * - so3_to_su2 / su2_to_so3 round trips including the theta = pi branches
 * - collinear action classification
 * - integration with a real dense orbital transform from the MgO fixture
 * - Identity fast path
 */
#include "../core/symmetry_spin_kernel.h"
#include "../core/atomic_basis.h"
#include "../core/pbc.h"
#include "../core/symmetry_context.h"
#include "../math/rsh.h"
#include "../math/symmetry.h"
#include "../utils/constants.h"

#include <array>
#include <cassert>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdio>
#include <map>
#include <vector>

#include "testutils.h"

namespace {

using namespace librpa_int;

using cplx = std::complex<double>;
using su2_t = std::array<cplx, 4>;
using Blocks4 = SpinorBlocks4<ComplexMatrix>;

constexpr double kTol = 1e-12;

void assert_frob_below(const double value, const double tol, const char *label)
{
    if (!(value < tol))
    {
        std::fprintf(stderr, "%s: frob diff %.6e >= tol %.6e\n", label, value, tol);
        assert(false);
    }
}

//! Deterministic LCG (Numerical Recipes constants), no external dependency.
struct SimpleLcg
{
    uint64_t state;
    explicit SimpleLcg(uint64_t seed) : state(seed) {}
    uint64_t next_u64()
    {
        state = state * 6364136223846793005ULL + 1442695040888963407ULL;
        return state;
    }
    //! Uniform in [0, 1).
    double next_uniform() { return (next_u64() >> 11) * (1.0 / 9007199254740992.0); }
    //! Uniform in [-1, 1).
    double next_symmetric() { return 2.0 * next_uniform() - 1.0; }
};

std::array<double, 3> random_unit_axis(SimpleLcg &rng)
{
    std::array<double, 3> axis;
    double norm = 0.0;
    do
    {
        axis = {rng.next_symmetric(), rng.next_symmetric(), rng.next_symmetric()};
        norm = std::sqrt(axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]);
    } while (norm < 1e-3);
    return {axis[0] / norm, axis[1] / norm, axis[2] / norm};
}

su2_t su2_from_axis_angle(const std::array<double, 3> &axis, const double theta)
{
    const double u0 = std::cos(0.5 * theta);
    const double s = std::sin(0.5 * theta);
    return {cplx(u0, -s * axis[2]), cplx(-s * axis[1], -s * axis[0]),
            cplx(s * axis[1], -s * axis[0]), cplx(u0, s * axis[2])};
}

su2_t random_su2(SimpleLcg &rng)
{
    return su2_from_axis_angle(random_unit_axis(rng), 2.0 * M_PI * rng.next_uniform());
}

su2_t su2_multiply(const su2_t &a, const su2_t &b)
{
    return {a[0] * b[0] + a[1] * b[2], a[0] * b[1] + a[1] * b[3],
            a[2] * b[0] + a[3] * b[2], a[2] * b[1] + a[3] * b[3]};
}

su2_t su2_negated(const su2_t &u)
{
    return {-u[0], -u[1], -u[2], -u[3]};
}

//! Rodrigues formula for the active rotation by theta about axis.
Matrix3 axis_angle_rotation(const std::array<double, 3> &axis_in, const double theta)
{
    const double norm = std::sqrt(axis_in[0] * axis_in[0] + axis_in[1] * axis_in[1] +
                                  axis_in[2] * axis_in[2]);
    const double nx = axis_in[0] / norm;
    const double ny = axis_in[1] / norm;
    const double nz = axis_in[2] / norm;
    const double c = std::cos(theta);
    const double s = std::sin(theta);
    const double C = 1.0 - c;
    return Matrix3(c + nx * nx * C, nx * ny * C - nz * s, nx * nz * C + ny * s,
                   ny * nx * C + nz * s, c + ny * ny * C, ny * nz * C - nx * s,
                   nz * nx * C - ny * s, nz * ny * C + nx * s, c + nz * nz * C);
}

ComplexMatrix random_unitary(SimpleLcg &rng, const int n)
{
    ComplexMatrix a(n, n, false);
    for (int i = 0; i != n; ++i)
        for (int j = 0; j != n; ++j)
            a(i, j) = cplx(rng.next_symmetric(), rng.next_symmetric());
    // classical Gram-Schmidt on columns
    for (int j = 0; j != n; ++j)
    {
        for (int k = 0; k != j; ++k)
        {
            cplx proj = 0.0;
            for (int i = 0; i != n; ++i)
                proj += std::conj(a(i, k)) * a(i, j);
            for (int i = 0; i != n; ++i)
                a(i, j) -= proj * a(i, k);
        }
        double col_norm = 0.0;
        for (int i = 0; i != n; ++i)
            col_norm += std::norm(a(i, j));
        col_norm = std::sqrt(col_norm);
        for (int i = 0; i != n; ++i)
            a(i, j) /= col_norm;
    }
    return a;
}

ComplexMatrix random_hermitian(SimpleLcg &rng, const int n)
{
    ComplexMatrix m(n, n, false);
    for (int i = 0; i != n; ++i)
        for (int j = 0; j != n; ++j)
            m(i, j) = cplx(rng.next_symmetric(), rng.next_symmetric());
    return m + transpose(m, true);
}

//! Channel-outermost (C8) split: block (a, b) holds X[a*N+i][b*N+j].
Blocks4 split4_channel_outer(const ComplexMatrix &x, const int n)
{
    Blocks4 out{ComplexMatrix(n, n, false), ComplexMatrix(n, n, false),
                ComplexMatrix(n, n, false), ComplexMatrix(n, n, false)};
    ComplexMatrix *blocks[2][2] = {{&out.b00, &out.b01}, {&out.b10, &out.b11}};
    for (int a = 0; a != 2; ++a)
        for (int b = 0; b != 2; ++b)
            for (int i = 0; i != n; ++i)
                for (int j = 0; j != n; ++j)
                    (*blocks[a][b])(i, j) = x(a * n + i, b * n + j);
    return out;
}

ComplexMatrix merge4_channel_outer(const Blocks4 &blocks)
{
    const int n = blocks.b00.nr;
    ComplexMatrix out(2 * n, 2 * n, false);
    const ComplexMatrix *b[2][2] = {{&blocks.b00, &blocks.b01}, {&blocks.b10, &blocks.b11}};
    for (int a = 0; a != 2; ++a)
        for (int c = 0; c != 2; ++c)
            for (int i = 0; i != n; ++i)
                for (int j = 0; j != n; ++j)
                    out(a * n + i, c * n + j) = (*b[a][c])(i, j);
    return out;
}

//! ABACUS spin-fast interleave: element (i, s) sits at I = 2*iw + s.
ComplexMatrix shuffle_to_spin_fast(const ComplexMatrix &x, const int n)
{
    ComplexMatrix out(2 * n, 2 * n, false);
    for (int i = 0; i != n; ++i)
        for (int j = 0; j != n; ++j)
            for (int s = 0; s != 2; ++s)
                for (int t = 0; t != 2; ++t)
                    out(2 * i + s, 2 * j + t) = x(s * n + i, t * n + j);
    return out;
}

//! Dense D_full in the channel-outermost layout: D[a*N+i][b*N+j] = U_ab A_ij.
ComplexMatrix build_dense_full(const ComplexMatrix &A, const su2_t &U)
{
    const int n = A.nr;
    ComplexMatrix d(2 * n, 2 * n, true);
    for (int a = 0; a != 2; ++a)
        for (int b = 0; b != 2; ++b)
            for (int i = 0; i != n; ++i)
                for (int j = 0; j != n; ++j)
                    d(a * n + i, b * n + j) = U[2 * a + b] * A(i, j);
    return d;
}

//! Same operator in the spin-fast layout: D[2*i+a][2*j+b] = A_ij U_ab.
ComplexMatrix build_dense_full_spin_fast(const ComplexMatrix &A, const su2_t &U)
{
    const int n = A.nr;
    ComplexMatrix d(2 * n, 2 * n, true);
    for (int a = 0; a != 2; ++a)
        for (int b = 0; b != 2; ++b)
            for (int i = 0; i != n; ++i)
                for (int j = 0; j != n; ++j)
                    d(2 * i + a, 2 * j + b) = A(i, j) * U[2 * a + b];
    return d;
}

//! J_AO = I_N orb (i sigma_y), channel-outermost; i sigma_y = [[0, 1], [-1, 0]].
ComplexMatrix build_j_ao(const int n)
{
    const double s[4] = {0.0, 1.0, -1.0, 0.0};
    ComplexMatrix j(2 * n, 2 * n, true);
    for (int a = 0; a != 2; ++a)
        for (int b = 0; b != 2; ++b)
            for (int i = 0; i != n; ++i)
                j(a * n + i, b * n + i) = s[2 * a + b];
    return j;
}

double frob_diff(const ComplexMatrix &x, const ComplexMatrix &y)
{
    return std::sqrt(abs2(x - y));
}

double blocks_frob_diff(const Blocks4 &x, const Blocks4 &y)
{
    return std::sqrt(abs2(x.b00 - y.b00) + abs2(x.b01 - y.b01) +
                     abs2(x.b10 - y.b10) + abs2(x.b11 - y.b11));
}

SymmetrySpinOperation make_op(const su2_t &U, const bool antiunitary,
                              const SymmetrySpinActionSource source =
                                  SymmetrySpinActionSource::ExplicitSpinSpace)
{
    SymmetrySpinOperation op;
    op.spatial_id = 0;
    op.spin_u = U;
    op.antiunitary = antiunitary;
    op.spin_source = source;
    return op;
}

//! Kronecker reference, planning report section 11.2.
void run_kronecker_check(const ComplexMatrix &A, const su2_t &U, const bool antiunitary,
                         SimpleLcg &rng, const char *label)
{
    const int n = A.nr;
    const ComplexMatrix X = random_hermitian(rng, 2 * n);

    const ComplexMatrix D = build_dense_full(A, U);
    ComplexMatrix X_ref = D * X * transpose(D, true);
    if (antiunitary)
    {
        const ComplexMatrix J = build_j_ao(n);
        X_ref = J * conj(X_ref) * transpose(J, true);
    }

    const auto orbit = [&A](std::size_t, const ComplexMatrix &x) {
        return A * x * transpose(A, true);
    };
    // The optional numerical U is authoritative for every public source tag.
    for (const auto source : {SymmetrySpinActionSource::Identity,
                              SymmetrySpinActionSource::ExplicitSpinSpace,
                              SymmetrySpinActionSource::DerivedFromSpatialSOC})
    {
        const auto op = make_op(U, antiunitary, source);
        const auto result = transform_spinor_bilinear(
            op, split4_channel_outer(X, n), orbit, BilinearConvention::SourceToTarget_DXDdag);
        assert_frob_below(frob_diff(merge4_channel_outer(result), X_ref), kTol, label);
    }
}

void test_kronecker_reference_random()
{
    SimpleLcg rng(20260726ULL);
    const int n = 6;
    for (int rep = 0; rep != 3; ++rep)
    {
        const ComplexMatrix A = random_unitary(rng, n);
        const su2_t U = random_su2(rng);
        run_kronecker_check(A, U, false, rng, "kronecker unitary");
        run_kronecker_check(A, U, true, rng, "kronecker antiunitary");
    }
}

//! Same reference in the ABACUS spin-fast layout I = 2*iw + s: the shuffled
//! channel-outermost result must match the directly interleaved Kronecker
//! product, guarding against a flipped direct-product order.
void test_layout_spin_fast_interleave()
{
    SimpleLcg rng(424242ULL);
    const int n = 6;
    const ComplexMatrix A = random_unitary(rng, n);
    const su2_t U = random_su2(rng);
    const ComplexMatrix X = random_hermitian(rng, 2 * n);

    const auto op = make_op(U, false);
    const auto orbit = [&A](std::size_t, const ComplexMatrix &x) {
        return A * x * transpose(A, true);
    };
    const auto result = transform_spinor_bilinear(
        op, split4_channel_outer(X, n), orbit, BilinearConvention::SourceToTarget_DXDdag);
    const ComplexMatrix merged = merge4_channel_outer(result);

    // (a) channel-outermost result shuffled to spin-fast
    const ComplexMatrix shuffled = shuffle_to_spin_fast(merged, n);
    // (b) reference built directly in the spin-fast layout
    const ComplexMatrix D_sf = build_dense_full_spin_fast(A, U);
    const ComplexMatrix X_sf = shuffle_to_spin_fast(X, n);
    const ComplexMatrix ref_sf = D_sf * X_sf * transpose(D_sf, true);

    assert_frob_below(frob_diff(shuffled, ref_sf), kTol, "spin-fast layout");
}

void test_pm_u_invariance()
{
    SimpleLcg rng(777ULL);
    const int n = 6;
    const ComplexMatrix A = random_unitary(rng, n);
    const su2_t U = random_su2(rng);
    const ComplexMatrix X = random_hermitian(rng, 2 * n);
    const auto blocks = split4_channel_outer(X, n);
    const auto orbit = [&A](std::size_t, const ComplexMatrix &x) {
        return A * x * transpose(A, true);
    };

    for (const bool anti : {false, true})
    {
        const auto plus = transform_spinor_bilinear(
            make_op(U, anti), blocks, orbit, BilinearConvention::SourceToTarget_DXDdag);
        const auto minus = transform_spinor_bilinear(
            make_op(su2_negated(U), anti), blocks, orbit,
            BilinearConvention::SourceToTarget_DXDdag);
        assert_frob_below(blocks_frob_diff(plus, minus), 1e-14, "+-U invariance");
    }
}

//! Theta remap is an involution: two antiunitary transforms with U = I and the
//! identity spatial action return the original blocks.
void test_tr_involution()
{
    SimpleLcg rng(1357ULL);
    const int n = 6;
    const ComplexMatrix X = random_hermitian(rng, 2 * n);
    const auto identity_orbit = [](std::size_t, const ComplexMatrix &x) { return x; };
    const su2_t identity_u{1.0, 0.0, 0.0, 1.0};

    const auto once = transform_spinor_bilinear(
        make_op(identity_u, true, SymmetrySpinActionSource::Identity),
        split4_channel_outer(X, n), identity_orbit,
        BilinearConvention::SourceToTarget_DXDdag);
    const auto twice = transform_spinor_bilinear(
        make_op(identity_u, true, SymmetrySpinActionSource::Identity), once,
        identity_orbit, BilinearConvention::SourceToTarget_DXDdag);

    assert_frob_below(blocks_frob_diff(twice, split4_channel_outer(X, n)), 1e-14,
                      "TR involution");
}

/*!
 * Group multiplication closure. With the antiunitary represented as
 * T(X) = J (D X D^dag)* J^dag, composing (A1, U1, eta1) then (A2, U2, eta2)
 * gives the kernel operation
 *   (conj(A2)^eta1 A1, U2 U1, eta1 XOR eta2)
 * (the spin part is always the plain SU(2) product; the conjugation lands on
 * the left factor when eta1 = 1, and J conj(U) J = -U removes it there).
 */
void test_group_multiplication()
{
    SimpleLcg rng(987654ULL);
    const int n = 6;
    const ComplexMatrix X = random_hermitian(rng, 2 * n);
    const auto blocks = split4_channel_outer(X, n);

    for (int rep = 0; rep != 3; ++rep)
    {
        const ComplexMatrix A1 = random_unitary(rng, n);
        const ComplexMatrix A2 = random_unitary(rng, n);
        const su2_t U1 = random_su2(rng);
        const su2_t U2 = random_su2(rng);

        for (const bool eta1 : {false, true})
        {
            for (const bool eta2 : {false, true})
            {
                if (eta1 && !eta2)
                {
                    // closure holds as well, but keep the spec-mandated set:
                    // all-unitary, Theta g2 o g1, and Theta g2 o Theta g1.
                    continue;
                }
                const auto orbit1 = [&A1](std::size_t, const ComplexMatrix &x) {
                    return A1 * x * transpose(A1, true);
                };
                const auto orbit2 = [&A2](std::size_t, const ComplexMatrix &x) {
                    return A2 * x * transpose(A2, true);
                };
                const auto step1 = transform_spinor_bilinear(
                    make_op(U1, eta1), blocks, orbit1,
                    BilinearConvention::SourceToTarget_DXDdag);
                const auto step2 = transform_spinor_bilinear(
                    make_op(U2, eta2), step1, orbit2,
                    BilinearConvention::SourceToTarget_DXDdag);

                const ComplexMatrix A_tot =
                    (eta1 ? conj(A2) : A2) * A1;
                const su2_t U_tot = su2_multiply(U2, U1);
                const auto orbit_tot = [&A_tot](std::size_t, const ComplexMatrix &x) {
                    return A_tot * x * transpose(A_tot, true);
                };
                const auto direct = transform_spinor_bilinear(
                    make_op(U_tot, eta1 != eta2), blocks, orbit_tot,
                    BilinearConvention::SourceToTarget_DXDdag);

                assert_frob_below(blocks_frob_diff(step2, direct), kTol,
                                  "group multiplication");
            }
        }
    }
}

//! Improper rotations must enter through the axial (proper) part det(Q) Q.
void test_improper_rotation_uses_axial_part()
{
    // mirror z: det = -1; axial part is the pi rotation about z.
    const Matrix3 mirror_z(1.0, 0.0, 0.0,
                           0.0, 1.0, 0.0,
                           0.0, 0.0, -1.0);
    assert(std::abs(mirror_z.Det() + 1.0) < 1e-14);

    const su2_t u_axial = so3_to_su2(axial_rotation_of(mirror_z));
    // pi about z: U = diag(exp(-i pi/2), exp(i pi/2)) = diag(-i, i)
    assert(fequal(u_axial[0], cplx(0.0, -1.0), cplx(kTol, 0.0)));
    assert(fequal(u_axial[3], cplx(0.0, 1.0), cplx(kTol, 0.0)));
    assert(fequal(u_axial[1], cplx(0.0, 0.0), cplx(kTol, 0.0)));
    assert(fequal(u_axial[2], cplx(0.0, 0.0), cplx(kTol, 0.0)));

    // Forgetting the proper part and feeding the improper matrix directly must
    // NOT give the same result (regression guard for risk R2).
    const su2_t u_naive = so3_to_su2(mirror_z);
    double diff = 0.0;
    for (int i = 0; i != 4; ++i)
        diff += std::norm(u_axial[i] - u_naive[i]);
    assert(std::sqrt(diff) > 1e-6);

    // A generic improper operation: inversion times C3 about (1,1,1).
    const std::array<double, 3> axis111{1.0, 1.0, 1.0};
    const Matrix3 improper = axis_angle_rotation(axis111, 2.0 * M_PI / 3.0) * (-1.0);
    assert(std::abs(improper.Det() + 1.0) < 1e-12);
    const Matrix3 proper = axial_rotation_of(improper);
    assert(std::abs(proper.Det() - 1.0) < 1e-12);
    assert(is_same_matrix(su2_to_so3(so3_to_su2(proper)), proper, kTol));
}

//! so3_to_su2 / su2_to_so3 round trips, covering every theta = pi branch.
void test_so3_su2_round_trip()
{
    const double third_pi = M_PI / 3.0;
    const std::array<std::array<double, 3>, 6> axes{{
        {{1.0, 0.0, 0.0}},
        {{0.0, 1.0, 0.0}},
        {{0.0, 0.0, 1.0}},
        {{1.0, 1.0, 1.0}},
        {{0.2, 1.0, 0.5}},
        {{0.3, -0.7, 0.1}},
    }};
    const std::array<double, 4> angles{{third_pi, M_PI, 1.234, 1e-13}};

    for (const auto &axis : axes)
    {
        for (const double theta : angles)
        {
            const Matrix3 R = axis_angle_rotation(axis, theta);
            const su2_t U = so3_to_su2(R);
            // U must be unitary
            const su2_t UUdag = su2_multiply(U, {std::conj(U[0]), std::conj(U[2]),
                                                 std::conj(U[1]), std::conj(U[3])});
            assert(fequal(UUdag[0], cplx(1.0, 0.0), cplx(1e-12, 0.0)));
            assert(fequal(UUdag[3], cplx(1.0, 0.0), cplx(1e-12, 0.0)));
            assert(fequal(UUdag[1], cplx(0.0, 0.0), cplx(1e-12, 0.0)));
            assert(fequal(UUdag[2], cplx(0.0, 0.0), cplx(1e-12, 0.0)));
            assert(is_same_matrix(su2_to_so3(U), R, kTol));
        }
    }

    // theta = pi branches: dominant x, y, z and a tie (1,1,1).
    assert(is_same_matrix(su2_to_so3(so3_to_su2(axis_angle_rotation({{1.0, 0.0, 0.0}}, M_PI))),
                          axis_angle_rotation({{1.0, 0.0, 0.0}}, M_PI), kTol));
    assert(is_same_matrix(su2_to_so3(so3_to_su2(axis_angle_rotation({{0.0, 1.0, 0.0}}, M_PI))),
                          axis_angle_rotation({{0.0, 1.0, 0.0}}, M_PI), kTol));
    assert(is_same_matrix(su2_to_so3(so3_to_su2(axis_angle_rotation({{0.0, 0.0, 1.0}}, M_PI))),
                          axis_angle_rotation({{0.0, 0.0, 1.0}}, M_PI), kTol));

    // identity maps to identity
    const su2_t u_id = so3_to_su2(Matrix3::IDENTITY);
    assert(fequal(u_id[0], cplx(1.0, 0.0), cplx(kTol, 0.0)));
    assert(fequal(u_id[3], cplx(1.0, 0.0), cplx(kTol, 0.0)));
    assert(fequal(u_id[1], cplx(0.0, 0.0), cplx(kTol, 0.0)));
    assert(fequal(u_id[2], cplx(0.0, 0.0), cplx(kTol, 0.0)));
}

void test_classify_collinear_action()
{
    const double tol = 1e-10;
    // diagonal up to phase -> Keep
    const su2_t u_keep{cplx(std::cos(0.7), std::sin(0.7)), 0.0,
                       0.0, cplx(std::cos(0.7), -std::sin(0.7))};
    assert(classify_collinear_action(u_keep, false, tol) == CollinearChannelAction::Keep);
    assert(classify_collinear_action(u_keep, true, tol) == CollinearChannelAction::Keep);
    // i sigma_y -> Swap
    const su2_t u_swap{0.0, 1.0, -1.0, 0.0};
    assert(classify_collinear_action(u_swap, false, tol) == CollinearChannelAction::Swap);
    assert(classify_collinear_action(u_swap, true, tol) == CollinearChannelAction::Swap);
    // +-U must classify identically
    assert(classify_collinear_action(su2_negated(u_swap), false, tol) ==
           CollinearChannelAction::Swap);
    // rotation about x by pi/3 -> Incompatible
    const su2_t u_mix = su2_from_axis_angle({{1.0, 0.0, 0.0}}, M_PI / 3.0);
    assert(classify_collinear_action(u_mix, false, tol) ==
           CollinearChannelAction::Incompatible);
    assert(classify_collinear_action(u_mix, true, tol) ==
           CollinearChannelAction::Incompatible);
}

void add_mgo_fractional_symmetry_operations(SymmetryContext &ctx)
{
    const std::array<std::array<int, 3>, 6> permutations{{
        {{0, 1, 2}},
        {{0, 2, 1}},
        {{1, 0, 2}},
        {{1, 2, 0}},
        {{2, 0, 1}},
        {{2, 1, 0}},
    }};
    const std::array<int, 2> signs{{-1, 1}};

    for (const auto &permutation : permutations)
    {
        for (const int sx : signs)
        {
            for (const int sy : signs)
            {
                for (const int sz : signs)
                {
                    const std::array<int, 3> sign{{sx, sy, sz}};
                    std::array<int, 9> rotation{};
                    for (int col = 0; col != 3; ++col)
                    {
                        rotation[3 * permutation[col] + col] = sign[col];
                    }
                    const Matrix3 col_cartesian_rotation(
                        rotation[0], rotation[1], rotation[2],
                        rotation[3], rotation[4], rotation[5],
                        rotation[6], rotation[7], rotation[8]);
                    const Matrix3 row_cartesian_rotation = col_cartesian_rotation.Transpose();
                    SymmetryOperation op;
                    op.rotation = ctx.lattice_available
                                      ? ctx.lattice_vectors * row_cartesian_rotation *
                                            ctx.lattice_vectors.Inverse()
                                      : row_cartesian_rotation;
                    op.translation = {0.0, 0.0, 0.0};
                    op.use_row_convention = true;
                    ctx.rspace_operations.push_back(op);
                }
            }
        }
    }
    assert(ctx.rspace_operations.size() == 48);
}

/*!
 * Integration with the production orbital transform: take the dense A of one
 * real MgO star member (build_symmetry_kspace_operator_transform_matrix) and
 * the SU(2) lifted from that member's Cartesian rotation, then run the
 * Kronecker reference against the kernel fed with A X A^dagger as the orbit
 * callback. This validates the kernel/orbit interface, not new physics.
 */
void test_mgo_fixture_orbital_transform_matches_kronecker()
{
    PeriodicBoundaryData pbc;
    pbc.set_latvec({0.0, 0.5, 0.5,
                    0.5, 0.0, 0.5,
                    0.5, 0.5, 0.0});

    std::vector<double> kvecs;
    for (const auto &kfrac : build_uniform_kmesh_frac({3, 3, 3}))
    {
        const auto kvec = kfrac * pbc.G;
        kvecs.push_back(kvec.x * TWO_PI);
        kvecs.push_back(kvec.y * TWO_PI);
        kvecs.push_back(kvec.z * TWO_PI);
    }
    pbc.set_kgrids_kvec(3, 3, 3, kvecs);

    SymmetryContext ctx;
    ctx.set_crystal_structure(pbc.latvec, pbc.G,
                              {{0, 0}, {1, 1}},
                              {{0, {0.0, 0.0, 0.0}}, {1, {0.5, 0.5, 0.5}}});
    add_mgo_fractional_symmetry_operations(ctx);
    ctx.set_available();
    ctx.build_periodic_mappings(pbc, pbc.Rlist);
    ctx.build_rsh_rotations({-1,
                             0,
                             LIBRPA_ANGULAR_ORDER_NATURAL,
                             LIBRPA_RSH_COEFF_1_M,
                             LIBRPA_RSH_COEFF_1_M},
                            1);
    ctx.build_kstar_member_rotations(1);

    AtomicBasis basis(std::vector<std::size_t>{4, 4});
    basis.set_l_shells({{0, 1}, {0, 1}});
    const auto layouts = basis.build_species_basis_layouts(ctx.atom_to_type);
    const std::map<atom_t, size_t> atom_nw{{0, 4}, {1, 4}};

    // pick a unitary member with a non-identity spatial rotation
    const SymmetryKStarMember *picked = nullptr;
    const SymmetryKStar *picked_star = nullptr;
    for (const auto &star : ctx.kstars)
    {
        for (const auto &member : star.members)
        {
            if (!member.time_reversal &&
                !ctx.rspace_operations[member.spatial_isym].is_identity_rotation())
            {
                picked = &member;
                picked_star = &star;
                break;
            }
        }
        if (picked != nullptr)
        {
            break;
        }
    }
    assert(picked != nullptr);

    const ComplexMatrix A = build_symmetry_kspace_operator_transform_matrix(
        ctx, layouts, *picked, atom_nw, picked_star->k_ibz, false, &picked->k_bz);
    assert(A.nr == 8 && A.nc == 8);

    const Matrix3 cartesian = fractional_rotation_to_cartesian(
        ctx.rspace_operations[picked->spatial_isym], ctx.lattice_vectors);
    const Matrix3 proper = axial_rotation_of(cartesian);
    const su2_t U = so3_to_su2(proper);
    // sanity: the lift round-trips the member's proper rotation
    assert(is_same_matrix(su2_to_so3(U), proper, kTol));

    SimpleLcg rng(20260727ULL);
    run_kronecker_check(A, U, false, rng, "MgO fixture unitary");
    run_kronecker_check(A, U, true, rng, "MgO fixture antiunitary");
}

//! Identity spin source with eta = 0 must return the plain scalar-orbit
//! result block by block (fast path preserving the old scalar behavior).
void test_identity_fast_path_matches_scalar_orbit()
{
    SimpleLcg rng(555000ULL);
    const int n = 6;
    const ComplexMatrix A = random_unitary(rng, n);
    const ComplexMatrix X = random_hermitian(rng, 2 * n);
    const auto blocks = split4_channel_outer(X, n);
    const auto orbit = [&A](std::size_t, const ComplexMatrix &x) {
        return A * x * transpose(A, true);
    };

    const su2_t identity_u{1.0, 0.0, 0.0, 1.0};
    const auto fast = transform_spinor_bilinear(
        make_op(identity_u, false, SymmetrySpinActionSource::Identity), blocks, orbit,
        BilinearConvention::SourceToTarget_DXDdag);

    const Blocks4 expected{orbit(0, blocks.b00), orbit(0, blocks.b01),
                           orbit(0, blocks.b10), orbit(0, blocks.b11)};
    assert_frob_below(blocks_frob_diff(fast, expected), 1e-15, "identity fast path");

    // an explicit U = I must agree numerically with the fast path
    const auto explicit_identity = transform_spinor_bilinear(
        make_op(identity_u, false, SymmetrySpinActionSource::ExplicitSpinSpace),
        blocks, orbit, BilinearConvention::SourceToTarget_DXDdag);
    assert_frob_below(blocks_frob_diff(fast, explicit_identity), kTol,
                      "explicit identity spin");
}

} // namespace

int main()
{
    test_kronecker_reference_random();
    test_layout_spin_fast_interleave();
    test_pm_u_invariance();
    test_tr_involution();
    test_group_multiplication();
    test_improper_rotation_uses_axial_part();
    test_so3_su2_round_trip();
    test_classify_collinear_action();
    test_mgo_fixture_orbital_transform_matches_kronecker();
    test_identity_fast_path_matches_scalar_orbit();
    return 0;
}
