#include "headwing_velocity.h"

#include <algorithm>
#include <cmath>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>

#include "../../utils/constants.h"

namespace librpa_int::qsgw
{
namespace
{
constexpr double two_pi = 6.283185307179586476925286766559;

Matz zero_matrix(const int n)
{
    Matz result(n, n, MAJOR::ROW);
    std::fill(result.ptr(), result.ptr() + result.size(), cplxdb{});
    return result;
}

double maximum_entry(const Matz& matrix)
{
    double result = 0.0;
    for (std::size_t i = 0; i < matrix.size(); ++i)
    {
        if (!std::isfinite(matrix.ptr()[i].real()) || !std::isfinite(matrix.ptr()[i].imag()))
            throw std::runtime_error("QSGW live velocity contains non-finite matrix data");
        result = std::max(result, std::abs(matrix.ptr()[i]));
    }
    return result;
}

cplxdb phase(const Vector3_Order<double>& k, const Vector3_Order<int>& r)
{
    return std::polar(1.0, two_pi * (k.x * r.x + k.y * r.y + k.z * r.z));
}

Matz fourier(const LocalizedMatrix& blocks, const Vector3_Order<double>& k, const int n,
             const Matrix3& lattice, const int derivative = -1)
{
    Matz result = zero_matrix(n);
    for (const auto& [r, matrix] : blocks)
    {
        if (matrix.nr() != n || matrix.nc() != n)
            throw std::invalid_argument("QSGW live velocity localized basis dimensions differ");
        cplxdb weight = phase(k, r);
        if (derivative >= 0)
        {
            const std::array<double, 3> distance{
                r.x * lattice.e11 + r.y * lattice.e21 + r.z * lattice.e31,
                r.x * lattice.e12 + r.y * lattice.e22 + r.z * lattice.e32,
                r.x * lattice.e13 + r.y * lattice.e23 + r.z * lattice.e33};
            weight *= cplxdb(0.0, distance.at(derivative));
        }
        maximum_entry(matrix);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) result(i, j) += weight * matrix(i, j);
    }
    return result;
}

// S/r use the audited PyATB upper-triangle storage and completion convention.
Matz complete_hermitian_upper(Matz upper)
{
    for (int i = 0; i < upper.nr(); ++i)
    {
        upper(i, i) = upper(i, i).real();
        for (int j = i + 1; j < upper.nc(); ++j) upper(j, i) = std::conj(upper(i, j));
    }
    return upper;
}

struct FourierVelocityModel
{
    Matz hamiltonian;
    Matz overlap;
    std::array<Matz, 3> hamiltonian_derivative;
    std::array<Matz, 3> overlap_derivative;
    std::array<Matz, 3> position;
};

FourierVelocityModel fourier_velocity_model(const LocalizedMatrix& hamiltonian,
                                            const HeadwingAoOperators& operators,
                                            const Vector3_Order<double>& k, const Matrix3& lattice)
{
    const int n = operators.dimension;
    if (n <= 0 || hamiltonian.empty() || operators.overlap.empty())
        throw std::invalid_argument("Localized velocity requires a nonempty square H/S model");
    const std::array<double, 12> coordinates{k.x,         k.y,         k.z,         lattice.e11,
                                             lattice.e12, lattice.e13, lattice.e21, lattice.e22,
                                             lattice.e23, lattice.e31, lattice.e32, lattice.e33};
    for (const double value : coordinates)
        if (!std::isfinite(value))
            throw std::invalid_argument("Localized velocity requires finite k and lattice vectors");
    FourierVelocityModel result;
    result.hamiltonian = fourier(hamiltonian, k, n, lattice);
    result.overlap = complete_hermitian_upper(fourier(operators.overlap, k, n, lattice));
    for (int a = 0; a < 3; ++a)
    {
        result.hamiltonian_derivative[a] = fourier(hamiltonian, k, n, lattice, a);
        result.overlap_derivative[a] =
            complete_hermitian_upper(fourier(operators.overlap, k, n, lattice, a));
        Matz r = fourier(operators.position[a], k, n, lattice);
        const Matz raw_r = r.copy();
        for (int i = 0; i < n; ++i)
            for (int j = 0; j <= i; ++j)
                r(i, j) =
                    std::conj(raw_r(j, i) + cplxdb(0.0, 1.0) * result.overlap_derivative[a](j, i));
        result.position[a] = std::move(r);
    }
    return result;
}

LocalizedVelocityResult project_velocity(const FourierVelocityModel& model,
                                         const std::vector<double>& eigenvalues,
                                         const Matz& coefficients)
{
    const int n = model.hamiltonian.nr();
    if (static_cast<int>(eigenvalues.size()) != n || coefficients.nr() != n ||
        coefficients.nc() != n)
        throw std::invalid_argument("Localized velocity requires complete matching E/C dimensions");
    for (const double energy : eigenvalues)
        if (!std::isfinite(energy))
            throw std::invalid_argument("Localized velocity contains a non-finite eigenvalue");
    maximum_entry(coefficients);
    LocalizedVelocityResult result;
    result.eigenvalues_ha = eigenvalues;
    result.ket_coefficients = coefficients.copy();
    result.ket_coefficients.swap_to_row_major();
    result.overlap = model.overlap.copy();
    const Matz& u = result.ket_coefficients;
    const Matz ud = transpose(u, true);
    Matz metric = ud * model.overlap * u;
    for (int i = 0; i < n; ++i) metric(i, i) -= 1.0;
    result.metric_residual = maximum_entry(metric);
    Matz state_h = ud * model.hamiltonian * u;
    for (int i = 0; i < n; ++i) state_h(i, i) -= eigenvalues[i];
    result.hamiltonian_frame_residual_ha = maximum_entry(state_h);
    result.hamiltonian_hermiticity_residual_ha =
        maximum_entry(model.hamiltonian - transpose(model.hamiltonian, true));
    for (int a = 0; a < 3; ++a)
    {
        const Matz dh = ud * model.hamiltonian_derivative[a] * u;
        const Matz ds = ud * model.overlap_derivative[a] * u;
        const Matz r = ud * model.position[a] * u;
        Matz v = zero_matrix(n);
        for (int i = 0; i < n; ++i)
            for (int j = i; j < n; ++j)
            {
                const cplxdb value = dh(i, j) - eigenvalues[i] * ds(i, j) +
                                     cplxdb(0.0, eigenvalues[i] - eigenvalues[j]) * r(i, j);
                if (!std::isfinite(value.real()) || !std::isfinite(value.imag()))
                    throw std::runtime_error("Localized velocity contains non-finite values");
                if (i == j)
                    result.velocity_upper_diagonal_imaginary_ha_bohr = std::max(
                        result.velocity_upper_diagonal_imaginary_ha_bohr, std::abs(value.imag()));
                v(i, j) = i == j ? cplxdb(value.real(), 0.0) : value;
                v(j, i) = std::conj(v(i, j));
            }
        result.velocity_ha_bohr[a] = std::move(v);
    }
    // Retain the existing live-path gates. The Hermiticity residual is reported
    // independently; an inconsistent full H also enters the eigenframe check.
    if (result.metric_residual > 1e-6 || result.hamiltonian_frame_residual_ha > 1e-4 / HA2EV ||
        result.velocity_upper_diagonal_imaginary_ha_bohr > 1e-8)
    {
        std::ostringstream message;
        message << std::scientific
                << "Localized velocity check failed: metric=" << result.metric_residual
                << " frame_Ha=" << result.hamiltonian_frame_residual_ha
                << " diagonal_imag_HaBohr=" << result.velocity_upper_diagonal_imaginary_ha_bohr;
        throw std::runtime_error(message.str());
    }
    return result;
}

void diagonalize_localized_model(const FourierVelocityModel& model,
                                 std::vector<double>& eigenvalues, Matz& coefficients)
{
    const int n = model.hamiltonian.nr();
    coefficients = model.hamiltonian.copy();
    coefficients.swap_to_col_major();
    Matz overlap = model.overlap.copy();
    overlap.swap_to_col_major();
    eigenvalues.resize(n);
    const int itype = 1;
    const char jobz = 'V', uplo = 'U';
    int info = 0;
    int lwork = -1;
    cplxdb query;
    std::vector<double> real_work(std::max(1, 3 * n - 2));
    // Use the existing LAPACK ABI on column-major copies. The older
    // ComplexMatrix zhegv connector passes info by value, losing failure status.
    zhegv_(&itype, &jobz, &uplo, &n, coefficients.ptr(), &n, overlap.ptr(), &n, eigenvalues.data(),
           &query, &lwork, real_work.data(), &info);
    if (info != 0)
        throw std::runtime_error("Localized generalized eigensolve workspace query failed");
    lwork = std::max(std::max(1, 2 * n - 1), static_cast<int>(std::ceil(query.real())));
    std::vector<cplxdb> work(lwork);
    zhegv_(&itype, &jobz, &uplo, &n, coefficients.ptr(), &n, overlap.ptr(), &n, eigenvalues.data(),
           work.data(), &lwork, real_work.data(), &info);
    if (info != 0)
    {
        std::ostringstream message;
        message << "Localized generalized eigensolve failed (LAPACK info=" << info << ")";
        if (info > n) message << ": overlap is not positive definite";
        throw std::runtime_error(message.str());
    }
    coefficients.swap_to_row_major();
}

}  // namespace

LocalizedVelocityResult solve_localized_velocity(const LocalizedMatrix& hamiltonian,
                                                 const HeadwingAoOperators& operators,
                                                 const Vector3_Order<double>& kpoint_fractional,
                                                 const Matrix3& lattice_bohr)
{
    const auto model =
        fourier_velocity_model(hamiltonian, operators, kpoint_fractional, lattice_bohr);
    std::vector<double> eigenvalues;
    Matz coefficients;
    diagonalize_localized_model(model, eigenvalues, coefficients);
    return project_velocity(model, eigenvalues, coefficients);
}

LocalizedVelocityResult compute_localized_velocity_in_frame(
    const LocalizedMatrix& hamiltonian, const HeadwingAoOperators& operators,
    const Vector3_Order<double>& kpoint_fractional, const Matrix3& lattice_bohr,
    const std::vector<double>& eigenvalues_ha, const Matz& ket_coefficients)
{
    return project_velocity(
        fourier_velocity_model(hamiltonian, operators, kpoint_fractional, lattice_bohr),
        eigenvalues_ha, ket_coefficients);
}

}  // namespace librpa_int::qsgw
