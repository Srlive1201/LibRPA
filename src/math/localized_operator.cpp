#include "localized_operator.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace librpa_int::localized_operator
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
            throw std::runtime_error("Localized operator contains non-finite matrix data");
        result = std::max(result, std::abs(matrix.ptr()[i]));
    }
    return result;
}

cplxdb phase(const Vector3_Order<double>& k, const Vector3_Order<int>& r)
{
    return std::polar(1.0, two_pi * (k.x * r.x + k.y * r.y + k.z * r.z));
}

Matz fourier(const LocalizedMatrix& blocks, const Vector3_Order<double>& k, const int n,
             const Matrix3& lattice, const int derivative)
{
    Matz result = zero_matrix(n);
    for (const auto& [r, matrix] : blocks)
    {
        if (matrix.nr() != n || matrix.nc() != n)
            throw std::invalid_argument("Localized operator basis dimensions differ");
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

}  // namespace librpa_int::localized_operator
