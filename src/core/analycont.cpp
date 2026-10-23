/*!
 * @author    Min-Ye Zhang
 * @date      2024-04-25
 */
#include "analycont.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "../utils/constants.h"
// #include <iostream>

#include "../math/complexmatrix.h"
// #include "../io/stl_io_helper.h"

namespace librpa_int
{

namespace
{

bool finite_complex(const cplxdb &value)
{
    return std::isfinite(value.real()) && std::isfinite(value.imag());
}

using cplxld = std::complex<long double>;

std::pair<cplxdb, cplxdb> evaluate_continuants(const std::vector<cplxdb> &coefficients,
                                               const std::vector<cplxdb> &nodes, const cplxdb &x)
{
    // Forward convergents remain regular at an internal nested-fraction pole.
    const cplxld z(x.real(), x.imag());
    cplxld previous_a = 0.0L, a(coefficients[0].real(), coefficients[0].imag());
    cplxld previous_b = 1.0L, b = 1.0L;
    cplxld previous_da = 0.0L, da = 0.0L, previous_db = 0.0L, db = 0.0L;
    for (int i = 1; i < static_cast<int>(coefficients.size()); ++i)
    {
        const cplxld coefficient(coefficients[i].real(), coefficients[i].imag());
        const cplxld term = coefficient * (z - cplxld(nodes[i - 1].real(), nodes[i - 1].imag()));
        const cplxld next_a = a + term * previous_a;
        const cplxld next_b = b + term * previous_b;
        const cplxld next_da = da + coefficient * previous_a + term * previous_da;
        const cplxld next_db = db + coefficient * previous_b + term * previous_db;
        previous_a = a;
        previous_b = b;
        previous_da = da;
        previous_db = db;
        a = next_a;
        b = next_b;
        da = next_da;
        db = next_db;
        const long double scale =
            std::max({std::abs(a), std::abs(b), std::abs(previous_a), std::abs(previous_b),
                      std::abs(da), std::abs(db), std::abs(previous_da), std::abs(previous_db)});
        if (scale > 0.0L)
            for (auto *value :
                 {&a, &b, &previous_a, &previous_b, &da, &db, &previous_da, &previous_db})
                *value /= scale;
    }
    const cplxld value = a / b;
    const cplxld derivative = (da - value * db) / b;
    return {{static_cast<double>(value.real()), static_cast<double>(value.imag())},
            {static_cast<double>(derivative.real()), static_cast<double>(derivative.imag())}};
}

// Complete pivoting avoids normal equations for this small homogeneous system.
std::vector<std::vector<cplxld>> rational_nullspace(std::vector<std::vector<cplxld>> matrix)
{
    const int rows = matrix.size();
    const int columns = matrix.front().size();
    std::vector<int> permutation(columns);
    std::vector<long double> scales(columns, 1.0L);
    for (auto &row : matrix)
    {
        long double scale = 0.0L;
        for (const auto value : row) scale = std::max(scale, std::abs(value));
        if (scale > 0.0L)
            for (auto &value : row) value /= scale;
    }
    for (int column = 0; column < columns; ++column)
    {
        permutation[column] = column;
        long double scale = 0.0L;
        for (const auto &row : matrix) scale = std::max(scale, std::abs(row[column]));
        if (scale > 0.0L)
        {
            scales[column] = scale;
            for (auto &row : matrix) row[column] /= scale;
        }
    }

    int rank = 0;
    for (int pivot = 0; pivot < rows; ++pivot)
    {
        int pivot_row = pivot, pivot_column = pivot;
        long double size = 0.0L;
        for (int row = pivot; row < rows; ++row)
            for (int column = pivot; column < columns; ++column)
                if (std::abs(matrix[row][column]) > size)
                {
                    size = std::abs(matrix[row][column]);
                    pivot_row = row;
                    pivot_column = column;
                }
        if (size <= 64.0L * std::numeric_limits<long double>::epsilon()) break;
        std::swap(matrix[pivot], matrix[pivot_row]);
        for (auto &row : matrix) std::swap(row[pivot], row[pivot_column]);
        std::swap(permutation[pivot], permutation[pivot_column]);
        std::swap(scales[pivot], scales[pivot_column]);
        for (int row = pivot + 1; row < rows; ++row)
        {
            const cplxld factor = matrix[row][pivot] / matrix[pivot][pivot];
            matrix[row][pivot] = 0.0L;
            for (int column = pivot + 1; column < columns; ++column)
                matrix[row][column] -= factor * matrix[pivot][column];
        }
        ++rank;
    }

    std::vector<std::vector<cplxld>> basis;
    for (int free_column = rank; free_column < columns; ++free_column)
    {
        std::vector<cplxld> vector(columns, 0.0L), weights(columns);
        vector[free_column] = 1.0L;
        for (int row = rank - 1; row >= 0; --row)
        {
            cplxld sum = 0.0L;
            for (int column = row + 1; column < columns; ++column)
                sum += matrix[row][column] * vector[column];
            vector[row] = -sum / matrix[row][row];
        }
        for (int column = 0; column < columns; ++column)
            weights[permutation[column]] = vector[column] / scales[column];
        basis.push_back(std::move(weights));
    }
    // A lower-degree interpolant can have several homogeneous representations.
    // Their sum can remove spurious support zeros of individual basis vectors.
    if (basis.size() > 1)
    {
        std::vector<cplxld> weights(columns, 0.0L);
        for (const auto &vector : basis)
            for (int column = 0; column < columns; ++column) weights[column] += vector[column];
        basis.push_back(std::move(weights));
    }
    return basis;
}

}  // namespace

AnalyContPade::AnalyContPade(int n_pars_in, const std::vector<cplxdb> &xs, const std::vector<cplxdb> &data)
    : n_pars(n_pars_in)
{
    const int n_data = data.size();

    // Use all data if requested number of parameters is invalid
    if (n_pars < 0) n_pars = n_data;
    if (xs.size() != data.size() || n_data == 0 || n_pars == 0)
        throw std::invalid_argument(
            "Padé requires nonempty matching data and a nonzero parameter count");

    if (n_data <= n_pars)
    {
        // Use all data points
        n_pars = n_data;
        source_xs = xs;
        source_data = data;
    }
    else
    {
        // Select the data points evenly, when number of parameters are fewer than data points
        source_xs.resize(n_pars);
        source_data.resize(n_pars);
        int step = n_pars == 1 ? 0 : n_data / (n_pars - 1);
        for (int ipar = 0; ipar < n_pars - 1; ipar++)
        {
            source_xs[ipar] = xs[ipar * step];
            source_data[ipar] = data[ipar * step];
        }
        source_xs[n_pars - 1] = xs[n_pars == 1 ? 0 : n_data - 1];
        source_data[n_pars - 1] = data[n_pars == 1 ? 0 : n_data - 1];
    }

    // Band exclusion can make the correlation self-energy identically zero.
    if (!source_data.empty() && std::all_of(source_data.begin(), source_data.end(),
                                          [](const cplxdb &value) { return value == cplxdb(0.0, 0.0); }))
    {
        n_pars = 1;
        par_y = {cplxdb(0.0, 0.0)};
        return;
    }

    // Calculate the continuation coefficients, using Thiel's reciprocal difference method
    ComplexMatrix g(n_pars, n_pars);
    double data_scale = 0.0;
    for (int i_par = 0; i_par < n_pars; i_par++)
    {
        if (!finite_complex(source_xs[i_par]) || !finite_complex(source_data[i_par]))
            throw std::invalid_argument("Padé source data must be finite");
        g(i_par, 0) = source_data[i_par];
        data_scale = std::max(data_scale, std::abs(source_data[i_par]));
    }
    par_y.resize(n_pars);
    par_y[0] = g(0, 0);
    const double epsilon = std::numeric_limits<double>::epsilon();
    const auto lower_order_fits = [&](int n_terms)
    {
        for (int i = 0; i != static_cast<int>(source_xs.size()); ++i)
        {
            cplxdb tmp{1.0, 0.0};
            for (int p = n_terms - 1; p > 0; --p)
                tmp = 1.0 + par_y[p] * (source_xs[i] - source_xs[p - 1]) / tmp;
            const cplxdb value = par_y[0] / tmp;
            if (!finite_complex(value) ||
                std::abs(value - source_data[i]) > 64.0 * epsilon * data_scale)
                return false;
        }
        return true;
    };
    for (int i_par = 1; i_par < n_pars; i_par++)
    {
        const cplxdb pivot = g(i_par - 1, i_par - 1);
        double column_scale = std::abs(pivot);
        double deviation = 0.0;
        bool zero_pivot = pivot == cplxdb{0.0, 0.0};
        for (int i = i_par; i < n_pars; ++i)
        {
            const cplxdb previous = g(i, i_par - 1);
            column_scale = std::max(column_scale, std::abs(previous));
            deviation = std::max(deviation, std::abs(previous - pivot));
            zero_pivot = zero_pivot || previous == cplxdb{0.0, 0.0};
        }
        // Terminate a lower-order fraction only after checking all selected data.
        if ((zero_pivot || deviation <= 32.0 * epsilon * column_scale) && lower_order_fits(i_par))
        {
            n_pars = i_par;
            par_y.resize(n_pars);
            return;
        }
        if (zero_pivot)
        {
            if (initialize_rational()) return;
            throw std::runtime_error(
                "Padé reciprocal differences have a zero pivot without a valid lower-order fit");
        }
        for (int i = i_par; i < n_pars; i++)
        {
            const cplxdb denominator = (source_xs[i] - source_xs[i_par - 1]) * g(i, i_par - 1);
            if (denominator == cplxdb{0.0, 0.0} || !finite_complex(denominator))
            {
                if (initialize_rational()) return;
                throw std::runtime_error("Padé reciprocal-difference denominator is singular");
            }
            g(i, i_par) = (pivot - g(i, i_par - 1)) / denominator;
            if (!finite_complex(g(i, i_par)))
            {
                if (initialize_rational()) return;
                throw std::runtime_error(
                    "Padé reciprocal differences produced non-finite coefficients");
            }
        }
        par_y[i_par] = g(i_par, i_par);
    }
}

bool AnalyContPade::initialize_rational()
{
    // The same type [(N-1)/2, N/2] interpolant in barycentric coordinates.
    // Unlike the multiplicative fraction, this representation permits zeros.
    const int count = source_xs.size();
    const int supports = count / 2 + 1;
    std::vector<cplxld> xs, values;
    rational_y_scale = 0.0L;
    for (int i = 0; i < count; ++i)
    {
        xs.emplace_back(source_xs[i].real(), source_xs[i].imag());
        values.emplace_back(source_data[i].real(), source_data[i].imag());
        rational_x_scale = std::max(rational_x_scale, std::abs(xs.back()));
        rational_y_scale = std::max(rational_y_scale, std::abs(values.back()));
    }
    if (count < 2 || rational_y_scale == 0.0L) return false;
    for (auto &x : xs) x /= rational_x_scale;
    for (auto &value : values) value /= rational_y_scale;
    for (int i = 0; i < count; ++i)
        for (int j = 0; j < i; ++j)
            if (xs[i] == xs[j]) return false;

    for (int shift = 0; shift < count; ++shift)
    {
        rational_supports.clear();
        std::vector<bool> is_support(count, false);
        for (int j = 0; j < supports; ++j)
        {
            const int index = (j * (count - 1) / (supports - 1) + shift) % count;
            rational_supports.push_back(index);
            is_support[index] = true;
        }
        std::vector<std::vector<cplxld>> equations;
        for (int i = 0; i < count; ++i)
        {
            if (is_support[i]) continue;
            std::vector<cplxld> row;
            for (const int j : rational_supports)
                row.push_back((values[i] - values[j]) / (xs[i] - xs[j]));
            equations.push_back(std::move(row));
        }
        if (count % 2 == 0)
        {
            // Cancel the leading numerator coefficient: degree <= N/2-1.
            std::vector<cplxld> row;
            for (const int j : rational_supports) row.push_back(values[j]);
            equations.push_back(std::move(row));
        }
        rational_anchor = -1;
        if (count % 2 == 0)
        {
            long double largest = 0.0L;
            for (int j = 0; j < supports; ++j)
                if (std::abs(values[rational_supports[j]]) > largest)
                {
                    largest = std::abs(values[rational_supports[j]]);
                    rational_anchor = j;
                }
            if (rational_anchor < 0) continue;
        }
        const auto candidates = rational_nullspace(std::move(equations));
        for (auto weights : candidates)
        {
            long double scale = 0.0L;
            for (const auto weight : weights) scale = std::max(scale, std::abs(weight));
            if (scale == 0.0L || !std::isfinite(scale)) continue;
            for (auto &weight : weights) weight /= scale;
            if (rational_anchor >= 0)
            {
                cplxld leading = 0.0L;
                for (int j = 0; j < supports; ++j)
                    if (j != rational_anchor) leading += weights[j] * values[rational_supports[j]];
                weights[rational_anchor] = -leading / values[rational_supports[rational_anchor]];
            }
            rational_weights = std::move(weights);
            bool fits = true;
            for (int i = 0; i < count; ++i)
            {
                const auto value = evaluate_rational(source_xs[i]).first;
                if (!finite_complex(value) ||
                    std::abs(cplxld(value.real(), value.imag()) / rational_y_scale - values[i]) >
                        64.0L * std::numeric_limits<double>::epsilon())
                {
                    fits = false;
                    break;
                }
            }
            if (fits) return true;
        }
    }
    rational_weights.clear();
    return false;
}

std::pair<cplxdb, cplxdb> AnalyContPade::evaluate_rational(const cplxdb &x) const
{
    const cplxld z = cplxld(x.real(), x.imag()) / rational_x_scale;
    cplxld numerator = 0.0L, denominator = 0.0L;
    cplxld d_numerator = 0.0L, d_denominator = 0.0L;
    int hit = -1;
    cplxld anchor = 0.0L;
    if (rational_anchor >= 0)
    {
        const auto point = source_xs[rational_supports[rational_anchor]];
        anchor = cplxld(point.real(), point.imag()) / rational_x_scale;
        denominator = rational_weights[rational_anchor];
    }
    for (int j = 0; j < static_cast<int>(rational_supports.size()); ++j)
    {
        const int index = rational_supports[j];
        if (j == rational_anchor) continue;
        if (x == source_xs[index])
        {
            hit = j;
            continue;
        }
        const cplxld node =
            cplxld(source_xs[index].real(), source_xs[index].imag()) / rational_x_scale;
        const cplxld delta = z - node;
        const cplxld value =
            cplxld(source_data[index].real(), source_data[index].imag()) / rational_y_scale;
        const cplxld term = rational_weights[j] / delta;
        if (rational_anchor >= 0)
        {
            // Explicitly cancel the even-N leading numerator coefficient.
            const cplxld distance = node - anchor;
            numerator += term * value * distance;
            denominator += term * (z - anchor);
            d_numerator -= term * value * distance / delta;
            d_denominator -= term * distance / delta;
        }
        else
        {
            numerator += term * value;
            denominator += term;
            d_numerator -= term * value / delta;
            d_denominator -= term / delta;
        }
    }
    cplxld value, derivative;
    if (hit >= 0 && rational_weights[hit] != cplxld{0.0L, 0.0L})
    {
        const int index = rational_supports[hit];
        value = cplxld(source_data[index].real(), source_data[index].imag()) / rational_y_scale;
        if (rational_anchor >= 0)
        {
            denominator += rational_weights[hit];
            derivative = (numerator - value * denominator) / (rational_weights[hit] * (z - anchor));
        }
        else
            derivative = (numerator - value * denominator) / rational_weights[hit];
    }
    else
    {
        // Zero-weight supports must pass the actual ratio check as well.
        value = numerator / denominator;
        derivative =
            (d_numerator * denominator - numerator * d_denominator) / (denominator * denominator);
    }
    value *= rational_y_scale;
    derivative *= rational_y_scale / rational_x_scale;
    return {{static_cast<double>(value.real()), static_cast<double>(value.imag())},
            {static_cast<double>(derivative.real()), static_cast<double>(derivative.imag())}};
}

cplxdb
AnalyContPade::get(const cplxdb &x) const
{
    if (!rational_weights.empty()) return evaluate_rational(x).first;
    cplxdb tmp = {1.0, 0.0};

    for (int i_par = n_pars - 1; i_par > 0; i_par--)
    {
        tmp = 1.0 + par_y[i_par] * (x - source_xs[i_par-1]) / tmp;
    }
    const cplxdb value = par_y[0] / tmp;
    return finite_complex(value) ? value : evaluate_continuants(par_y, source_xs, x).first;
}

cplxdb
AnalyContPade::get_derivative(const cplxdb &x) const
{
    if (!rational_weights.empty()) return evaluate_rational(x).second;
    if (n_pars <= 1)
    {
        return {0.0, 0.0};
    }

    std::vector<cplxdb> g(n_pars);
    std::vector<cplxdb> dg(n_pars, {0.0, 0.0});

    g[n_pars-1] = par_y[n_pars-1];
    for (int i_par = n_pars - 2; i_par >= 0; i_par--)
    {
        const cplxdb denominator = 1.0 + (x - source_xs[i_par]) * g[i_par+1];
        g[i_par] = par_y[i_par] / denominator;
        dg[i_par] = -par_y[i_par]
                    * (g[i_par+1] + (x - source_xs[i_par]) * dg[i_par+1])
                    / (denominator * denominator);
    }
    return finite_complex(dg[0]) ? dg[0] : evaluate_continuants(par_y, source_xs, x).second;
}

const std::vector<double> get_specfunc(const AnalyCont &ac, const std::vector<cplxdb> &omegas,
                                       const double &ref, const double &e_ks, const double &v_xc,
                                       const double &v_exx,
                                       const double &sigc_omega_imag_shift,
                                       const double &gf_omega_imag_shift)
{
    const int n_freq = omegas.size();
    std::vector<double> sf(n_freq, 0.0);
    for (int i = 0; i < n_freq; i++)
    {
        const auto omega_gf = omegas[i] + cplxdb(0.0, gf_omega_imag_shift);
        const auto omega_sigc = omegas[i] + cplxdb(0.0, sigc_omega_imag_shift);
        cplxdb sf_c = omega_gf - e_ks + v_xc - v_exx - ac.get(omega_sigc - ref);
        sf[i] = - ((1.0 / PI) / sf_c).imag();
    }
    return sf;
}

}
