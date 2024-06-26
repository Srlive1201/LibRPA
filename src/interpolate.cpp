#include "interpolate.h"
#include <cmath>
#include <cstdio>
#include <stdexcept>

namespace LIBRPA
{
namespace utils
{
CubicSpline::CubicSpline(const std::vector<double>& xs, const std::vector<double>& ys): d_xs(xs), d_ys(ys)
{
    // Get the number of data points
    int n = d_xs.size();
    d_a.resize(n);
    d_b.resize(n);
    d_c.resize(n);
    d_d.resize(n);
    if (n < 2)
        throw std::logic_error("CubicSpline require at least 2 points to interpolate");

    // Create and initialize the vectors h, alpha, l, and mu
    std::vector<double> h(n - 1);
    std::vector<double> alpha(n - 1);
    std::vector<double> l(n);
    std::vector<double> mu(n - 1);

    for (int i = 0; i < n - 1; ++i)
    {
        h[i] = d_xs[i + 1] - d_xs[i];
        alpha[i] = 3 * (d_ys[i + 1] - d_ys[i]) / h[i];
    }

    // Initialize l and mu
    l[0] = 1;
    mu[0] = 0;

    // Calculate l and mu using the Thomas algorithm
    for (int i = 1; i < n - 1; ++i)
    {
        l[i] = 2 * (d_xs[i + 1] - d_xs[i - 1]) - h[i - 1] * mu[i - 1];
        mu[i] = h[i] / l[i];
    }

    l[n - 1] = 1;

    // Create and initialize the vectors z and c
    std::vector<double> z(n);
    d_c[n - 1] = 0;

    // Calculate z and c using the Thomas algorithm
    for (int i = 1; i < n; ++i)
    {
        z[i] = (alpha[i - 1] - h[i - 1] * z[i - 1]) / l[i];
        d_c[n - i - 1] = z[n - i] - mu[n - i - 1] * d_c[n - i];
    }

    // Calculate the coefficients of the natural cubic spline interpolating polynomial
    for (int i = 0; i < n - 1; ++i)
    {
        d_a[i] = d_ys[i];
        d_b[i] = (d_ys[i + 1] - d_ys[i]) / h[i] - h[i] * (d_c[i + 1] + 2 * d_c[i]) / 3;
        d_d[i] = (d_c[i + 1] - d_c[i]) / (3 * h[i]);
    }

    d_a[n - 1] = d_ys[n - 1];
    d_b[n - 1] = 0;
    d_d[n - 1] = 0;
}

// This function evaluates the natural cubic spline interpolating polynomial at a given point x
double CubicSpline::operator()(double x) const
{
    int n = d_xs.size();

    // Check if the point is outside the range of the available data points
    if (x < d_xs[0] || x > d_xs[n - 1])
    {
        // Extrapolate the natural cubic spline polynomial
        double x0, x1, y0, y1, m0, m1, h;
        int i;

        if (x < d_xs[0])
        {
            // Extrapolate to the left of the available data points
            x0 = d_xs[0];
            x1 = d_xs[1];
            y0 = d_a[0];
            y1 = d_a[1];
            m0 = d_b[0];
            m1 = d_b[1];
            h = x1 - x0;
            i = 0;
        }
        else
        {
            // Extrapolate to the right of the available data points
            x0 = d_xs[n - 2];
            x1 = d_xs[n - 1];
            y0 = d_a[n - 2];
            y1 = d_a[n - 1];
            m0 = d_b[n - 2];
            m1 = d_b[n - 1];
            h = x1 - x0;
            i = n - 2;
        }

        double x_minus_xi = x - x0;
        double y_minus_yi = y1 - y0;
        double m_minus_mi = m1 - m0;

        double result = y0 + m0 * x_minus_xi + ((3 * y_minus_yi / std::pow(h, 2)) - (2 * m0 / h) - (m_minus_mi / h)) * std::pow(x_minus_xi, 2) + ((-2 * y_minus_yi / std::pow(h, 3)) + (m0 / std::pow(h, 2)) + (m_minus_mi / std::pow(h, 2))) * std::pow(x_minus_xi, 3);

        return result;
    }
    else
    {
        // Evaluate the natural cubic spline polynomial for the given point
        for (int i = 0; i < n - 1; ++i)
        {
            if (x >= d_xs[i] && x <= d_xs[i + 1])
            {
                double xi = d_xs[i];
                double hi = d_xs[i + 1] - d_xs[i];
                double yi = d_a[i];
                double bi = d_b[i];
                double ci = d_c[i];
                double di = d_d[i];

                double x_minus_xi = x - xi;
                double result = yi + bi * x_minus_xi + ci * std::pow(x_minus_xi, 2) + di * std::pow(x_minus_xi, 3);

                return result;
            }
        }
    }
    return 0.0;
}

std::vector<double> CubicSpline::operator()(const std::vector<double>& xs) const
{
    std::vector<double> ys;
    for (auto x: xs)
    {
        ys.push_back((*this)(x));
    }
    return ys;
}

std::vector<double> interp_cubic_spline(const std::vector<double> &xs, const std::vector<double> &ys, const std::vector<double> &xs_new)
{
    const CubicSpline cs(xs, ys);
    return cs(xs_new);
}

} /* end of namespace utils */
} /* end of namespace LIBRPA */

