/*!
 * @file interpolate.h
 * @brief Component to perform interpolation
 * @details This component includes basic functions to interpolate values from
 *          available grid data. The implementation is probably naive. For better
 *          performance and validity, it should be replaced by GSL or other GPL-licensed
 *          library when license issue is resolved.
 */
#pragma once
#include <vector>

namespace LIBRPA
{
namespace utils
{

//! Class to perform cubic spline interolation
class CubicSpline
{
private:
    const std::vector<double> d_xs;
    const std::vector<double> d_ys;
    std::vector<double> d_a;
    std::vector<double> d_b;
    std::vector<double> d_c;
    std::vector<double> d_d;
public:
    // default constructor
    CubicSpline(const std::vector<double>& xs, const std::vector<double>& ys);
    // destructor
    ~CubicSpline() {};
    // copy & move constructor
    CubicSpline(const CubicSpline& s) = delete;
    CubicSpline(CubicSpline&& s) = delete;
    // copy & move assignment operator
    CubicSpline& operator=(const CubicSpline& s) = delete;
    CubicSpline& operator=(CubicSpline&& s) = delete;

    //! interpolate on individual x
    double operator()(double x) const;
    //! interpolate on a list of x
    std::vector<double> operator()(const std::vector<double>& xs) const;
}; /* end of class CubicSpline */

/*!
 * @brief Helper function to performance cubic spline interpolation
 * @param[in] x values of avaiable data points
 * @param[in] y values of avaiable data points
 * @param[in] x values of points to interpolate
 * @return y values of points to interpolate
 */
std::vector<double> interp_cubic_spline(const std::vector<double> &xs,
                                        const std::vector<double> &ys,
                                        const std::vector<double> &xs_new);

} /* end of namepsace utils */

} /* end of namespace LIBRPA */
