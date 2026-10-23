/*!
 * @file      analycont.h
 * @brief     Utilities for analytic continuation
 * @author    Min-Ye Zhang
 * @date      2024-04-23
 */
#pragma once
#include <utility>
#include <vector>

#include "../utils/base_utility.h"

namespace librpa_int
{

class AnalyCont
{
protected:
    std::vector<cplxdb> source_xs;
    std::vector<cplxdb> source_data;

public:
    virtual ~AnalyCont() = default;
    virtual cplxdb get(const cplxdb &x) const = 0;

    /*!
     * @brief get the complex source points used by the analytic continuation
     */
    const std::vector<cplxdb> &get_source_xs() const { return source_xs; };

    /*!
     * @brief get the complex source values used by the analytic continuation
     */
    const std::vector<cplxdb> &get_source_data() const { return source_data; };
};

class AnalyContPade: public AnalyCont
{
private:
    int n_pars;
    std::vector<cplxdb> par_y;
    std::vector<int> rational_supports;
    std::vector<std::complex<long double>> rational_weights;
    long double rational_x_scale = 1.0L;
    long double rational_y_scale = 1.0L;
    int rational_anchor = -1;

    bool initialize_rational();
    std::pair<cplxdb, cplxdb> evaluate_rational(const cplxdb &x) const;

public:
    AnalyContPade(int n_pars_in,
                  const std::vector<cplxdb> &xs,
                  const std::vector<cplxdb> &data);

    /*!
     * @brief get the value of continued function at complex number
     *
     * @param [in]    x    complex argument of function
     *
     * @return    a complex double, the value of function at x
     */
    cplxdb get(const cplxdb &x) const;

    /*!
     * @brief get the analytic derivative of continued function at complex number
     *
     * @param [in]    x    complex argument of function
     *
     * @return    a complex double, the derivative of function at x
     */
    cplxdb get_derivative(const cplxdb &x) const;
};

/*!
 * @brief get the spectral function from Pade approximant, based on diagonal approximation
 */
const std::vector<double> get_specfunc(const AnalyCont &ac, const std::vector<cplxdb> &omegas,
                                       const double &ref, const double &e_ks, const double &v_xc,
                                       const double &v_exx,
                                       const double &sigc_omega_imag_shift,
                                       const double &gf_omega_imag_shift);
}
