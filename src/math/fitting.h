/*!
 * \file fitting.h
 * \brief functionality for fitting
 */
#pragma once
#include <functional>
#include <vector>

namespace LIBRPA
{

namespace utils
{

//! Non-linear fitting using the Levenberg-Marquardt algorithm
struct LevMarqFitting
{
    //! Maximal number of iterations
    int    d_maxiter;
    double d_init_lambda;
    double d_up_factor;
    double d_down_factor;
    //! Target difference betwee errors of adjacent parameters estimate
    double d_target_derr;

    int    d_final_it;
    double d_final_err;
    double d_final_derr;

    // default constructor
    LevMarqFitting();
    // destructor
    ~LevMarqFitting() {};
    // copy & move constructor
    LevMarqFitting(const LevMarqFitting& s) = delete;
    LevMarqFitting(LevMarqFitting&& s) = delete;
    // copy & move assignment operator
    LevMarqFitting& operator=(const LevMarqFitting& s) = delete;
    LevMarqFitting& operator=(LevMarqFitting&& s) = delete;

    //! perform the fitting
    /*!
     * \param[inout] pars: initial parameters array, optimized after fitting
     * \param[in]    xs: abscissa
     * \param[in]    ys: values
     * \param[in]    func: function to fit the data against
     * \param[in]    grad: gradient with respect to parameters
     */
    void fit(std::vector<double> &pars, const std::vector<double> &xs,
             const std::vector<double> &ys,
             const std::function<double(double, const std::vector<double> &)> &func,
             const std::function<void(std::vector<double> &, double, const std::vector<double> &)> &grad);

    //! perform the fitting and evaluate the functin on a set of abscissa points
    std::vector<double> fit_eval(std::vector<double> &pars, const std::vector<double> &xs,
                  const std::vector<double> &ys,
                  const std::function<double(double, const std::vector<double> &)> &func,
                  const std::function<void(std::vector<double> &, double, const std::vector<double> &)> &grad,
                  const std::vector<double> &xs_eval);
};

} /* end of namespace utils */

} /* end of namespace LIBRPA */

