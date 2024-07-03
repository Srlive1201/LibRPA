#include "fitting.h"

#include <cmath>

#include "utils_io.h"

/*
 * The fitting code is adapted from C version by Ron Babich
 * See https://gist.github.com/rbabich/3539146
 */

static double error_func(const std::vector<double> &par,
                         const std::vector<double> &xs,
                         const std::vector<double> &ys,
                         const std::function<double(double, const std::vector<double>&)> &func)
{
    int i;
    double res, e = 0;
    const std::size_t n = ys.size();

    for (i = 0; i < n; i++)
    {
        res = func(xs[i], par) - ys[i];
        e += res * res;
    }
    return e;
}

/*
 * This function takes a symmetric, positive-definite matrix "a" and returns
 * its (lower-triangular) Cholesky factor in "l".  Elements above the
 * diagonal are neither used nor modified.  The same array may be passed
 * as both l and a, in which case the decomposition is performed in place.
 */
static int cholesky_decomp(const int &n,
                           std::vector<std::vector<double>> &l,
                           std::vector<std::vector<double>> &a)
{
    int i, j, k;
    double sum;
    const double TOL = 1e-30;

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < i; j++)
        {
            sum = 0;
            for (k = 0; k < j; k++) sum += l[i][k] * l[j][k];
            l[i][j] = (a[i][j] - sum) / l[j][j];
        }

        sum = 0;
        for (k = 0; k < i; k++) sum += l[i][k] * l[i][k];
        sum = a[i][i] - sum;
        if (sum < TOL) return 1; /* not positive-definite */
        l[i][i] = std::sqrt(sum);
    }
    return 0;
}

/*
 * solve the equation Ax=b for a symmetric positive-definite matrix A,
 * using the Cholesky decomposition A=LL^T.  The matrix L is passed in "l".
 * Elements above the diagonal are ignored.
 */
static void solve_axb_cholesky(int n, const std::vector<std::vector<double>> &l, std::vector<double> &x,
                        const std::vector<double> &b)
{
    int i, j;
    double sum;

    /* solve L*y = b for y (where x[] is used to store y) */
    for (i = 0; i < n; i++)
    {
        sum = 0;
        for (j = 0; j < i; j++) sum += l[i][j] * x[j];
        x[i] = (b[i] - sum) / l[i][i];
    }

    /* solve L^T*x = y for x (where x[] is used to store both y and x) */
    for (i = n - 1; i >= 0; i--)
    {
        sum = 0;
        for (j = i + 1; j < n; j++) sum += l[j][i] * x[j];
        x[i] = (x[i] - sum) / l[i][i];
    }
}

namespace LIBRPA
{

namespace utils
{

LevMarqFitting::LevMarqFitting()
{
    d_maxiter = 10000;
    d_init_lambda = 0.1;
    d_up_factor = 10;
    d_down_factor = 10;
    d_target_derr = 1e-8;
}

void LevMarqFitting::fit(
    std::vector<double> &pars, const std::vector<double> &xs, const std::vector<double> &ys,
    const std::function<double(double, const std::vector<double> &)> &func,
    const std::function<void(std::vector<double> &, double, const std::vector<double> &)> &grad)
{
    const int npars = pars.size();
    const std::size_t n = xs.size();
    if (xs.size() != ys.size())
    {
        lib_printf("Warning: inconsistent x and y size of data\n");
    }
    int nit = d_maxiter;
    double lambda = d_init_lambda;
    double up = d_up_factor;
    double down = 1 / d_down_factor;
    double target_derr = d_target_derr;
    double weight = 1;
    double derr;
    double newerr = 0; /* to avoid compiler warnings */

    double err = error_func(pars, xs, ys, func);

    std::vector<std::vector<double>> ch(npars, std::vector<double>(npars));
    std::vector<std::vector<double>> h(npars, std::vector<double>(npars));
    std::vector<double> g(npars);
    std::vector<double> newpar(npars);
    std::vector<double> d(npars);
    std::vector<double> delta(npars);

    /* main iteration */
    for (int it = 0; it < nit; it++)
    {
        /* calculate the approximation to the Hessian and the "derivative" d */
        for (int i = 0; i < npars; i++)
        {
            d[i] = 0;
            for (int j = 0; j <= i; j++) h[i][j] = 0;
        }
        for (int ix = 0; ix < n; ix++)
        {
            grad(g, xs[ix], pars);
            for (int i = 0; i < npars; i++)
            {
                d[i] += (ys[ix] - func(xs[ix], pars)) * g[i] * weight;
                for (int j = 0; j <= i; j++) h[i][j] += g[i] * g[j] * weight;
            }
        }

        /*
         * make a step "delta."  If the step is rejected, increase
         * lambda and try again
         */
        double mult = 1 + lambda;
        int ill = 1; /* ill-conditioned? */
        while (ill && (it < nit))
        {
            for (int i = 0; i < npars; i++) h[i][i] = h[i][i] * mult;

            // cholesky_decomp(int n, double *(*l), double *(*a));
            ill = cholesky_decomp(npars, ch, h);

            if (!ill)
            {
                solve_axb_cholesky(npars, ch, delta, d);
                for (int i = 0; i < npars; i++) newpar[i] = pars[i] + delta[i];
                newerr = error_func(newpar, xs, ys, func);
                derr = newerr - err;
                ill = (derr > 0);
            }
            // if (verbose)
            //     lib_printf(
            //         "it = %4d,   lambda = %10g,   err = %10g,   "
            //         "derr = %10g\n",
            //         it, lambda, err, derr);
            if (ill)
            {
                mult = (1 + lambda * up) / (1 + lambda);
                lambda *= up;
                it++;
            }
        }
        for (int i = 0; i < npars; i++) pars[i] = newpar[i];
        err = newerr;
        lambda *= down;
#ifdef PRINT_DEBUG
        lib_printf("Iteration %d: ", it);
        for (i = 0; i < npar; i++)
        {
            lib_printf("%f:", par[i]);
        }
        lib_printf("Error: %f\n", derr);
#endif
        if ((!ill) && (-derr < target_derr))
        {
#ifdef PRINT_DEBUG
            lib_printf("Converge after: %d cycles, final error: %f\n", it, derr);
#endif
            break;
        }
    }
}

std::vector<double> LevMarqFitting::fit_eval(std::vector<double> &pars, const std::vector<double> &xs, const std::vector<double> &ys, const std::function<double (double, const std::vector<double> &)> &func, const std::function<void (std::vector<double> &, double, const std::vector<double> &)> &grad, const std::vector<double> &xs_eval)
{
    fit(pars, xs, ys, func, grad);
    std::vector<double> ys_eval(xs_eval.size());
    for (int i = 0; i != ys_eval.size(); i++)
    {
        ys_eval[i] = func(xs_eval[i], pars);
    }
    return ys_eval;
}

} /* end of namespace utils */

} /* end of namespace LIBRPA */
