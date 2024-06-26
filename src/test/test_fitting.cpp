#include "../fitting.h"
#include "testutils.h"

#include <cmath>
#include <cassert>
#include <functional>
#include <iostream>

void test_levmarq_sin()
{
    LIBRPA::utils::LevMarqFitting levmarq;

    // sin(ax)
    auto f_sin =
        [](double x, const std::vector<double> &pars) { return std::sin(x*pars[0]); };
    // d sin(ax) / d a = x cos(ax)
    auto g_sin =
        [](std::vector<double> &gs, double x, const std::vector<double> &pars) { gs[0] = x * std::cos(x*pars[0]); };
    const int npars = 1;
    std::vector<double> pars(npars);

    int ndata = 100;
    std::vector<double> xs(ndata);
    std::vector<double> ys(ndata);

    std::vector<double> alphas { 0.5, 1.5, 2.0, -2.0, };
    for (const auto &a: alphas)
    {
        for (int i = 0; i != ndata; i++)
        {
            xs[i] = (double(i) / (ndata - 1) - 0.5) * M_PI;
            ys[i] = sin(a * xs[i]);
        }
        // initialize
        pars[0] = 1.0;
        levmarq.fit(pars, xs, ys, f_sin, g_sin);
        printf("Fitted params: %f\n", pars[0]);
        assert(fequal(pars[0], a, 1e-8));
    }
}

int main(int argc, char *argv[])
{
    test_levmarq_sin();
}
