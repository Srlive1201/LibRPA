#include "qpe_solver.h"
namespace LIBRPA
{

int qpe_linear_solver_pade(
        const AnalyContPade &pade,
        const double &e_mf,
        const double &e_fermi,
        const double &vxc,
        const double &sigma_x,
        double &e_qp)
{
    int info = 0;
    const double escale = 0.1;
    const double thres = 1e-6;
    const int n_iter_max = 200;
    int n_iter = 0;

    // initial guess of e_qp as input mean-field energy
    e_qp = e_mf;

    double diff = 1e-3;

    while (n_iter++ < n_iter_max)
    {
        e_qp += escale * diff;
        auto sigc = pade.get(static_cast<cplxdb>(e_qp - e_fermi)).real();
        diff = e_mf - vxc + sigma_x + sigc - e_qp;
        if (abs(diff) < thres)
        {
            break;
        }
    }

    if (n_iter > n_iter_max)
    {
        info = 1;
    }

    return info;
}

} /* end of namespace LIBRPA */
