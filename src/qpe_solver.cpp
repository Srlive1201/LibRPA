#include "qpe_solver.h"

#include <cmath>
#include <iostream>
#include <iomanip>

namespace LIBRPA
{

int qpe_solver_pade_self_consistent(
        const AnalyContPade &pade,
        const double &e_mf,
        const double &e_fermi,
        const double &vxc,
        const double &sigma_x,
        double &e_qp,
        cplxdb &sigc)
{
    int info = 0;
    const double escale = 0.1;
    const double thres = 1e-6;
    const int n_iter_max = 200;
    int n_iter = 0;

    // initial guess of e_qp as input mean-field energy
    e_qp = e_mf;

    double diff = 1e-3;

    // std::cout << "QPE: " << e_mf << " " << e_fermi << " " << vxc << " " << sigma_x << "\n";
    while (n_iter++ < n_iter_max)
    {
        e_qp = e_qp + escale * diff;
        sigc = pade.get(static_cast<cplxdb>(e_qp - e_fermi));
        diff = e_mf - vxc + sigma_x + sigc.real() - e_qp;
        if (std::abs(diff) < thres)
        {
            break;
        }
    }
    // std::cout << "Finished QPE solve: " << n_iter << " " << std::scientific << diff << " " << std::abs(diff) << "\n";

    if (n_iter > n_iter_max)
    {
        info = 1;
    }

    return info;
}

} /* end of namespace LIBRPA */
