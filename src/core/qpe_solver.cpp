#include "qpe_solver.h"

#include <cmath>
// #include <iostream>
// #include <iomanip>

namespace librpa_int
{

int qpe_solver_pade_self_consistent(
        const AnalyContPade &pade,
        const double &e_mf,
        const double &e_fermi,
        const double &vxc,
        const double &sigma_x,
        double &e_qp,
        cplxdb &sigc,
        double thres)
{
    int info = 0;
    const double escale = 0.1;
    const int n_iter_max = 200;
    int n_iter = 0;


    // initial guess of e_qp as input mean-field energy
    double e_qp_last = e_mf, e_qp_this;
    cplxdb sigc_this, sigc_last;

    double diff = 1e-3;

    // std::cout << "QPE: " << e_mf << " " << e_fermi << " " << vxc << " " << sigma_x << "\n";
    while (n_iter++ < n_iter_max)
    {
        e_qp_this = e_qp_last + escale * diff;
        sigc_this = pade.get(static_cast<cplxdb>(e_qp_this - e_fermi));
        diff = e_mf - vxc + sigma_x + sigc_this.real() - e_qp_this;
        e_qp_last = e_mf - vxc + sigma_x + sigc_this.real();
        sigc_last = sigc_this;
        if (std::abs(diff) < thres)
        {
            break;
        }
    }
    // std::cout << "Finished QPE solve: " << n_iter << " " << std::scientific << diff << " " << std::abs(diff) << "\n";
    sigc = sigc_last;
    e_qp = e_qp_last;
    // std::cout << "Compare QPE solver: " << std::setw(16) << std::setprecision(8) << e_mf - vxc + sigma_x + sigc.real() << " " << e_qp << std::endl;

    if (n_iter > n_iter_max)
    {
        info = 1;
    }

    return info;
}

} /* end of namespace librpa_int */
