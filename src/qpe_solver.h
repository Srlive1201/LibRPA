#pragma once

#include "analycont.h"

namespace LIBRPA
{

/*!
 * \brief Solve quasi-particle equation self-consistently by computing self-energy on real freqeuncy through Pade Analytic continuation
 *
 * \param [in]     pade       AnalyContPade object, constructed from correlation self-energy at imaginary frequency
 * \param [in]     e_mf       Energy of the state of meanf-field calculation
 * \param [in]     e_fermi    Fermi energy
 * \param [in]     vxc        Exchange-correlation potential in the mean-field calculation
 * \param [in]     sigma_x    Exchange self-energy
 * \param [out]    e_qp       Quasi-particle energy as the solution of QPE
 * \param [out]    sigc       Correlation self-energy of the quasi-particle
 * \retval         info       0 if QPE is solved successfully
 */
int qpe_solver_pade_self_consistent(
        const AnalyContPade &pade,
        const double &e_mf,
        const double &e_fermi,
        const double &vxc,
        const double &sigma_x,
        double &e_qp,
        cplxdb &sigc,
        double thres = 1.0e-5);

}
