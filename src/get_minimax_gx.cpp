#include "get_minimax.h"

#include "gx_minimax_wrp.h"

bool check_minimax_available(int ngrids)
{
    bool avail = false;
    return avail;
}

void get_minimax_grid_frequency(int ngrids, double e_min, double e_max,
                                double *omega_points, double *omega_weights,
                                int &ierr)
{
    gx_minimax_grid_frequency_wrp(&ngrids, &e_min, &e_max, omega_points, omega_weights, &ierr);
}

void get_minimax_grid(int ngrids, double e_min, double e_max,
                      double *tau_points, double *tau_weights,
                      double *omega_points, double *omega_weights,
                      double *cosft_wt, double *cosft_tw, double *sinft_wt,
                      double max_errors[3], double &cosft_duality_error,
                      int &ierr)
{
    gx_minimax_grid_wrp(&ngrids, &e_min, &e_max, tau_points, tau_weights,
                        omega_points, omega_weights, cosft_wt, cosft_tw,
                        sinft_wt, max_errors, &cosft_duality_error, &ierr);
}
