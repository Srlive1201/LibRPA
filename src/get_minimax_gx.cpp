#include "get_minimax.h"

#include "gx_minimax_wrp.h"

// NOTE(MYZ): rescale frequency weights for number of grids below
//            some point due to inconsistency in greenx
// NOTE(MYZ): Fixed upstream, so there is no need to do this any more
// static const int ngrids_below_freqweight_inconsistent = 26;
// static const double scale_inconsistent = 0.25;

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
    // if (ngrids < ngrids_below_freqweight_inconsistent)
    //     for (int i = 0; i != ngrids; i++)
    //         omega_weights[i] *= scale_inconsistent;
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
    // if (ngrids < ngrids_below_freqweight_inconsistent)
    //     for (int i = 0; i != ngrids; i++)
    //         omega_weights[i] *= scale_inconsistent;
}
