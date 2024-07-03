#pragma once

//! \brief get the minimax frequency grids
void get_minimax_grid_frequency(int ngrids, double e_min, double e_max,
                                double *omega_points, double *omega_weights,
                                int &ierr);

//! \brief get the minimax time and frequency grids and related transformation matrices
void get_minimax_grid(int ngrids, double e_min, double e_max,
                      double *tau_points, double *tau_weights,
                      double *omega_points, double *omega_weights,
                      double *cosft_wt, double *cosft_tw, double *sinft_wt,
                      double max_errors[3], double &cosft_duality_error, int &ierr);
