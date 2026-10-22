#pragma once

#include <vector>

#include "../meanfield.h"
#include "../timefreq.h"
#include "librpa_options.h"

namespace librpa_int::qsgw
{
//! Solve the entering QSGW state's diagonal QP diagnostic, retaining root status.
void solve_qsgw_qpe_diagonal(
    const MeanField& mf, const TFGrids& tfg, const LibrpaOptions& opts,
    const std::vector<cplxdb>& sigc_diag, const std::vector<int>& iks_collect, int n_kpts_this,
    const int* iks_this, int i_state_low, int n_states_calc, const double* vxc, const double* vexx,
    double* sigc_re, double* sigc_im, double* eqp, const char* context, int* solver_info = nullptr);
}  // namespace librpa_int::qsgw
