#pragma once

#include <librpa_options.h>

#include <vector>

#include "../core/atom.h"
#include "../core/geometry.h"
#include "../core/meanfield.h"
#include "../core/pbc.h"
#include "../core/timefreq.h"
#include "../mpi/base_mpi.h"

namespace librpa_int::api
{

//! Collect the union of k-point indices requested by all ranks in an extraction API call.
//! The result gives a compact, deterministic k-list for collective result handoff.
std::vector<int> collect_requested_iks(const MpiCommHandler &comm_h, int n_kpts_this,
                                       const int *iks_this, int n_kpoints);

//! Build the atom-pair BvK remap requested by a band interpolation API call.
AtomPairBvKRemap<atom_t> build_band_bvk_remap(const Atoms &atoms,
                                              const PeriodicBoundaryData &pbc,
                                              int remap_convention);

//! Shared G0W0 continuation and QPE evaluation, without building/projection or MPI.
//! Sigma diagonal layout: [spin][collected k][frequency][selected state].
//! Potentials and outputs: [spin][requested k][selected state], in Hartree.
//! All inputs refer to the same mean-field eigenbasis; solver_info retains failures
//! even when override_qpe_solver_nan requests the final unconverged value.
void solve_g0w0_qpe_diagonal(const MeanField &mf, const TFGrids &tfg, const LibrpaOptions &opts,
                             const std::vector<cplxdb> &sigc_diag,
                             const std::vector<int> &iks_collect, int n_kpts_this,
                             const int *iks_this, int i_state_low, int n_states_calc,
                             const double *vxc, const double *vexx, double *sigc_re,
                             double *sigc_im, double *eqp, const char *context,
                             int *solver_info = nullptr);
}
