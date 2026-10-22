#pragma once

#include <vector>

#include "../core/atom.h"
#include "../core/geometry.h"
#include "../core/pbc.h"
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

}
