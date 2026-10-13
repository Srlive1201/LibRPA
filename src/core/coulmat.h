/*!
 * @file coulmat.h
 * @brief utilities for handling with Coulomb matrix
 */
#pragma once
#include "../mpi/base_mpi.h"
#include "atomic_basis.h"
#include "symmetry_context.h"
#include "pbc.h"
#include "ri.h"

namespace librpa_int {

/*!
 * @brief perform Fourier transform of q/k-space Coulomb matrix to R-space
 * @param comm_h communicator used for distributed symmetry collection
 * @param coulmat  q/k-space Coulomb matrix by atom-pair index
 * @param Rlist  vectors of R where the R-space Coulomb matrix should be computed
 * @param return_ordered_atom_pair  ensure that the returned type
 *        will always have keys of ordered atom pairs,
 *        i.e. there always exists map[I][J] and map[J][I] for I not equal to J
 * @param result_is_replicated optional ownership receipt: true only when every
 *        rank receives the same irreducible-sector map. Such maps must be
 *        uniquely partitioned before an additive downstream communication.
 * @return  atpair_R_mat_t
 */
atpair_R_mat_t FT_Vq(const MpiCommHandler &comm_h,
                     const AtomicBasis &basis_abf,
                     const SymmetryContext &symmetry_context,
                     const atpair_k_cplx_mat_t &coulmat_k,
                     const PeriodicBoundaryData &pbc,
                     bool return_ordered_atom_pair,
                     bool use_symmetry_context = true,
                     bool* result_is_replicated = nullptr);

}
