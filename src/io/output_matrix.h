#pragma once

#include <string>

#include "../core/atomic_basis.h"
#include "../core/pbc.h"
#include "../math/matrix_m.h"
#include "../mpi/base_blacs.h"

namespace librpa_int
{

//! Export one set of complete real-space AO blocks, without threshold filtering.
//! Native binary layout: size_t block count, then for each block size_t[5]
//! (R index in pbc.Rlist, I, J, rows, columns) and row-major complex doubles.
//! This is the existing Sigma_c(R,iw) checkpoint format. Returns false on I/O
//! failure; callers handle MPI error propagation. Blocks must be row-major.
bool write_rspace_matrices_binary(const ap_p_map<std::map<Vector3_Order<int>, Matz>> &blocks,
                                  const AtomicBasis &basis, const PeriodicBoundaryData &pbc,
                                  const std::string &fn);

//! Collectively export a square complex matrix, or the same half-open index
//! range in both dimensions. A negative index_end selects all remaining indices.
//! Local blocks must use column-major ScaLAPACK storage. The descriptor's source
//! rank writes the file; I/O failure is reported on every rank in desc.comm().
//! Format (native byte order): int32 dimension, int32 sizeof(double), then
//! row-major pairs of real/imaginary doubles, with no threshold filtering.
void write_matrix_binary_parallel(const Matz &mat_loc, const ArrayDesc &desc, const std::string &fn,
                                  int index_start = 0, int index_end = -1);

}  // namespace librpa_int
