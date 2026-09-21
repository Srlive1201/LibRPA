#pragma once

#include <string>

#include "../math/matrix_m.h"
#include "../mpi/base_blacs.h"

namespace librpa_int
{

//! Collectively export a square complex matrix, or the same half-open index
//! range in both dimensions. A negative index_end selects all remaining indices.
//! Local blocks must use column-major ScaLAPACK storage. The descriptor's source
//! rank writes the file; I/O failure is reported on every rank in desc.comm().
//! Format (native byte order): int32 dimension, int32 sizeof(double), then
//! row-major pairs of real/imaginary doubles, with no threshold filtering.
void write_matrix_binary_parallel(const Matz &mat_loc, const ArrayDesc &desc, const std::string &fn,
                                  int index_start = 0, int index_end = -1);

}  // namespace librpa_int
