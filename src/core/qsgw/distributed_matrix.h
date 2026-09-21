#pragma once

#include "matrix_map.h"

#include "../../mpi/base_mpi.h"

namespace librpa_int
{
namespace qsgw
{

void broadcast_spin_k_matrix_map(SpinKMatrixMap& values,
                                 int root,
                                 const MpiCommHandler& communicator);

} // namespace qsgw
} // namespace librpa_int
