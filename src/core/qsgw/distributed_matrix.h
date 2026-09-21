#pragma once

#include "../../mpi/base_mpi.h"
#include "matrix_map.h"

namespace librpa_int
{
namespace qsgw
{

void broadcast_spin_k_matrix_map(SpinKMatrixMap& values,
                                 int root,
                                 const MpiCommHandler& communicator);

} // namespace qsgw
} // namespace librpa_int
