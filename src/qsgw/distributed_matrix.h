#pragma once

#include "matrix_map.h"

#include "../math/matrix_m.h"
#include "../mpi/base_blacs.h"
#include "../mpi/base_mpi.h"

namespace librpa_int
{
namespace qsgw
{

Matz collect_blacs_matrix_root(const Matz& local,
                               const ArrayDesc& distributed_descriptor);

void broadcast_spin_k_matrix_map(SpinKMatrixMap& values,
                                 int root,
                                 const MpiCommHandler& communicator);

} // namespace qsgw
} // namespace librpa_int
