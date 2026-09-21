#pragma once

#include "../mpi/base_blacs.h"
#include "../mpi/base_mpi.h"
#include "complexmatrix.h"
#include "matrix.h"
#include "matrix_m.h"

namespace librpa_int
{

void reduce_matrix(const matrix& mat_send, matrix& mat_recv, int root, MPI_Comm mpi_comm);
void reduce_ComplexMatrix(const ComplexMatrix& cmat_send, ComplexMatrix& cmat_recv, int root, MPI_Comm mpi_comm);

void allreduce_matrix(const matrix& mat_send, matrix& mat_recv, MPI_Comm mpi_comm);
void allreduce_ComplexMatrix(const ComplexMatrix& cmat_send, ComplexMatrix & cmat_recv, MPI_Comm mpi_comm);

void broadcast_matrix(matrix &mat, const int root, MPI_Comm mpi_comm);
void broadcast_ComplexMatrix(ComplexMatrix &cmat, const int root, MPI_Comm mpi_comm);

// Gather on the descriptor source process; return an empty matrix elsewhere.
// All processes in the BLACS context must participate with column-major data.
Matz collect_blacs_matrix_root(const Matz& local, const ArrayDesc& distributed_descriptor);

} /* end of namespace librpa_int */
