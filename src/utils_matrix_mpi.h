#pragma once

#include "matrix.h"
#include "complexmatrix.h"
#include "base_mpi.h"

namespace LIBRPA
{

void reduce_matrix(const matrix& mat_send, matrix& mat_recv, int root, MPI_Comm mpi_comm);
void reduce_ComplexMatrix(const ComplexMatrix& cmat_send, ComplexMatrix& cmat_recv, int root, MPI_Comm mpi_comm);

void allreduce_matrix(const matrix& mat_send, matrix& mat_recv, MPI_Comm mpi_comm);
void allreduce_ComplexMatrix(const ComplexMatrix& cmat_send, ComplexMatrix & cmat_recv, MPI_Comm mpi_comm);


} /* end of namespace LIBRPA */
