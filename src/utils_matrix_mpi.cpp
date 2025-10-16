#include "utils_matrix_mpi.h"

#include <cassert>

namespace LIBRPA
{

void allreduce_matrix(const matrix &mat_send, matrix &mat_recv, MPI_Comm mpi_comm)
{
    assert(mat_send.nr==mat_recv.nr);
    assert(mat_send.nc==mat_recv.nc);
    int mat_size=mat_recv.nr*mat_recv.nc;
    MPI_Allreduce(mat_send.c, mat_recv.c, mat_size, MPI_DOUBLE, MPI_SUM, mpi_comm);
}

void reduce_matrix(const matrix &mat_send, matrix &mat_recv, int root, MPI_Comm mpi_comm)
{
    assert(mat_send.nr==mat_recv.nr);
    assert(mat_send.nc==mat_recv.nc);
    int mat_size = mat_recv.nr * mat_recv.nc;
    MPI_Reduce(mat_send.c, mat_recv.c, mat_size, MPI_DOUBLE, MPI_SUM, root, mpi_comm);
}

void allreduce_ComplexMatrix(const ComplexMatrix &cmat_send, ComplexMatrix &cmat_recv, MPI_Comm mpi_comm)
{
    assert(cmat_send.nr==cmat_recv.nr);
    assert(cmat_send.nc==cmat_recv.nc);
    MPI_Allreduce(cmat_send.c, cmat_recv.c, cmat_recv.size, MPI_DOUBLE_COMPLEX, MPI_SUM, mpi_comm);
}

void reduce_ComplexMatrix(const ComplexMatrix &cmat_send, ComplexMatrix &cmat_recv, int root, MPI_Comm mpi_comm)
{
    assert(cmat_send.nr==cmat_recv.nr);
    assert(cmat_send.nc==cmat_recv.nc);
    MPI_Reduce(cmat_send.c, cmat_recv.c, cmat_recv.size, MPI_DOUBLE_COMPLEX, MPI_SUM, root, mpi_comm);
}

} /* end of namespace LIBRPA */

