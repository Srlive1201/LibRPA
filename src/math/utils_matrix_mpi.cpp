#include "utils_matrix_mpi.h"

#include "../utils/base_utility.h"
#include "complexmatrix.h"
#include <cassert>
#include <complex>

namespace librpa_int
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

void broadcast_matrix(matrix &mat, int root, MPI_Comm mpi_comm)
{
    int dims[2];
    dims[0] = mat.nr;
    dims[1] = mat.nc;
    MPI_Bcast(dims, 2, MPI_INT, root, mpi_comm);

    int myrank;
    MPI_Comm_rank(mpi_comm, &myrank);
    if (myrank != root)
    {
        mat.create(dims[0], dims[1]);
    }

    size_t count = as_size(dims[0]) * as_size(dims[1]);
    MPI_Bcast(mat.c, count, mpi_datatype<double>::value, root, mpi_comm);
}

void broadcast_ComplexMatrix(ComplexMatrix &cmat, int root, MPI_Comm mpi_comm)
{
    int dims[2];
    dims[0] = cmat.nr;
    dims[1] = cmat.nc;
    MPI_Bcast(dims, 2, MPI_INT, root, mpi_comm);

    int myrank;
    MPI_Comm_rank(mpi_comm, &myrank);
    if (myrank != root)
    {
        cmat.create(dims[0], dims[1]);
    }

    size_t count = as_size(dims[0]) * as_size(dims[1]);
    MPI_Bcast(cmat.c, count, mpi_datatype<std::complex<double>>::value, root, mpi_comm);
}

} /* end of namespace librpa_int */

