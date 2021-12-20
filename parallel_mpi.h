#ifndef PARALLEL_MPI_H
#define PARALLEL_MPI_H

#include "matrix.h"
#include "complexmatrix.h"
#include "vector3_order.h"
#include <mpi.h>
#include <vector>
#include <map>
#include <memory>
class Parallel_MPI
{
    int myid, size;
    public:
    Parallel_MPI();
    ~Parallel_MPI();
    
    void mpi_init(int argc, char **argv);
    int get_myid(){return myid;}
    int get_size(){return size;}
    void allreduce_matrix(matrix &cmat_loc, matrix & cmat_glo);
    void allreduce_ComplexMatrix(ComplexMatrix &cmat_loc, ComplexMatrix & cmat_glo);
    void reduce_ComplexMatrix(ComplexMatrix &cmat_loc, ComplexMatrix & cmat_glo);
    vector<double> pack_mat(const map<size_t,map<size_t,map<Vector3_Order<int>,shared_ptr<matrix>>>> &Cs_m);
    map<size_t,map<size_t,map<Vector3_Order<int>,shared_ptr<matrix>>>> unpack_mat(vector<double> &pack);
};
#endif