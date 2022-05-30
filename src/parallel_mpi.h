#ifndef PARALLEL_MPI_H
#define PARALLEL_MPI_H

#include "matrix.h"
#include "complexmatrix.h"
#include "vector3_order.h"
#include <mpi.h>
#include <string>
#include <cassert>
#include <vector>
#include <utility>
#include <map>
#include <memory>

using std::vector;
using std::map;
using std::pair;

class Parallel_MPI
{
    int myid, size;
    public:
    enum parallel_type { ATOM_PAIR, R_TAU };
    parallel_type chi_parallel_type;
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

//! task dispatchers, implemented single index and double indices versions. Ending indices are excluded
vector<int> dispatcher(int ist, int ied, unsigned myid, unsigned size, bool sequential);
vector<pair<int, int>> dispatcher(int i1st, int i1ed, int i2st, int i2ed,
                                  unsigned myid, unsigned size, bool sequential, bool favor_1st);

template <typename T>
vector<T> dispatch_vector(vector<T> world_vec, unsigned myid, unsigned size, bool sequential)
{
    vector<int> ids = dispatcher(0, world_vec.size(), myid, size, sequential);
    vector<T> local_vec;
    for ( auto id: ids )
        local_vec.push_back(world_vec[id]);
    return local_vec;
}

extern Parallel_MPI para_mpi;
#endif
