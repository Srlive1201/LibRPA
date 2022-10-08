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
    int my_blacs_ctxt;
    int myprow, mypcol, nprow, npcol;
    char BLACS_LAYOUT;
    enum parallel_type { ATOM_PAIR, R_TAU };
    parallel_type chi_parallel_type;
    Parallel_MPI();
    ~Parallel_MPI();
    
    void mpi_init(int argc, char **argv);
    void mpi_barrier();
    int get_myid(){return myid;}
    int get_size(){return size;}
    void allreduce_matrix(matrix &cmat_loc, matrix & cmat_glo);
    void allreduce_ComplexMatrix(ComplexMatrix &cmat_loc, ComplexMatrix & cmat_glo);
    void reduce_ComplexMatrix(ComplexMatrix &cmat_loc, ComplexMatrix & cmat_glo);
    vector<double> pack_mat(const map<size_t,map<size_t,map<Vector3_Order<int>,shared_ptr<matrix>>>> &Cs_m);
    map<size_t,map<size_t,map<Vector3_Order<int>,shared_ptr<matrix>>>> unpack_mat(vector<double> &pack);

    void set_blacs_parameters();
    void set_blacs_mat(
        int *desc, int &loc_row, int &loc_col, 
        const int tot_row, const int tot_col,
        const int row_blk=1, const int col_blk=1 );
    int globalIndex(int localIndex, int nblk, int nprocs, int myproc);
    int localIndex(int globalIndex, int nblk, int nprocs, int myproc);
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

template <typename T1, typename T2>
vector<pair<T1, T2>> dispatch_vector_prod(vector<T1> vec1, vector<T2> vec2, unsigned myid, unsigned size, bool sequential, bool favor_1st)
{
    vector<pair<int, int>> ids = dispatcher(0, int(vec1.size()), 0, int(vec2.size()), myid, size, sequential, favor_1st);
    vector<pair<T1, T2>> local_vec;
    for ( auto id: ids )
        local_vec.push_back({vec1[id.first], vec2[id.second]});
    return local_vec;
}

extern Parallel_MPI para_mpi;
#endif
