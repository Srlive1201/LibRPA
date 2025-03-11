#ifndef PARALLEL_MPI_H
#define PARALLEL_MPI_H

#include "matrix.h"
#include "complexmatrix.h"
#include "vector3_order.h"
#include <mpi.h>
#include <string>
#include <cassert>
#include <set>
#include <vector>
#include <utility>
#include <map>
#include <memory>

using std::vector;
using std::map;
using std::pair;

namespace LIBRPA {

enum ParallelRouting {
    ATOM_PAIR,
    R_TAU,
    LIBRI,
    COUNT
};

extern const string parallel_routing_notes[ParallelRouting::COUNT];

extern ParallelRouting parallel_routing;

void set_parallel_routing(const string &option, const int &atpais_num, const int &Rt_num, ParallelRouting &routing);

namespace MPI_Wrapper
{
    // 新增广播声明
    void broadcast_matrix(matrix &mat, int root, MPI_Comm mpi_comm);
    void broadcast_ComplexMatrix(ComplexMatrix &cmat, int root, MPI_Comm mpi_comm);
    void allreduce_matrix(const matrix& mat_send, matrix& mat_recv, MPI_Comm mpi_comm);
    void allreduce_ComplexMatrix(const ComplexMatrix& cmat_send, ComplexMatrix & cmat_recv, MPI_Comm mpi_comm);
    void reduce_matrix(const matrix& mat_send, matrix& mat_recv, int root, MPI_Comm mpi_comm);
    void reduce_ComplexMatrix(const ComplexMatrix& cmat_send, ComplexMatrix& cmat_recv, int root, MPI_Comm mpi_comm);
};

class MPI_COMM_handler
{
private:
    bool initialized_;
    bool comm_set_;
public:
    MPI_Comm comm;
    int myid;
    int nprocs;
    std::string procname;
public:
    MPI_COMM_handler();
    MPI_COMM_handler(MPI_Comm comm_in);
    ~MPI_COMM_handler() {};
    void init();
    void reset_comm();
    void reset_comm(MPI_Comm comm_in);
    bool is_root() const { return this->myid == 0; }
    void barrier() const;
    string str() const;
    void check_initialized() const;
    void allreduce_matrix(matrix &mat_send, matrix &mat_recv) const;
    void allreduce_ComplexMatrix(ComplexMatrix &cmat_send, ComplexMatrix & cmat_recv) const;
    void reduce_matrix(matrix &mat_send, matrix & cmat_recv, int root) const;
    void reduce_ComplexMatrix(ComplexMatrix &cmat_send, ComplexMatrix & cmat_recv, int root) const;
    template <typename T>
    void broadcast(T& data, int root = 0) const {
        MPI_Bcast(&data, sizeof(T), MPI_BYTE, root, comm);
    }
    
    // 添加这两个广播函数的声明
    void broadcast_matrix(matrix &mat, const int root = 0) const;
    void broadcast_ComplexMatrix(ComplexMatrix &cmat, const int root = 0) const;
};

} // namespace LIBRPA

class Parallel_MPI
{
public:
    static vector<double> pack_mat(const map<size_t,map<size_t,map<Vector3_Order<int>,shared_ptr<matrix>>>> &Cs_m);
    // static map<size_t,map<size_t,map<Vector3_Order<int>,shared_ptr<matrix>>>> unpack_mat(vector<double> &pack);

    Parallel_MPI();
    ~Parallel_MPI();

    // void set_blacs_parameters();
    // void set_blacs_mat(
    //     int *desc, int &loc_row, int &loc_col, 
    //     const int tot_row, const int tot_col,
    //     const int row_blk=1, const int col_blk=1 );
    static int globalIndex(int localIndex, int nblk, int nprocs, int myproc);
    static int localIndex(int globalIndex, int nblk, int nprocs, int myproc);
};

//! task dispatchers, implemented single index and double indices versions. Ending indices are excluded
vector<int> dispatcher(int ist, int ied, unsigned myid, unsigned size, bool sequential);
vector<pair<int, int>> dispatcher(int i1st, int i1ed, int i2st, int i2ed,
                                  unsigned myid, unsigned size, bool sequential, bool favor_1st);

vector<pair<int,int>> pick_upper_trangular_tasks(vector<int> list_row, vector<int> list_col);
vector<pair<int,int>> dispatch_upper_trangular_tasks(const int &natoms, const int &myid, const int &nprows, const int &npcols, const int &myprow, const int &mypcol);

/*!
 * @brief find duplicate pairs across all processes
 * @param n: the maximum index in the ordered pair
 * @param ordered_pairs: the pairs on each process
 * @return the pairs to be removed in each process to ensure that one pair appears only once across the processes in comm
 */
vector<pair<int, int>> find_duplicate_ordered_pair(int n, const vector<pair<int, int>>& ordered_pairs, const MPI_Comm &comm);

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

#endif
