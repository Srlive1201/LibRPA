#pragma once
#include "librpa_enums.h"

#include "../interface/mpi.h"

#include <string>
#include <vector>
#include <complex>

namespace librpa_int
{

// traits to decide MPI_Datatype for communication
template <typename T> struct mpi_datatype;
template <> struct mpi_datatype<int>
{ static constexpr MPI_Datatype value = MPI_INT; };
template <> struct mpi_datatype<float>
{ static constexpr MPI_Datatype value = MPI_FLOAT; };
template <> struct mpi_datatype<double>
{ static constexpr MPI_Datatype value = MPI_DOUBLE; };
template <> struct mpi_datatype<long>
{ static constexpr MPI_Datatype value = MPI_LONG; };
template <> struct mpi_datatype<std::complex<float>>
{ static constexpr MPI_Datatype value = MPI_C_FLOAT_COMPLEX; };
template <> struct mpi_datatype<std::complex<double>>
{ static constexpr MPI_Datatype value = MPI_C_DOUBLE_COMPLEX; };

class MpiCommHandler
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
    MpiCommHandler();
    MpiCommHandler(MPI_Comm comm_in, bool init_on_construct = false);
    ~MpiCommHandler() {};
    void init();
    void reset_comm();
    void reset_comm(MPI_Comm comm_in, bool init_on_reset = false);
    bool is_root() const { return this->myid == 0; }
    void barrier() const;
    std::string str() const;
    void check_initialized() const;

    // void allreduce_matrix(matrix &mat_send, matrix &mat_recv) const;
    // void allreduce_ComplexMatrix(ComplexMatrix &cmat_send, ComplexMatrix & cmat_recv) const;
    // void reduce_matrix(matrix &mat_send, matrix & cmat_recv, int root) const;
    // void reduce_ComplexMatrix(ComplexMatrix &cmat_send, ComplexMatrix & cmat_recv, int root) const;
};

// extern const std::string parallel_routing_notes[LIBRPA_ROUTING_COUNT];

// extern ParallelRouting parallel_routing;

// void set_parallel_routing(const std::string &option, const int &atpais_num, const int &Rt_num, ParallelRouting &routing);

//! Return the actual routing inside LibRPA when auto is selected
LibrpaParallelRouting decide_auto_routing(const int n_atoms, const int Rt_num);

//! Wrapper of MPI_Comm_rank
int get_mpi_rank(const MPI_Comm &comm);

//! Wrapper of MPI_Comm_size
int get_mpi_size(const MPI_Comm &comm);


//! task dispatchers, implemented single index and double indices versions. Ending indices are excluded
std::vector<int> dispatcher(int ist, int ied, unsigned myid, unsigned size, bool sequential);
std::vector<std::pair<int, int>> dispatcher(int i1st, int i1ed, int i2st, int i2ed,
                                  unsigned myid, unsigned size, bool sequential, bool favor_1st);

std::vector<std::pair<int,int>> pick_upper_triangular_tasks(std::vector<int> list_row, std::vector<int> list_col);
std::vector<std::pair<int,int>> dispatch_upper_triangular_tasks(const int &natoms, const int &myid, const int &nprows, const int &npcols, const int &myprow, const int &mypcol);

/*!
 * @brief find duplicate pairs across all processes
 * @param n: the maximum index in the ordered pair
 * @param ordered_pairs: the pairs on each process
 * @return the pairs to be removed in each process to ensure that one pair appears only once across the processes in comm
 */
std::vector<std::pair<int, int>> find_duplicate_ordered_pair(int n, const std::vector<std::pair<int, int>> &ordered_pairs, const MPI_Comm &comm);

template <typename T>
std::vector<T> dispatch_vector(std::vector<T> world_vec, unsigned myid, unsigned size, bool sequential)
{
    std::vector<int> ids = dispatcher(0, world_vec.size(), myid, size, sequential);
    std::vector<T> local_vec;
    for ( auto id: ids )
        local_vec.push_back(world_vec[id]);
    return local_vec;
}

template <typename T1, typename T2>
std::vector<std::pair<T1, T2>> dispatch_vector_prod(const std::vector<T1> &vec1, const std::vector<T2> &vec2, unsigned myid, unsigned size, bool sequential, bool favor_1st)
{
    std::vector<std::pair<int, int>> ids = dispatcher(0, int(vec1.size()), 0, int(vec2.size()), myid, size, sequential, favor_1st);
    std::vector<std::pair<T1, T2>> local_vec;
    for ( auto id: ids )
        local_vec.push_back({vec1[id.first], vec2[id.second]});
    return local_vec;
}

} /* end of namespace librpa_int */
