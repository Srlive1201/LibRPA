#ifndef PARALLEL_MPI_H
#define PARALLEL_MPI_H

#include "atomic_basis.h"
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

enum parallel_type {
    ATOM_PAIR,
    R_TAU,
    LIBRI_USED,
    COUNT
};

extern const string parallel_types_note[parallel_type::COUNT];

extern parallel_type chi_parallel_type;
extern parallel_type exx_parallel_type;

extern ofstream fout_para;

void set_parallel_type(const string &option, parallel_type &ptype);
void set_chi_parallel_type(const string &option, const int &atpais_num, const int Rt_num, const bool use_libri);
void set_exx_parallel_type(const string &option, const int &atpais_num, const int Rt_num, const bool use_libri);
void check_parallel_type();

namespace MPI_Wrapper
{
    extern std::string procname;
    extern bool initialized;
    extern int myid_world;
    extern int nprocs_world;
    bool is_root_world();
    void init(int argc, char **argv);
    void init(MPI_Comm comm_in);
    void finalize();
    void barrier_world();
    void allreduce_matrix(const matrix& mat_send, matrix& mat_recv, MPI_Comm mpi_comm);
    void allreduce_ComplexMatrix(const ComplexMatrix& cmat_send, ComplexMatrix & cmat_recv, MPI_Comm mpi_comm);
    void reduce_matrix(const matrix& mat_send, matrix& mat_recv, int root, MPI_Comm mpi_comm);
    void reduce_ComplexMatrix(const ComplexMatrix& cmat_send, ComplexMatrix& cmat_recv, int root, MPI_Comm mpi_comm);
};

class MPI_COMM_handler
{
public:
    MPI_Comm comm;
    int myid;
    int nprocs;
    bool initialized;
    void check_initialized() const;
public:
    MPI_COMM_handler(MPI_Comm comm_in);
    ~MPI_COMM_handler() {};
    void init();
    bool is_root() const { return this->myid == 0; }
    void barrier() const;
    string str() const;
    void allreduce_matrix(matrix &mat_send, matrix &mat_recv) const;
    void allreduce_ComplexMatrix(ComplexMatrix &cmat_send, ComplexMatrix & cmat_recv) const;
    void reduce_matrix(matrix &mat_send, matrix & cmat_recv, int root) const;
    void reduce_ComplexMatrix(ComplexMatrix &cmat_send, ComplexMatrix & cmat_recv, int root) const;
};

extern MPI_COMM_handler mpi_comm_world_h;

enum class CTXT_LAYOUT {R, C};
enum class CTXT_SCOPE {R, C, A};

void CTXT_barrier(int ictxt, CTXT_SCOPE scope = CTXT_SCOPE::A);

class BLACS_CTXT_handler
{
private:
    MPI_COMM_handler mpi_comm_h;
    char layout_ch;
    bool initialized_;
    bool pgrid_set;
public:
    int ictxt;
    CTXT_LAYOUT layout;
    int myid;
    int nprocs;
    int nprows;
    int npcols;
    int mypcol;
    int myprow;
    BLACS_CTXT_handler(MPI_Comm comm_in): mpi_comm_h(comm_in) { pgrid_set = initialized_ = false; }
    ~BLACS_CTXT_handler() {};

    void init();
    void set_grid(const int &nprows_in, const int &npcols_in, CTXT_LAYOUT layout_in = CTXT_LAYOUT::R);
    void set_square_grid(bool more_rows = true, CTXT_LAYOUT layout_in = CTXT_LAYOUT::R);
    void set_horizontal_grid();
    void set_vertical_grid();
    std::string info() const;
    int get_pnum(int prow, int pcol) const;
    void get_pcoord(int pid, int &prow, int &pcol) const;
    void barrier(CTXT_SCOPE scope = CTXT_SCOPE::A) const;
    //! call gridexit to reset process grid
    void exit();
    bool initialized() const { return initialized_; }
};

extern BLACS_CTXT_handler blacs_ctxt_world_h;

class Array_Desc
{
private:
    // BLACS parameters obtained upon construction
    int ictxt_;
    int nprocs_;
    int myid_;
    int nprows_;
    int myprow_;
    int npcols_;
    int mypcol_;
    void set_blacs_params_(int ictxt, int nprocs, int myid, int nprows, int myprow, int npcols, int mypcol);
    int set_desc_(const int &m, const int &n, const int &mb, const int &nb,
                  const int &irsrc, const int &icsrc);

    // Array dimensions
    int m_;
    int n_;
    int mb_;
    int nb_;
    int irsrc_;
    int icsrc_;
    int lld_;
    int m_local_;
    int n_local_;

    //! flag to indicate that the current process should contain no data of local matrix, but for scalapack routines, it will generate a dummy matrix of size 1, nrows = ncols = 1
    bool empty_local_mat_ = false;

    //! flag for initialization
    bool initialized_ = false;

public:
    int desc[9];
    Array_Desc(const BLACS_CTXT_handler &blacs_ctxt_h);
    Array_Desc(const int &ictxt);
    //! initialize the array descriptor
    int init(const int &m, const int &n,
             const int &mb, const int &nb,
             const int &irsrc, const int &icsrc);
    //! initialize the array descriptor such that each process has exactly one block
    int init_1b1p(const int &m, const int &n,
                  const int &irsrc, const int &icsrc);
    int init_square_blk(const int &m, const int &n,
                        const int &irsrc, const int &icsrc);
    int indx_g2l_r(int gindx) const;
    int indx_g2l_c(int gindx) const;
    int indx_l2g_r(int lindx) const;
    int indx_l2g_c(int lindx) const;
    const int& ictxt() const { return ictxt_; }
    const int& m() const { return m_; }
    const int& n() const { return n_; }
    const int& mb() const { return mb_; }
    const int& nb() const { return nb_; }
    const int& lld() const { return lld_; }
    const int& irsrc() const { return irsrc_; }
    const int& icsrc() const { return icsrc_; }
    const int& m_loc() const { return m_local_; }
    const int& n_loc() const { return n_local_; }
    const int& myprow() const { return myprow_; }
    const int& mypcol() const { return mypcol_; }
    const int& nprows() const { return nprows_; }
    const int& npcols() const { return npcols_; }
    std::string info() const;
    std::string info_desc() const;
    bool is_src() const;
    void barrier(CTXT_SCOPE scope = CTXT_SCOPE::A);
};

//! prepare array descriptors for distributing(collecting) submatrices
//! from(to) a full matrix on source process with p?gemr2d
std::pair<Array_Desc, Array_Desc> prepare_array_desc_mr2d_src_and_all(
    const BLACS_CTXT_handler &ctxt_h, const int &m, const int &n, const int &mb,
    const int &nb, const int &irsrc, const int &icsrc);

//! obtain the necessary atom pair of atomic basis to build the block-cyclic submatrix
std::set<std::pair<int, int>> get_necessary_IJ_from_block_2D(const AtomicBasis &atbasis_row, const AtomicBasis &atbasis_col, const Array_Desc& arrdesc);

std::set<std::pair<int, int>> get_necessary_IJ_from_block_2D_sy(const char &uplo, const AtomicBasis &atbasis, const Array_Desc& arrdesc);

} // namespace LIBRPA

class Parallel_MPI
{
public:
    enum parallel_type { ATOM_PAIR, R_TAU, LIBRI_USED };
    parallel_type chi_parallel_type;
    
    static vector<double> pack_mat(const map<size_t,map<size_t,map<Vector3_Order<int>,shared_ptr<matrix>>>> &Cs_m);
    static map<size_t,map<size_t,map<Vector3_Order<int>,shared_ptr<matrix>>>> unpack_mat(vector<double> &pack);

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
