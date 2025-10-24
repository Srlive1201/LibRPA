#pragma once

#include <string>
#include <valarray>
#include <vector>

#include "base_mpi.h"
#include "scalapack_connector.h"

namespace LIBRPA
{

enum class CTXT_LAYOUT {R, C};
enum class CTXT_SCOPE {R, C, A};

void CTXT_barrier(int ictxt, CTXT_SCOPE scope = CTXT_SCOPE::A);

class BLACS_CTXT_handler
{
private:
    MPI_COMM_handler mpi_comm_h;
    char layout_ch;
    bool initialized_;
    bool pgrid_set_;
    bool comm_set_;
public:
    int ictxt;
    CTXT_LAYOUT layout;
    int myid;
    int nprocs;
    int nprows;
    int npcols;
    int mypcol;
    int myprow;
    BLACS_CTXT_handler() { comm_set_ = pgrid_set_ = initialized_ = false; }
    BLACS_CTXT_handler(MPI_Comm comm_in): mpi_comm_h(comm_in) { comm_set_ = true; pgrid_set_ = initialized_ = false; }
    ~BLACS_CTXT_handler() {};

    void init();
    void reset_comm(MPI_Comm comm_in);
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
    const MPI_Comm &get_comm() const { return this->mpi_comm_h.comm; };
};

/*!
 * @brief Get indices of process to which the distributed elements of global matrix belong
 *
 * @param  [in]  m, n            Number of rows and columns of the global matrix
 * @param  [in]  mb, nb          Block size along row and column direction.
 * @param  [in]  irsrc, icsrc    Source process along row and column.
 * @param  [in]  ictxt           BLACS context, must be initialized
 * @param  [in]  row_fast        Flag to set the row basis index goes faster.
 *
 * @retval       indices         List of process indices (size m * n)
 */
std::vector<int> get_proc_indices_blacs(const int &m, const int &n, const int &mb, const int &nb,
                                        const int &irsrc, const int &icsrc,
                                        const int &ictxt, bool row_fast);

class Array_Desc
{
private:
    // BLACS parameters obtained upon construction
    MPI_Comm comm_;
    int ictxt_;
    int nprocs_;
    int myid_;
    int nprows_;
    int myprow_;
    int npcols_;
    int mypcol_;
    void set_blacs_params_(MPI_Comm comm, int ictxt, int nprocs, int myid, int nprows, int myprow, int npcols, int mypcol);
    int set_desc_indices_(const int &m, const int &n, const int &mb, const int &nb,
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

    // Precomputed indices
    std::valarray<int> g2l_r;
    std::valarray<int> g2l_c;
    std::valarray<int> l2g_r;
    std::valarray<int> l2g_c;

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
    inline int indx_g2l_r(int gindx) const noexcept { return g2l_r[gindx]; };
    inline int indx_g2l_c(int gindx) const noexcept { return g2l_c[gindx]; };
    inline int indx_l2g_r(int lindx) const noexcept { return l2g_r[lindx]; };
    inline int indx_l2g_c(int lindx) const noexcept { return l2g_c[lindx]; };
    const int& myid() const { return myid_; }
    const int& ictxt() const { return ictxt_; }
    const MPI_Comm& comm() const { return comm_; }
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
    const int& nprocs() const { return nprocs_; }
    const int& nprows() const { return nprows_; }
    const int& npcols() const { return npcols_; }
    std::string info() const;
    std::string info_desc() const;
    bool is_src() const;
    bool initialized() const { return this->initialized_; };
    void barrier(CTXT_SCOPE scope = CTXT_SCOPE::A);
};


int get_globalIndex(int localIndex, int nblk, int nprocs, int myproc);
int get_localIndex(int globalIndex, int nblk, int nprocs, int myproc);

/*!
 * @brief Get indices of process to which the distributed elements of global matrix belong
 *
 * @param  [in]  ad              Array_Desc object, must be initialized beforehand.
 * @param  [in]  row_fast        Flag to set the row basis index goes faster.
 *
 * @warning      Memory intensive for large system
 *
 * @retval       indices         List of process indices (size m * n)
 */
std::vector<int> get_proc_indices_blacs(const Array_Desc &ad, bool row_fast);

/*!
 * @brief Get 2D indices of elements of submatrix in BLACS 2D block-cyclic format
 *
 * @param  [in]  m, n            Size of rows and columns.
 * @param  [in]  mb, nb          Block size along row and column direction.
 * @param  [in]  irsrc, icsrc    Source process along row and column.
 * @param  [in]  nprows, npcols  Number of processes along row and column,
 *                               i.e. the shape of proccess grid.
 * @param  [in]  myprow, mypcol  Coordinate of current process in a process grid.
 * @param  [in]  row_fast        Flag to set the row basis index faster.
 *
 * @retval       indices
 */
std::vector<std::pair<size_t, size_t>>
get_2d_mat_indices_blacs(const int &m, const int &n,
                         const int &mb, const int &nb,
                         const int &irsrc, const int &icsrc,
                         const int &nprows, const int &npcols,
                         const int &myprow, const int &mypcol, bool row_fast);

/*!
 * @brief Get 2D indices of elements of submatrix in BLACS 2D block-cyclic format
 *        from the Array_Desc object.
 *
 * @param  [in]  ad         Array_Desc object, must be initialized before parsing
 * @param  [in]  row_fast   Flag to set the row basis index faster.
 *
 * @retval indices
 */
std::vector<std::pair<size_t, size_t>>
get_2d_mat_indices_blacs(const Array_Desc &ad, const int &myid, bool row_fast);

/*!
 * @brief Get 1D indices of elements of submatrix in BLACS 2D block-cyclic format
 *
 * @param  [in]  m, n            Size of rows and columns.
 * @param  [in]  mb, nb          Block size along row and column direction.
 * @param  [in]  irsrc, icsrc    Source process along row and column.
 * @param  [in]  nprows, npcols  Number of processes along row and column,
 *                               i.e. the shape of proccess grid.
 * @param  [in]  myprow, mypcol  Coordinate of current process in a process grid.
 * @param  [in]  row_fast        Flag to set the row basis index faster.
 * @param  [in]  row_major       Flag to compute the 1D indices in row-major (C-style)
 *
 * @retval       indices
 */
std::vector<size_t>
get_1d_mat_indices_blacs(const int &m, const int &n,
                         const int &mb, const int &nb,
                         const int &irsrc, const int &icsrc,
                         const int &nprows, const int &npcols,
                         const int &myprow, const int &mypcol,
                         bool row_fast, bool row_major);

/*!
 * @brief Get 1D indices of elements of submatrix in BLACS 2D block-cyclic format
 *        from the Array_Desc object.
 *
 * @param  [in]  ad         Array_Desc object, must be initialized before parsing
 * @param  [in]  row_fast   Flag to set the row basis index faster.
 * @param  [in]  row_major  Flag to compute the 1D indices in row-major (C-style)
 *
 * @retval indices
 */
std::vector<size_t>
get_1d_mat_indices_blacs(const Array_Desc &ad, const int &myid, bool row_fast, bool row_major);

} /* end of namespace LIBRPA */
