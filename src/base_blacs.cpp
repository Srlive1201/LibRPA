#include "base_blacs.h"
#include "base_utility.h"

#include "interface/blacs_scalapack.h"
#include "utils_io.h"


namespace LIBRPA
{

void CTXT_barrier(int ictxt, CTXT_SCOPE scope)
{
    char scope_ch;
    switch (scope)
    {
        case (CTXT_SCOPE::R): scope_ch = 'R';
        case (CTXT_SCOPE::C): scope_ch = 'C';
        case (CTXT_SCOPE::A): scope_ch = 'A';
    }
    Cblacs_barrier(ictxt, &scope_ch);
}


void BLACS_CTXT_handler::init()
{
    this->mpi_comm_h.init();
    this->ictxt = Csys2blacs_handle(this->mpi_comm_h.comm);
    Cblacs_pinfo(&this->myid, &this->nprocs);
    this->initialized_ = true;
}

void BLACS_CTXT_handler::reset_comm(MPI_Comm comm_in)
{
    this->mpi_comm_h.reset_comm(comm_in);
    this->comm_set_ = true;
    this->pgrid_set_ = false;
    this->initialized_ = false;
}

void BLACS_CTXT_handler::set_grid(const int &nprows_in, const int &npcols_in,
                                  CTXT_LAYOUT layout_in)
{
    // if the grid has been set, exit it first
    if (pgrid_set_) exit();
    if (nprocs != nprows_in * npcols_in)
        throw std::invalid_argument("nprocs != nprows * npcols");
    layout = layout_in;
    if (layout == CTXT_LAYOUT::C)
        layout_ch = 'C';
    else
        layout_ch = 'R';
    Cblacs_gridinit(&ictxt, &layout_ch, nprows_in, npcols_in);
    Cblacs_gridinfo(ictxt, &nprows, &npcols, &myprow, &mypcol);
    pgrid_set_ = true;
}

void BLACS_CTXT_handler::set_square_grid(bool more_rows, CTXT_LAYOUT layout_in)
{
    int nroc;
    layout = layout_in;
    for (nroc = int(sqrt(double(nprocs))); nroc >= 2; --nroc)
    {
        if ((nprocs) % nroc == 0) break;
    }
    more_rows ? nprows = nprocs / (npcols = nroc)
              : npcols = nprocs / (nprows = nroc);
    set_grid(nprows, npcols, layout_in);
}

void BLACS_CTXT_handler::set_horizontal_grid()
{
    set_grid(1, nprocs, CTXT_LAYOUT::R);
}

void BLACS_CTXT_handler::set_vertical_grid()
{
    set_grid(nprocs, 1, CTXT_LAYOUT::C);
}

void BLACS_CTXT_handler::exit()
{
    if (pgrid_set_)
    {
        Cblacs_gridexit(ictxt);
        // recollect the system context
        ictxt = Csys2blacs_handle(mpi_comm_h.comm);
        pgrid_set_ = false;
    }
}

std::string BLACS_CTXT_handler::info() const
{
    std::string info;
    info = std::string("BLACS_CTXT_handler: ")
         + "ICTXT " + std::to_string(ictxt) + " "
         + "PSIZE " + std::to_string(nprocs) + " "
         + "PID " + std::to_string(myid) + " "
         + "PGRID (" + std::to_string(nprows) + "," + std::to_string(npcols) + ") "
         + "PCOOD (" + std::to_string(myprow) + "," + std::to_string(mypcol) +")";
    return info;
}

int BLACS_CTXT_handler::get_pnum(int prow, int pcol) const
{
    return Cblacs_pnum(ictxt, prow, pcol);
}

void BLACS_CTXT_handler::get_pcoord(int pid, int &prow, int &pcol) const
{
    Cblacs_pcoord(ictxt, pid, &prow, &pcol);
}

void BLACS_CTXT_handler::barrier(CTXT_SCOPE scope) const
{
    CTXT_barrier(ictxt, scope);
}

std::vector<int> get_proc_indices_blacs(const int &m, const int &n, const int &mb, const int &nb,
                                        const int &irsrc, const int &icsrc,
                                        const int &ictxt, bool row_fast)
{
    std::vector<int> procids;
    procids.reserve(as_size(m * n));
    // myprow/mypcol are not used here, declared for indxg2p interface
    int nprows, npcols, myprow, mypcol;
    Cblacs_gridinfo(ictxt, &nprows, &npcols, &myprow, &mypcol);

    if (row_fast)
    {
        for (auto j = 0; j < n; j++)
        {
            int pcol = ScalapackConnector::indxg2p(j, nb, mypcol, icsrc, npcols);
            for (auto i = 0; i < m; i++)
            {
                int prow = ScalapackConnector::indxg2p(i, mb, myprow, irsrc, nprows);
                int rank = Cblacs_pnum(ictxt, prow, pcol);
                procids.push_back(rank);
            }
        }
    }
    else
    {
        for (auto i = 0; i < m; i++)
        {
            int prow = ScalapackConnector::indxg2p(i, mb, myprow, irsrc, nprows);
            for (auto j = 0; j < n; j++)
            {
                int pcol = ScalapackConnector::indxg2p(j, nb, mypcol, icsrc, npcols);
                int rank = Cblacs_pnum(ictxt, prow, pcol);
                procids.push_back(rank);
            }
        }
    }

    return procids;
}

void Array_Desc::set_blacs_params_(int comm, int ictxt, int nprocs, int myid, int nprows,
                                  int myprow, int npcols, int mypcol)
{
    assert(myid < nprocs && myprow < nprows && mypcol < npcols);
    comm_ = comm;
    ictxt_ = ictxt;
    nprocs_ = nprocs;
    myid_ = myid;
    nprows_ = nprows;
    myprow_ = myprow;
    npcols_ = npcols;
    mypcol_ = mypcol;
}

int Array_Desc::set_desc_(const int &m, const int &n, const int &mb, const int &nb,
                          const int &irsrc, const int &icsrc)
{
    int info = 0;
    m_local_ = ScalapackConnector::numroc(m, mb, myprow_, irsrc, nprows_);
    // leading dimension
    lld_ = std::max(m_local_, 1);
    n_local_ = ScalapackConnector::numroc(n, nb, mypcol_, icsrc, npcols_);
    if (m_local_ < 1 || n_local_ < 1)
    {
        empty_local_mat_ = true;
    }

    ScalapackConnector::descinit(this->desc, m, n, mb, nb, irsrc, icsrc, ictxt_, lld_, info);
    if (info)
    {
        LIBRPA::utils::lib_printf(
            "ERROR DESCINIT! PROC %d (%d,%d) PARAMS: DESC %d %d %d %d %d %d %d %d\n",
            myid_, myprow_, mypcol_, m, n, mb, nb, irsrc, icsrc, ictxt_, m_local_);
    }
    // else
    //     LIBRPA::utils::lib_printf("SUCCE DESCINIT! PROC %d (%d,%d) PARAMS: DESC %d %d %d %d %d %d %d %d\n", myid_, myprow_, mypcol_, m, n, mb, nb, irsrc, icsrc, ictxt_, m_local_);
    m_ = desc[2];
    n_ = desc[3];
    mb_ = desc[4];
    nb_ = desc[5];
    irsrc_ = desc[6];
    icsrc_ = desc[7];
    lld_ = desc[8];
    initialized_ = true;
    return info;
}

Array_Desc::Array_Desc(const BLACS_CTXT_handler &blacs_h)
    : ictxt_(0), nprocs_(0), myid_(0),
      nprows_(0), myprow_(0), npcols_(0), mypcol_(0),
      m_(0), n_(0), mb_(0), nb_(0), irsrc_(0), icsrc_(0),
      lld_(0), m_local_(0), n_local_(0),
      empty_local_mat_(false), initialized_(false)
{
    if (!blacs_h.initialized())
        throw std::logic_error("BLACS context is not initialized before creating ArrayDesc");
    set_blacs_params_(blacs_h.get_comm(), blacs_h.ictxt, blacs_h.nprocs, blacs_h.myid,
                      blacs_h.nprows, blacs_h.myprow, blacs_h.npcols,
                      blacs_h.mypcol);
}

Array_Desc::Array_Desc(const int &ictxt)
    : ictxt_(0), nprocs_(0), myid_(0),
      nprows_(0), myprow_(0), npcols_(0), mypcol_(0),
      m_(0), n_(0), mb_(0), nb_(0), irsrc_(0), icsrc_(0),
      lld_(0), m_local_(0), n_local_(0),
      empty_local_mat_(false), initialized_(false)
{
    // TODO: how to check if ictxt is a valid context?
    int nprocs, myid, nprows, npcols, myprow, mypcol;
    MPI_Comm comm = Cblacs2sys_handle(ictxt);
    Cblacs_gridinfo(ictxt, &nprows, &npcols, &myprow, &mypcol);
    myid = Cblacs_pnum(ictxt, myprow, mypcol);
    nprocs = nprows * npcols;
    set_blacs_params_(comm, ictxt, nprocs, myid,
                      nprows, myprow, npcols,
                      mypcol);
}

int Array_Desc::init(const int &m, const int &n, const int &mb, const int &nb,
                    const int &irsrc, const int &icsrc)
{
    return set_desc_(m, n, mb, nb, irsrc, icsrc);
}

int Array_Desc::init_1b1p(const int &m, const int &n,
                          const int &irsrc, const int &icsrc)
{
    int mb = 1, nb = 1;
    mb = std::ceil(double(m)/nprows_);
    nb = std::ceil(double(n)/npcols_);
    return set_desc_(m, n, mb, nb, irsrc, icsrc);
}

int Array_Desc::init_square_blk(const int &m, const int &n,
                                    const int &irsrc, const int &icsrc)
{
    int mb = 1, nb = 1, minblk = 1;
    mb = std::ceil(double(m)/nprows_);
    nb = std::ceil(double(n)/npcols_);
    minblk = std::min(mb, nb);
    return set_desc_(m, n, minblk, minblk, irsrc, icsrc);
}

std::string Array_Desc::info() const
{
    std::string info;
    info = std::string("ArrayDesc: ")
         + "ICTXT " + std::to_string(ictxt_) + " "
         + "ID " + std::to_string(myid_) + " "
         + "PCOOR (" + std::to_string(myprow_) + "," + std::to_string(mypcol_) + ") "
         + "GSIZE (" + std::to_string(m_) + "," + std::to_string(n_) + ") "
         + "LSIZE (" + std::to_string(m_local_) + "," + std::to_string(n_local_) + ") "
         + "DUMMY? " + std::string(empty_local_mat_? "T" : "F");
    return info;
}

std::string Array_Desc::info_desc() const
{
    char s[100];
    sprintf(s, "DESC %d %d %d %d %d %d %d %d %d",
            desc[0], desc[1], desc[2],
            desc[3], desc[4], desc[5],
            desc[6], desc[7], desc[8]);
    return std::string(s);
}

void Array_Desc::barrier(CTXT_SCOPE scope)
{
    CTXT_barrier(ictxt_, scope);
}

int get_globalIndex(int localIndex, int nblk, int nprocs, int myproc)
{
    int iblock, gIndex;
    iblock = localIndex / nblk;
    gIndex = (iblock * nprocs + myproc) * nblk + localIndex % nblk;
    return gIndex;
}

int get_localIndex(int globalIndex, int nblk, int nprocs, int myproc)
{
    int inproc = int((globalIndex % (nblk * nprocs)) / nblk);
    if (myproc == inproc)
    {
        return int(globalIndex / (nblk * nprocs)) * nblk + globalIndex % nblk;
    }
    else
    {
        return -1;
    }
}

std::vector<int> get_proc_indices_blacs(const Array_Desc &ad, bool row_fast)
{
    assert(ad.initialized());
    return get_proc_indices_blacs(ad.m(), ad.n(), ad.mb(), ad.nb(), ad.irsrc(), ad.icsrc(),
                                  ad.ictxt(), row_fast);
}

std::vector<std::pair<size_t, size_t>>
get_2d_mat_indices_blacs(const int &m, const int &n,
                         const int &mb, const int &nb,
                         const int &irsrc, const int &icsrc,
                         const int &nprows, const int &npcols,
                         const int &myprow, const int &mypcol, bool row_fast)
{
    std::vector<std::pair<size_t, size_t>> indices;

    std::vector<size_t> row_ids;
    std::vector<size_t> col_ids;

    for (auto i = 0; i < m; i++)
    {
        if (myprow == ScalapackConnector::indxg2p(i, mb, myprow, irsrc, nprows))
        {
            row_ids.push_back(as_size(i));
        }
    }
    for (auto i = 0; i < n; i++)
    {
        if (mypcol == ScalapackConnector::indxg2p(i, nb, mypcol, icsrc, npcols))
        {
            col_ids.push_back(as_size(i));
        }
    }

    if (row_fast)
    {
        for (const auto &c: col_ids)
        {
            for (const auto &r: row_ids)
            {
                indices.push_back({r, c});
            }
        }
    }
    else
    {
        for (const auto &r: row_ids)
        {
            for (const auto &c: col_ids)
            {
                indices.push_back({r, c});
            }
        }
    }

    return indices;
}

std::vector<std::pair<size_t, size_t>> get_2d_mat_indices_blacs(const Array_Desc &ad, const int &myid, bool row_fast)
{
    assert(ad.initialized());

    int myprow, mypcol;
    Cblacs_pcoord(ad.ictxt(), myid, &myprow, &mypcol);

    return get_2d_mat_indices_blacs(ad.m(), ad.n(), ad.mb(), ad.nb(), ad.irsrc(), ad.icsrc(),
                                    ad.nprows(), ad.npcols(), myprow, mypcol, row_fast);
}

std::vector<size_t>
get_1d_mat_indices_blacs(const int &m, const int &n,
                         const int &mb, const int &nb,
                         const int &irsrc, const int &icsrc,
                         const int &nprows, const int &npcols,
                         const int &myprow, const int &mypcol,
                         bool row_fast, bool row_major)
{
    const auto indices_2d = get_2d_mat_indices_blacs(m, n, mb, nb, irsrc, icsrc, nprows, npcols, myprow, mypcol, row_fast);
    return flatten_2d_indices(indices_2d, as_size(m), as_size(n),
                              row_major);
}

std::vector<size_t>
get_1d_mat_indices_blacs(const Array_Desc &ad, const int &myid, bool row_fast, bool row_major)
{
    const auto indices_2d = get_2d_mat_indices_blacs(ad, myid, row_fast);
    return flatten_2d_indices(indices_2d, as_size(ad.m()), as_size(ad.n()),
                              row_major);
}

} /* end of namespace LIBRPA */
