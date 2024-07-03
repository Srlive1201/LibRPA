#include "parallel_mpi.h"

#include <algorithm>
#include <stdexcept>

#include "interface/blacs_scalapack.h"
#include "scalapack_connector.h"
#include "utils_io.h"

namespace LIBRPA {

const string parallel_routing_notes[ParallelRouting::COUNT] = {
    "atom-pair",
    "R-tau",
    "LibRI"
};

ParallelRouting parallel_routing = ParallelRouting::ATOM_PAIR;

void set_parallel_routing(const string &option, const int &atpais_num, const int &Rt_num, ParallelRouting &routing)
{
    if (option == "auto")
        routing = atpais_num < Rt_num ? ParallelRouting::R_TAU : ParallelRouting::ATOM_PAIR;
    else if (option == "atompair")
        routing = ParallelRouting::ATOM_PAIR;
    else if (option == "rtau")
        routing = ParallelRouting::R_TAU;
    else if (option == "libri")
        routing = ParallelRouting::LIBRI;
    else
        throw std::invalid_argument("Unsupported parallel type option: " + option);
}


namespace MPI_Wrapper {

void allreduce_matrix(matrix &mat_send, matrix &mat_recv, MPI_Comm mpi_comm)
{
    assert(mat_send.nr==mat_recv.nr);
    assert(mat_send.nc==mat_recv.nc);
    int mat_size=mat_recv.nr*mat_recv.nc;
    MPI_Allreduce(mat_send.c, mat_recv.c, mat_size, MPI_DOUBLE, MPI_SUM, mpi_comm);
}

void allreduce_ComplexMatrix(ComplexMatrix &cmat_send, ComplexMatrix &cmat_recv, MPI_Comm mpi_comm)
{
    assert(cmat_send.nr==cmat_recv.nr);
    assert(cmat_send.nc==cmat_recv.nc);
    MPI_Allreduce(cmat_send.c, cmat_recv.c, cmat_recv.size, MPI_DOUBLE_COMPLEX, MPI_SUM, mpi_comm);
}

void reduce_matrix(matrix &mat_send, matrix &mat_recv, int root, MPI_Comm mpi_comm)
{
    assert(mat_send.nr==mat_recv.nr);
    assert(mat_send.nc==mat_recv.nc);
    int mat_size = mat_recv.nr * mat_recv.nc;
    MPI_Reduce(mat_send.c, mat_recv.c, mat_size, MPI_DOUBLE, MPI_SUM, root, mpi_comm);
}

void reduce_ComplexMatrix(ComplexMatrix &cmat_send, ComplexMatrix &cmat_recv, int root, MPI_Comm mpi_comm)
{
    assert(cmat_send.nr==cmat_recv.nr);
    assert(cmat_send.nc==cmat_recv.nc);
    MPI_Reduce(cmat_send.c, cmat_recv.c, cmat_recv.size, MPI_DOUBLE_COMPLEX, MPI_SUM, root, mpi_comm);
}

}  // namespace MPI_Wraper

MPI_COMM_handler::MPI_COMM_handler()
{
    this->comm_set_ = false;
    this->initialized_ = false; 
}

MPI_COMM_handler::MPI_COMM_handler(MPI_Comm comm_in)
        : comm(comm_in)
{
    this->comm_set_ = true;
    this->initialized_ = false; 
}

void MPI_COMM_handler::reset_comm(MPI_Comm comm_in)
{
    this->comm = comm_in; 
    this->comm_set_ = true;
    this->initialized_ = false; 
}

void MPI_COMM_handler::check_initialized() const
{
    if (!initialized_)
        throw std::logic_error("MPI_COMM_handler not initialized");
}

void MPI_COMM_handler::init()
{
    MPI_Comm_rank(this->comm, &(this->myid));
    MPI_Comm_size(this->comm, &(this->nprocs));
    this->initialized_ = true;
}

void MPI_COMM_handler::barrier() const
{
#ifdef LIBRPA_DEBUG
    this->check_initialized();
#endif
    MPI_Barrier(this->comm);
}

string MPI_COMM_handler::str() const
{
    char s[80];
    sprintf(s, "Proc %4d of Size %4d", this->myid, this->nprocs);
    return string(s);
}

void MPI_COMM_handler::allreduce_matrix(matrix &mat_send,
                                        matrix &mat_recv) const
{
    this->check_initialized();
    MPI_Wrapper::allreduce_matrix(mat_send, mat_recv, this->comm);
}

void MPI_COMM_handler::reduce_matrix(matrix &mat_send, matrix &mat_recv,
                                     int root) const
{
    this->check_initialized();
    MPI_Wrapper::reduce_matrix(mat_send, mat_recv, root, this->comm);
}

void MPI_COMM_handler::allreduce_ComplexMatrix(ComplexMatrix &cmat_sent,
                                               ComplexMatrix &cmat_recv) const
{
    this->check_initialized();
    MPI_Wrapper::allreduce_ComplexMatrix(cmat_sent, cmat_recv, this->comm);
}

void MPI_COMM_handler::reduce_ComplexMatrix(ComplexMatrix &cmat_sent,
                                            ComplexMatrix &cmat_recv,
                                            int root) const
{
    this->check_initialized();
    MPI_Wrapper::reduce_ComplexMatrix(cmat_sent, cmat_recv, root, this->comm);
}

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

void Array_Desc::set_blacs_params_(int ictxt, int nprocs, int myid, int nprows,
                                  int myprow, int npcols, int mypcol)
{
    assert(myid < nprocs && myprow < nprows && mypcol < npcols);
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
    set_blacs_params_(blacs_h.ictxt, blacs_h.nprocs, blacs_h.myid,
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
    int nprocs, myid, nprows, npcols, myprow, mypcol;
    Cblacs_gridinfo(ictxt, &nprows, &npcols, &myprow, &mypcol);
    myid = Cblacs_pnum(ictxt, myprow, mypcol);
    nprocs = nprows * npcols;
    set_blacs_params_(ictxt, nprocs, myid,
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

bool Array_Desc::is_src() const { return myprow_ == irsrc_ && mypcol_ == icsrc_; }

void Array_Desc::barrier(CTXT_SCOPE scope)
{
    CTXT_barrier(ictxt_, scope);
}

int Array_Desc::indx_g2l_r(int gindx) const
{
    return myprow_ != ScalapackConnector::indxg2p(gindx, mb_, myprow_, irsrc_, nprows_) || gindx >= m_
               ? -1
               : ScalapackConnector::indxg2l(gindx, mb_, myprow_, irsrc_, nprows_);
	// int inproc = int((gindx % (mb_*nprows_)) / mb_);
	// if(myprow_==inproc)
	// {
	// 	return int(gindx / (mb_*nprows_))*mb_ + gindx % mb_;
	// }
	// else
	// {
	// 	return -1;
	// }
}

int Array_Desc::indx_g2l_c(int gindx) const
{
    return mypcol_ != ScalapackConnector::indxg2p(gindx, nb_, mypcol_, icsrc_, npcols_) || gindx >= n_
               ? -1
               : ScalapackConnector::indxg2l(gindx, nb_, mypcol_, icsrc_, npcols_);
	// int inproc = int((gindx % (nb_*npcols_)) / mb_);
	// if(mypcol_==inproc)
	// {
	// 	return int(gindx / (nb_*npcols_))*nb_ + gindx % nb_;
	// }
	// else
	// {
	// 	return -1;
	// }
}

int Array_Desc::indx_l2g_r(int lindx) const
{
    return ScalapackConnector::indxl2g(lindx, mb_, myprow_, irsrc_, nprows_);
	// int iblock, gIndex;
	// iblock = lindx / mb_;
	// gIndex = (iblock*nprows_ + myprow_)* mb_ + lindx % mb_;
	// return gIndex;
}

int Array_Desc::indx_l2g_c(int lindx) const
{
    return ScalapackConnector::indxl2g(lindx, nb_, mypcol_, icsrc_, npcols_);
	// int iblock, gIndex;
	// iblock = lindx / nb_;
	// gIndex = (iblock*npcols_ + mypcol_)* nb_ + lindx % nb_;
	// return gIndex;
}

std::pair<Array_Desc, Array_Desc> prepare_array_desc_mr2d_src_and_all(const BLACS_CTXT_handler &ctxt_h, const int &m, const int &n, const int &mb, const int &nb, const int &irsrc, const int &icsrc)
{
    Array_Desc desc_fullblk_on_src(ctxt_h), desc_all_procs(ctxt_h);
    desc_fullblk_on_src.init(m, n, m, n, irsrc, icsrc);
    desc_all_procs.init(m, n, mb, nb, irsrc, icsrc);
    return {desc_fullblk_on_src, desc_all_procs};
}

std::set<std::pair<int, int>> get_necessary_IJ_from_block_2D(const AtomicBasis &atbasis_row, const AtomicBasis &atbasis_col, const Array_Desc& arrdesc)
{
    std::set<std::pair<int, int>> IJs;
    if (arrdesc.m() != atbasis_row.nb_total)
        throw std::invalid_argument("row basis and array descriptor inconsistent");
    if (arrdesc.n() != atbasis_col.nb_total)
        throw std::invalid_argument("col basis and array descriptor inconsistent");

    const auto mlo = arrdesc.m_loc();
    const auto nlo = arrdesc.n_loc();
    size_t glo;
    int I, J;
    for (int ilo = 0; ilo != mlo; ilo++)
        for (int jlo = 0; jlo != nlo; jlo++)
        {
            glo = arrdesc.indx_l2g_r(ilo);
            I = atbasis_row.get_i_atom(glo);
            glo = arrdesc.indx_l2g_c(jlo);
            J = atbasis_col.get_i_atom(glo);
            IJs.insert({I, J});
        }
    return IJs;
}

std::set<std::pair<int, int>> get_necessary_IJ_from_block_2D_sy(const char &uplo, const AtomicBasis &atbasis, const Array_Desc& arrdesc)
{
    std::set<std::pair<int, int>> IJs;
    if (uplo != 'U' && uplo != 'L')
        throw std::invalid_argument("uplo should be U or L");
    if (arrdesc.m() != atbasis.nb_total)
        throw std::invalid_argument("row basis and array descriptor inconsistent");
    if (arrdesc.n() != atbasis.nb_total)
        throw std::invalid_argument("col basis and array descriptor inconsistent");

    const auto mlo = arrdesc.m_loc();
    const auto nlo = arrdesc.n_loc();
    size_t glo;
    int I, J;
    for (int ilo = 0; ilo != mlo; ilo++)
        for (int jlo = 0; jlo != nlo; jlo++)
        {
            glo = arrdesc.indx_l2g_r(ilo);
            I = atbasis.get_i_atom(glo);
            glo = arrdesc.indx_l2g_c(jlo);
            J = atbasis.get_i_atom(glo);
            if (uplo == 'L')
                I > J? IJs.insert({I, J}): IJs.insert({J, I});
            else
                I > J? IJs.insert({J, I}): IJs.insert({I, J});
        }
    return IJs;
}

} // namespace LIBRPA

Parallel_MPI::Parallel_MPI(void)
{
}

Parallel_MPI::~Parallel_MPI(void)
{
    // MPI_Finalize();
    // cout<<"Parallel_MPI Object is being deleted !"<<endl;
}

vector<double> Parallel_MPI::pack_mat(const map<size_t,map<size_t,map<Vector3_Order<int>,shared_ptr<matrix>>>> &Cs_m)
{
    vector<double> pack;
    pack.push_back(-999);		const size_t index_size = pack.size()-1;
    for(const auto &I_pair:Cs_m)
    {
        size_t I=I_pair.first;
        for(const auto &J_pair:I_pair.second )
        {
            size_t J=J_pair.first;
            for(const auto &R_pair:J_pair.second)
            {
                const Vector3_Order<int> &R=R_pair.first;
                const matrix &m= *R_pair.second;
                pack.push_back(I);
                pack.push_back(J);
                pack.push_back(R.x);pack.push_back(R.y);pack.push_back(R.z);
                pack.push_back(m.nr);pack.push_back(m.nc);
                pack.insert(pack.end(), m.c, m.c+m.nr*m.nc);
            }
        }
    }
    pack[index_size] = pack.size() - (index_size+1);
    return pack;
}


// map<size_t,map<size_t,map<Vector3_Order<int>,shared_ptr<matrix>>>>  Parallel_MPI::unpack_mat(vector<double> &pack)
// {
//     map<size_t,map<size_t,map<Vector3_Order<int>,shared_ptr<matrix>>>> Cs;
//     auto ptr=pack.begin();
//     const size_t pack_size=ptr[0];
//     ptr+=1;
//     auto ptr_end=ptr+pack_size;
//     while(ptr<ptr_end)
//     {
//         const size_t I=ptr[0], J=ptr[1];
//         const Vector3_Order<int> R={ptr[2],ptr[3],ptr[4]};
//         const int nr=ptr[5], nc=ptr[6];
//         matrix &m=*Cs[I][J][R];
//         m.create(nr,nc);
//         copy(ptr+7,ptr+7+nr*nc,m.c);
//         ptr+=7+nr*nc;
//     }
//     return Cs;
// }

int Parallel_MPI::globalIndex(int localIndex, int nblk, int nprocs, int myproc)
{
	int iblock, gIndex;
	iblock = localIndex / nblk;
	gIndex = (iblock*nprocs + myproc)*nblk + localIndex % nblk;
	return gIndex;
}


int Parallel_MPI::localIndex(int globalIndex, int nblk, int nprocs, int myproc)
{
	int inproc = int((globalIndex % (nblk*nprocs)) / nblk);
	if(myproc==inproc)
	{
		return int(globalIndex / (nblk*nprocs))*nblk + globalIndex % nblk;
	}
	else
	{
		return -1;
	}
}


vector<int> dispatcher(int ist, int ied, unsigned myid, unsigned size, bool sequential)
{
    vector<int> ilist;
    if (ist >= ied) return ilist;
    assert(size > 0);
    unsigned dist = ied - ist;
    int n = dist / size;
    int extra = dist % size;
    bool has_extra = myid < extra;
    /* LIBRPA::utils::lib_printf("%u %d\n", myid, size); */
    unsigned id;
    for ( int i = 0; i != n + has_extra; i++)
    {
        if (sequential)
            // sequential mode: ist, ist+1, ist+2, ...
            id = i + n * myid + min(extra, int(myid));
        else
            // even mode: ist, ist+size, ist+2*size, ...
            id = size * i + myid;
        ilist.push_back(ist + id);
    }
    return ilist;
}

vector<pair<int, int>> dispatcher(int i1st, int i1ed, int i2st, int i2ed,
                                  unsigned myid, unsigned size, bool sequential, bool favor_1st)
{
    vector<pair<int, int>> ilist;
    assert(size > 0);
    if ( (i1st >= i1ed) || (i2st >= i2ed) ) return ilist;
    unsigned dist1 = i1ed - i1st;
    unsigned dist2 = i2ed - i2st;
    int n = dist1 * dist2 / size;
    int extra = (dist1 * dist2) % size;
    bool has_extra = myid < extra;
    /* LIBRPA::utils::lib_printf("%u %d\n", myid, size); */
    unsigned id, id1, id2;
    for ( int i = 0; i != n + has_extra; i++)
    {
        if (sequential)
            // sequential mode: ist, ist+1, ist+2, ...
            id = i + n * myid + min(extra, int(myid));
        else
            // even mode: ist, ist+size, ist+2*size, ...
            id = size * i + myid;
        if ( favor_1st )
        {
            // id1 goes faster
            id1 = id % dist1;
            id2 = id / dist1;
        }
        else
        {
            // id2 goes faster
            id2 = id % dist2;
            id1 = id / dist2;
        }
        ilist.push_back({i1st+id1, i2st+id2});
    }
    return ilist;
}

vector<pair<int,int>> dispatch_upper_trangular_tasks(const int &natoms, const int &myid, const int &nprows, const int &npcols, const int &myprow, const int &mypcol)
{
    // int myid = blacs_ctxt_world_h.myid;
    // int nprows = blacs_ctxt_world_h.nprows;
    // int npcols = blacs_ctxt_world_h.npcols;
    // int myprow = blacs_ctxt_world_h.myprow;
    // int mypcol = blacs_ctxt_world_h.mypcol;

    int rev_myprow = nprows - 1 - myprow;
    int rev_mypcol = npcols - 1 - mypcol;

    bool flag_former=true;
    if(myprow>rev_myprow)
    {
        flag_former=false;
    }
    else if(myprow == rev_myprow && mypcol> rev_mypcol)
    {
        flag_former=false;
    }

    auto list_row = dispatcher(0, natoms, myprow, nprows, true );
    auto list_col = dispatcher(0, natoms, mypcol, npcols, true );
    auto list_rev_row = dispatcher(0, natoms, rev_myprow, nprows, true );
    auto list_rev_col = dispatcher(0, natoms, rev_mypcol, npcols, true );

    auto loc_task= pick_upper_trangular_tasks(list_row,list_col);
    auto rev_loc_task =pick_upper_trangular_tasks(list_rev_row,list_rev_col);

    vector<pair<int,int>> combine_task, final_loc_task;
    if(myprow == rev_myprow && mypcol==rev_mypcol)
    {
        return loc_task;
    }
    else if(flag_former)
    {
        combine_task.insert(combine_task.end(),loc_task.begin(),loc_task.end());
        combine_task.insert(combine_task.end(),rev_loc_task.begin(),rev_loc_task.end());
        int n_half_task= combine_task.size()/2;
        final_loc_task.insert(final_loc_task.end(),combine_task.begin(),combine_task.begin()+n_half_task);
    }
    else
    {
        combine_task.insert(combine_task.end(),rev_loc_task.begin(),rev_loc_task.end());
        combine_task.insert(combine_task.end(),loc_task.begin(),loc_task.end());
        int n_half_task= combine_task.size()/2;
        final_loc_task.insert(final_loc_task.end(),combine_task.begin()+n_half_task,combine_task.end());
    }
    // for(auto &iap:final_loc_task)
    //     LIBRPA::utils::lib_printf(" loc_task  myid: %d, myprow ,mypcol: %d, %d  task-pair ( %d, %d ) \n",myid, myprow,mypcol, iap.first, iap.second);
    return final_loc_task;
}

vector<pair<int,int>> pick_upper_trangular_tasks(vector<int> list_row, vector<int> list_col)
{
    vector<pair<int,int>> loc_task;
    for(auto &lr:list_row)
        for(auto &lc:list_col)
            if(lr<=lc)
                loc_task.push_back(pair<int,int>(lr,lc));
    return loc_task;
}

vector<pair<int, int>> find_duplicate_ordered_pair(int n, const vector<pair<int, int>>& ordered_pairs, const MPI_Comm &comm)
{
    int myid, nprocs;
    vector<pair<int, int>> pairs_duplicate;
    MPI_Comm_rank(comm, &myid);
    MPI_Comm_size(comm, &nprocs);
    const size_t npairs_total = n * n;
    const size_t maxbytes = 1000 * 1000 * 1000; // 1 GB
    // Communicate atom pair information per batch, reduce memory usage
    const size_t npairs_batch = min(maxbytes / size_t(nprocs), npairs_total);
    const size_t nbatch = npairs_total % npairs_batch ?
                          npairs_total / npairs_batch + 1 :
                          npairs_total / npairs_batch;
    size_t npairs_local = ordered_pairs.size();
    vector<size_t> npairs(nprocs, 0);
    MPI_Allgather(&npairs_local, 1, MPI_UNSIGNED_LONG, npairs.data(), 1, MPI_UNSIGNED_LONG, comm);
    // flatten the pair index, using ordered set
    // NOTE: assuming no duplicates in ordered_pairs
    std::set<size_t> ordered_pairs_flatten;
    for (const auto& op: ordered_pairs)
    {
        ordered_pairs_flatten.insert(op.first * n + op.second);
    }
    for (int ib = 0; ib != nbatch; ib++)
    {
        const size_t pair_start_batch = ib * npairs_batch;
        const size_t npairs_batch_current = ib == nbatch - 1?
            npairs_total - npairs_batch * ib : npairs_batch;
        vector<unsigned char> have_pair_all(npairs_batch_current * nprocs, 0);
        {
            vector<unsigned char> have_pair_this(npairs_batch_current * nprocs, 0);
            for (size_t i = 0; i < npairs_batch_current; i++)
            {
                size_t id_pair = pair_start_batch + i;
                if (std::find(ordered_pairs_flatten.cbegin(), ordered_pairs_flatten.cend(), id_pair) != ordered_pairs_flatten.cend())
                {
                    int ind = nprocs * i + myid;
                    have_pair_this[ind] = 49;
                }
            }
            MPI_Allreduce(have_pair_this.data(), have_pair_all.data(), npairs_batch * nprocs, MPI_UNSIGNED_CHAR, MPI_SUM, comm);
        }
        map<size_t, vector<int>> has_copies;
        for (size_t iap = 0; iap < npairs_batch_current; iap++)
        {
            for (int id = 0; id < nprocs; id++)
            {
                if (have_pair_all[iap * nprocs + id] == 49)
                {
                    has_copies[iap].push_back(id);
                }
            }
        }
        // now check the duplicates
        for (const auto &i_v: has_copies)
        {
            size_t id_pair = pair_start_batch + i_v.first;
            // keep the one with least elements at this stage, remove the other copies (short-sighted)
            // Assuming that there is at least one copy
            // if (i_v.second.size() == 0) cout << "Warning! " << i_v.first << " has no copy across\n";
            if (i_v.second.size() == 1) continue;
            vector<size_t> sizes(i_v.second.size());
            for (int i = 0; i < sizes.size(); i++)
            {
                sizes[i] = npairs[i_v.second[i]];
            }
            auto iter_id = std::min_element(sizes.cbegin(), sizes.cend());
            int id_keep = i_v.second[std::distance(sizes.cbegin(), iter_id)];
            if (myid != id_keep && std::find(i_v.second.cbegin(), i_v.second.cend(), myid) != i_v.second.cend())
            {
                pairs_duplicate.push_back({id_pair / n, id_pair % n});
            }
            for (const auto &id: i_v.second)
            {
                if (id != id_keep) npairs[id] -= 1;
            }
        }
    }
    return pairs_duplicate;
}
