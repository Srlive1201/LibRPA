#include "parallel_mpi.h"

#include <stdexcept>

#include "interface/blacs_scalapack.h"
#include "scalapack_connector.h"

namespace LIBRPA {

parallel_type chi_parallel_type = parallel_type::ATOM_PAIR;

void set_chi_parallel_type(const int &atpais_num, const int Rt_num,
                           const bool use_libri)
{
    if (MPI_Wrapper::is_root_world())
    {
        cout << "Chi parallel type" << endl;
        cout << "| Atom_pairs_num:  " << atpais_num << endl;
        cout << "| R_tau_num:  " << Rt_num << endl;
#ifdef __USE_LIBRI
        cout << "| USE_LibRI:  " << boolalpha << use_libri << endl;
#endif
    }
    if (use_libri)
    {
#ifdef __USE_LIBRI
        chi_parallel_type = parallel_type::LIBRI_USED;
        if (MPI_Wrapper::is_root_world())
            cout << "| Use LibRI for chi0" << endl;
#else
        cout << "LibRI routing requested, but the executable is not compiled "
                "with LibRI"
             << endl;
        cout << "Please recompiler libRPA with -DUSE_LIBRI and configure "
                "include path"
             << endl;
        MPI_Wrapper::barrier_world();
        throw std::logic_error("compilation");
#endif
    }
    else if (atpais_num < Rt_num)
    {
        chi_parallel_type = parallel_type::R_TAU;
        if (MPI_Wrapper::is_root_world()) cout << "| R_tau_routing" << endl;
    }
    else
    {
        chi_parallel_type = parallel_type::ATOM_PAIR;
        if (MPI_Wrapper::is_root_world()) cout << "| atom_pair_routing" << endl;
    }
}

namespace MPI_Wrapper {

bool initialized = false;
std::string procname;
int myid_world = -1;
int nprocs_world = -1;

bool is_root_world() { return PROCID_ROOT == myid_world; }

void init(int argc, char **argv)
{
    if (initialized) return;
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    if (MPI_THREAD_MULTIPLE != provided)
    {
        printf ("Warning: MPI_Init_thread provide %d != required %d", provided, MPI_THREAD_MULTIPLE);
    }
    MPI_Comm_rank(MPI_COMM_WORLD, &myid_world);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs_world);
    char name[MPI_MAX_PROCESSOR_NAME];
    int length;
    MPI_Get_processor_name (name, &length);
    procname = name;
    initialized = true;
}

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

void finalize() { MPI_Finalize(); }

void barrier_world() { MPI_Barrier(MPI_COMM_WORLD); }

}  // namespace MPI_Wraper

MPI_COMM_handler::MPI_COMM_handler(MPI_Comm comm_in)
        : comm(comm_in)
{
    this->initialized = false; 
}

void MPI_COMM_handler::check_initialized() const
{
    if (!initialized)
        throw std::logic_error("MPI_COMM_handler not initialized");
}

void MPI_COMM_handler::init()
{
    MPI_Comm_rank(this->comm, &(this->myid));
    MPI_Comm_size(this->comm, &(this->nprocs));
    this->initialized = true;
}

void MPI_COMM_handler::barrier() const
{
    this->check_initialized();
    MPI_Barrier(this->comm);
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

MPI_COMM_handler mpi_comm_world_h(MPI_COMM_WORLD);

void BLACS_CTXT_handler::init()
{
    this->mpi_comm_h.init();
    this->ictxt = Csys2blacs_handle(this->mpi_comm_h.comm);
    Cblacs_pinfo(&this->myid, &this->nprocs);
    this->initialized = true;
}

void BLACS_CTXT_handler::set_grid(const int &nprows_in, const int &npcols_in,
                                  LAYOUT layout_in)
{
    // if the grid has been set, exit it first
    if (pgrid_set) exit();
    if (nprocs != nprows_in * npcols_in)
        throw std::invalid_argument("nprocs != nprows * npcols");
    layout = layout_in;
    if (layout == LAYOUT::C)
        layout_ch = 'C';
    else
        layout_ch = 'R';
    Cblacs_gridinit(&ictxt, &layout_ch, nprows_in, npcols_in);
    Cblacs_gridinfo(ictxt, &nprows, &npcols, &myprow, &mypcol);
    pgrid_set = true;
}

void BLACS_CTXT_handler::set_square_grid(bool more_rows, LAYOUT layout_in)
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
    set_grid(1, nprocs, LAYOUT::R);
}

void BLACS_CTXT_handler::set_vertical_grid()
{
    set_grid(nprocs, 1, LAYOUT::C);
}

void BLACS_CTXT_handler::exit()
{
    if (pgrid_set)
    {
        Cblacs_gridexit(ictxt);
        // recollect the system context
        ictxt = Csys2blacs_handle(mpi_comm_h.comm);
        pgrid_set = false;
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

void BLACS_CTXT_handler::barrier(SCOPE scope) const
{
    char scope_ch;
    switch (scope)
    {
        case (SCOPE::R): scope_ch = 'R';
        case (SCOPE::C): scope_ch = 'C';
        case (SCOPE::A): scope_ch = 'A';
    }
    Cblacs_barrier(ictxt, &scope_ch);
}

BLACS_CTXT_handler blacs_ctxt_world_h(MPI_COMM_WORLD);

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
        printf(
            "ERROR DESCINIT! PROC %d (%d,%d) PARAMS: DESC %d %d %d %d %d %d %d %d\n",
            myid_, myprow_, mypcol_, m, n, mb, nb, irsrc, icsrc, ictxt_, m_local_);
    // else
    //     printf("SUCCE DESCINIT! PROC %d (%d,%d) PARAMS: DESC %d %d %d %d %d %d %d %d\n", myid_, myprow_, mypcol_, m, n, mb, nb, irsrc, icsrc, ictxt_, m_local_);
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

int Array_Desc::indx_g2l_r(int gindx) const
{
    return myprow_ != ScalapackConnector::indxg2p(gindx, mb_, myprow_, irsrc_, nprows_)
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
    return mypcol_ != ScalapackConnector::indxg2p(gindx, nb_, mypcol_, icsrc_, npcols_)
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
    // cout<<"Parallel_MPI Object is being created !"<<endl;
    chi_parallel_type = Parallel_MPI::parallel_type::ATOM_PAIR;
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


map<size_t,map<size_t,map<Vector3_Order<int>,shared_ptr<matrix>>>>  Parallel_MPI::unpack_mat(vector<double> &pack)
{
    map<size_t,map<size_t,map<Vector3_Order<int>,shared_ptr<matrix>>>> Cs;
    auto ptr=pack.begin();
    const size_t pack_size=ptr[0];
    ptr+=1;
    auto ptr_end=ptr+pack_size;
    while(ptr<ptr_end)
    {
        const size_t I=ptr[0], J=ptr[1];
        const Vector3_Order<int> R={ptr[2],ptr[3],ptr[4]};
        const int nr=ptr[5], nc=ptr[6];
        matrix &m=*Cs[I][J][R];
        m.create(nr,nc);
        copy(ptr+7,ptr+7+nr*nc,m.c);
        ptr+=7+nr*nc;
    }
    return Cs;
}

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
    /* printf("%u %d\n", myid, size); */
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
    /* printf("%u %d\n", myid, size); */
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

Parallel_MPI para_mpi;
