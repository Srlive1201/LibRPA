#include "parallel_mpi.h"
#include "scalapack_connector.h"
Parallel_MPI::Parallel_MPI(void)
{
    cout<<"Parallel_MPI Object is being created !"<<endl;
    chi_parallel_type = Parallel_MPI::parallel_type::ATOM_PAIR;
}

Parallel_MPI::~Parallel_MPI(void)
{
    MPI_Finalize();
    cout<<"Parallel_MPI Object is being deleted !"<<endl;
}

void Parallel_MPI::allreduce_matrix(matrix &cmat_loc, matrix & cmat_glo)
{
    assert(cmat_loc.nr==cmat_glo.nr);
    assert(cmat_loc.nc==cmat_glo.nc);
    int mat_size=cmat_glo.nr*cmat_glo.nc;
    MPI_Allreduce(cmat_loc.c,cmat_glo.c,mat_size,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
}


void Parallel_MPI::allreduce_ComplexMatrix(ComplexMatrix &cmat_loc, ComplexMatrix & cmat_glo)
{
    assert(cmat_loc.nr==cmat_glo.nr);
    assert(cmat_loc.nc==cmat_glo.nc);
    MPI_Allreduce(cmat_loc.c,cmat_glo.c,cmat_glo.size,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD);
}

void Parallel_MPI::reduce_ComplexMatrix(ComplexMatrix &cmat_loc, ComplexMatrix & cmat_glo)
{
    /* cout << "cmat_loc nr/nc: " << cmat_loc.nr << ", " << cmat_loc.nc << endl; */
    /* cout << "cmat_glo nr/nc: " << cmat_glo.nr << ", " << cmat_glo.nc << endl; */
    assert(cmat_loc.nr==cmat_glo.nr);
    assert(cmat_loc.nc==cmat_glo.nc);
    MPI_Reduce(cmat_loc.c,cmat_glo.c,cmat_glo.size,MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD);
}

void Parallel_MPI::mpi_init(int argc, char **argv)
{
    int provided;
	MPI_Init_thread (&argc, &argv, MPI_THREAD_FUNNELED, &provided);
	if (MPI_THREAD_FUNNELED != provided)
	{
		printf ("%d != required %d", MPI_THREAD_FUNNELED, provided);
    }
    MPI_Comm_rank (MPI_COMM_WORLD, &myid);
    MPI_Comm_size (MPI_COMM_WORLD, &size);
    char name[MPI_MAX_PROCESSOR_NAME];
	int length;
	MPI_Get_processor_name (name, &length);
 
	int omp_num_threads = 1;
	if (argc > 1)
	{
		omp_num_threads = atoi (argv[1]);
	}
	printf ("%s \n", name);
}

void Parallel_MPI::mpi_barrier()
{
    MPI_Barrier(MPI_COMM_WORLD);
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

void Parallel_MPI::set_blacs_parameters()
{
    int np_rows,np_cols; 
    int nprocs=this->size;
    for(np_cols=int(sqrt(double(nprocs))); np_cols>=2; --np_cols)
	{
		if((nprocs)%np_cols==0) break;
	}
	np_rows=nprocs/np_cols;
	my_blacs_ctxt = MPI_Comm_f2c(MPI_COMM_WORLD);
    
	blacs_gridinit_(&my_blacs_ctxt, &BLACS_LAYOUT, &np_rows, &np_cols);
	blacs_gridinfo_(&my_blacs_ctxt, &nprow, &npcol, &myprow, &mypcol);
    cout<<"blacs   myid: " <<myid<<"    np_rows: "<<np_rows<<"   np_cols: "<<np_cols
        <<"   myrow: "<<myprow<<"   mycol: "<<mypcol<<endl;

}

void Parallel_MPI::set_blacs_mat(
    int *desc, int &loc_row, int &loc_col, 
    const int tot_row, const int tot_col,
    const int row_blk, const int col_blk )
{
    int IRSRC=0;
    int info;
    loc_row=ScalapackConnector::numroc(tot_row, row_blk, this->myprow,IRSRC,this->nprow);
    loc_col=ScalapackConnector::numroc(tot_col, row_blk, this->mypcol,IRSRC,this->npcol);
    ScalapackConnector::descinit(desc,tot_row, tot_col, row_blk, col_blk, IRSRC, IRSRC, this->my_blacs_ctxt,loc_row,info);

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
