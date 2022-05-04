#include "parallel_mpi.h"

Parallel_MPI::Parallel_MPI(void)
{
    cout<<"Parallel_MPI Object is being created !"<<endl;
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

Parallel_MPI para_mpi;
