#include "../core/meanfield_mpi.h"

#include "../mpi/global_mpi.h"

using namespace librpa_int;
using namespace librpa_int::global;

static void test_dmat_cplx_Rs_kpara(int nk, int nb, int nocc)
{
    MeanField mf(1, nk, nb, nb);
}

int main (int argc, char *argv[])
{
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    init_global_mpi();

    if ( size_global != 4 )
        throw std::runtime_error("test imposes 4 MPI processes");

    test_dmat_cplx_Rs_kpara(4, 3, 2);

    finalize_global_mpi();
    MPI_Finalize();

    return 0;
}
