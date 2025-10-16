#include "../base_blacs.h"
#include "../envs_blacs.h"

#include "../envs_mpi.h"

#include <cassert>
#include <stdexcept>

void test_arraydesc()
{
    using namespace LIBRPA::envs;
    blacs_ctxt_global_h.set_square_grid();
    LIBRPA::Array_Desc ad(blacs_ctxt_global_h);
    const int m = 10, n = 10;
    // one-block per process, distribution as even as possible
    ad.init_1b1p(m, n, 0, 0);
    printf("%s\n", ad.info().c_str());
    ad.barrier();
    // check overflow
    printf("r%1d c%1d row(m)=%d\n", blacs_ctxt_global_h.myprow, blacs_ctxt_global_h.mypcol, ad.indx_g2l_r(m));
    printf("r%1d c%1d col(n)=%d\n", blacs_ctxt_global_h.myprow, blacs_ctxt_global_h.mypcol, ad.indx_g2l_c(n));
    assert(ad.indx_g2l_r(m) < 0);
    assert(ad.indx_g2l_c(n) < 0);
    // printf("r%1d c%1d row(m)=%d\n",blacs_ctxt_world_h.myprow, blacs_ctxt_world_h.mypcol, ad.indx_g2l_r(m-1));
    // printf("r%1d c%1d col(n)=%d\n",blacs_ctxt_world_h.myprow, blacs_ctxt_world_h.mypcol, ad.indx_g2l_c(n-1));
    blacs_ctxt_global_h.exit();
}

int main (int argc, char *argv[])
{
    using namespace LIBRPA::envs;
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

    initialize_mpi(MPI_COMM_WORLD);
    initialize_blacs(MPI_COMM_WORLD);

    if ( size_global != 4 )
        throw std::runtime_error("test imposes 4 MPI processes");

    test_arraydesc();

    finalize_blacs();
    finalize_mpi();
    MPI_Finalize();

    return 0;
}
