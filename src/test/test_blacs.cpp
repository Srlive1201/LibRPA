#include "../mpi/base_blacs.h"
#include "../mpi/envs_blacs.h"

#include "../mpi/envs_mpi.h"

#include "testutils.h"

#include <cassert>
#include <stdexcept>

void test_proc_indices()
{
    using namespace LIBRPA::envs;

    const int m = 4, n = 6;
    const int mb = m / 2, nb = n / 2;
    const int irsrc = 0, icsrc = 0;

    // row-major 2x2 process grid:
    //  0 1
    //  2 3
    blacs_ctxt_global_h.set_square_grid(true, LIBRPA::CTXT_LAYOUT::R);
    {
        const auto pids_rr = LIBRPA::get_proc_indices_blacs(m, n, mb, nb, irsrc, icsrc, blacs_ctxt_global_h.ictxt, true);
        assert(pids_rr.size() == m * n);
        const std::vector<int> pids_rr_ref(
                {0, 0, 2, 2,
                 0, 0, 2, 2,
                 0, 0, 2, 2,
                 1, 1, 3, 3,
                 1, 1, 3, 3,
                 1, 1, 3, 3,});
        assert(equal_array(m * n, pids_rr.data(), pids_rr_ref.data())); // add myid_global == 0 to print
    }
    blacs_ctxt_global_h.exit();

    // col-major 2x2 process grid:
    //  0 2
    //  1 3
    blacs_ctxt_global_h.set_square_grid(true, LIBRPA::CTXT_LAYOUT::C);
    {
        const auto pids_cc = LIBRPA::get_proc_indices_blacs(m, n, mb, nb, irsrc, icsrc, blacs_ctxt_global_h.ictxt, false);
        assert(pids_cc.size() == m * n);
        const std::vector<int> pids_cc_ref(
                {0, 0, 0, 1, 1, 1,
                 0, 0, 0, 1, 1, 1,
                 2, 2, 2, 3, 3, 3,
                 2, 2, 2, 3, 3, 3});
        assert(equal_array(m * n, pids_cc.data(), pids_cc_ref.data())); // add myid_global == 0 to print
    }
    blacs_ctxt_global_h.exit();
}

void test_arraydesc()
{
    using namespace LIBRPA::envs;
    blacs_ctxt_global_h.set_square_grid();
    LIBRPA::Array_Desc ad(blacs_ctxt_global_h);
    const int m = 10, n = 10;
    // one-block per process, distribution as even as possible
    ad.init_1b1p(m, n, 0, 0);
    assert(ad.initialized());
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

    test_proc_indices();
    test_arraydesc();

    finalize_blacs();
    finalize_mpi();
    MPI_Finalize();

    return 0;
}
