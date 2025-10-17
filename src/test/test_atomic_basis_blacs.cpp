#include "../envs_blacs.h"
#include "../envs_mpi.h"
#include "../utils_atomic_basis_blacs.h"
#include "testutils.h"

#include <stdexcept>
#include <cassert>

void test_ap_2d_indices_communicate()
{
    using namespace LIBRPA::envs;

    blacs_ctxt_global_h.set_square_grid(true, LIBRPA::CTXT_LAYOUT::R);

    LIBRPA::AtomicBasis ab(std::vector<size_t>{1, 3});
    LIBRPA::Array_Desc ad(blacs_ctxt_global_h);
    ad.init_1b1p(ab.nb_total, ab.nb_total, 0, 0);
    assert(ad.initialized());

    assert(ad.mb() == 2);
    assert(ad.nb() == 2);

    const auto pids = LIBRPA::get_proc_indices_blacs(ad, true);
    const std::vector<int> pids_ref(
            {
                0, 0, 2, 2,
                0, 0, 2, 2,
                1, 1, 3, 3,
                1, 1, 3, 3,
            });
    assert(pids.size() == pids_ref.size());
    assert(equal_array(pids.size(), pids.data(), pids_ref.data()));

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

    // test functions begin
    test_ap_2d_indices_communicate();
    // test functions end

    finalize_blacs();
    finalize_mpi();
    MPI_Finalize();

    return 0;
}
