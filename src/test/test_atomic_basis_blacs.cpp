#include "../envs_blacs.h"
#include "../envs_mpi.h"
#include "../utils_atomic_basis_blacs.h"
#include "testutils.h"

#include <stdexcept>
#include <cassert>
#include <unordered_set>

void test_ap_to_2d_indices_communicate()
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
    // assign atom-pairs to each process
    const std::vector<std::vector<atpair_t>> atpairs_all(
        {
            {{0, 0}}, // proc 0
            {{0, 1}}, // proc 1
            {{1, 0}}, // proc 2
            {{1, 1}}, // proc 3
        }
    );
    const auto atpairs = atpairs_all[myid_global];
    // compute task-local basis indices in atom-pair format
    const auto id_ap = LIBRPA::get_1d_indices_in_atpair_blocks(ab, ab, atpairs, true, false);
    const std::vector<std::vector<size_t>> id_ap_ref_all( // reference for column-major
        {
            {0,}, // proc 0
            {4, 8, 12}, // proc 1
            {1, 2, 3}, // proc 2
            {5, 6, 7, 9, 10, 11, 13, 14, 15}, // proc 3
        }
    );
    const auto id_ap_ref = id_ap_ref_all[myid_global];
    assert(id_ap_ref.size() == id_ap.size());
    assert(equal_array(id_ap_ref.size(), id_ap.data(), id_ap_ref.data()));

    // compute task-local basis indices in 2D format
    const auto id_2d = LIBRPA::get_1d_indices_blacs(ad, true, false);
    const std::vector<std::vector<size_t>> id_2d_ref_all( // reference for column-major
        {
            {0, 1, 4, 5}, // proc 0
            {8, 9, 12, 13}, // proc 1
            {2, 3, 6, 7}, // proc 2
            {10, 11, 14, 15}, // proc 3
        }
    );
    const auto id_2d_ref = id_2d_ref_all[myid_global];
    assert(id_2d_ref.size() == id_2d.size());
    assert(equal_array(id_2d_ref.size(), id_2d.data(), id_2d_ref.data()));

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
    test_ap_to_2d_indices_communicate();
    // test functions end

    finalize_blacs();
    finalize_mpi();
    MPI_Finalize();

    return 0;
}
