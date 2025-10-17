#include "../envs_blacs.h"
#include "../envs_mpi.h"
#include "../utils_atomic_basis_blacs.h"
#include "testutils.h"
#include "../stl_io_helper.h"

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

    const auto proc2idlist = LIBRPA::utils::get_communicate_ids_list_ap_to_blacs(
        myid_global,
        {{0, {{0, 0}}}, {1, {{0, 1}}}, {2, {{1, 0}}}, {3, {{1, 1}}}},
        ab, ab, ad, true, false);
    const std::vector<std::map<int, std::vector<size_t>>> sendlist_ref_all(
        {
            {},
            {{0, {4}}},
            {{0, {1}}},
            {{0, {5}}, {1, {9, 13}}, {2, {6, 7}}},
        }
    );
    const auto &sendlist_ref = sendlist_ref_all[myid_global];
    const std::vector<std::map<int, std::vector<size_t>>> recvlist_ref_all(
        {
            {{1, {4}}, {2, {1}}, {3, {5}}},
            {{3, {9, 13}}},
            {{3, {6, 7}}},
            {}
        }
    );
    const auto &recvlist_ref = recvlist_ref_all[myid_global];

    // for (int i = 0; i < 4; i++)
    // {
    //     blacs_ctxt_global_h.barrier();
    //     if (myid_global == i)
    //     {
    //         std::cout << "myid " << i << std::endl;
    //         std::cout << "Recv Ref:" << std::endl << recvlist_ref << std::endl;
    //         std::cout << "Recv Test:" << std::endl << proc2idlist.second << std::endl;
    //         std::cout << "Send Ref:" << std::endl << sendlist_ref << std::endl;
    //         std::cout << "Send Test:" << std::endl << proc2idlist.first << std::endl;
    //     }
    // }
    assert(equal_map_vector(sendlist_ref, proc2idlist.first));
    assert(equal_map_vector(recvlist_ref, proc2idlist.second));

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
