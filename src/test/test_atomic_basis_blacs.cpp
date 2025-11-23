#include "../mpi/envs_blacs.h"
#include "../mpi/envs_mpi.h"
#include "../utils_atomic_basis_blacs.h"

#include "../stl_io_helper.h"
#include "../mpi/base_blacs.h"
#include "testutils.h"

#include <stdexcept>
#include <cassert>
#include <unordered_map>

static void test_ap_to_blacs_global_indices_communicate()
{
    using namespace LIBRPA::envs;

    blacs_ctxt_global_h.set_square_grid(true, LIBRPA::CTXT_LAYOUT::R);
    // Process grid:
    //    0  1
    //    2  3
    LIBRPA::Array_Desc ad(blacs_ctxt_global_h);

    // m = n = 4
    size_t m = 4;
    size_t n = m;
    ad.init_1b1p(m, n, 0, 0);
    assert(ad.initialized());
    assert(ad.mb() == 2);
    assert(ad.nb() == 2);

    LIBRPA::AtomicBasis ab;

    // 2 atoms, atom 0 with 1 basis, atom 1 with 3
    ab.set(std::vector<size_t>{1, 3});
    assert(ab.nb_total == m);
    {
        const auto proc2idlist = LIBRPA::utils::get_communicate_global_ids_list_ap_to_blacs(
            myid_global,
            {{0, {{0, 0}}}, {1, {{0, 1}}}, {2, {{1, 0}}}, {3, {{1, 1}}}},
            ab, ab, ad, true, false);
        // flattened indices with column major and sorted with row index going fastest
        const std::vector<std::unordered_map<int, std::vector<size_t>>> sendlist_ref_all(
            {
                {},
                {{0, {4}}},
                {{0, {1}}},
                {{0, {5}}, {1, {9, 13}}, {2, {6, 7}}},
            }
        );
        const auto &sendlist_ref = sendlist_ref_all[myid_global];
        const std::vector<std::unordered_map<int, std::vector<size_t>>> recvlist_ref_all(
            {
                {{1, {4}}, {2, {1}}, {3, {5}}},
                {{3, {9, 13}}},
                {{3, {6, 7}}},
                {}
            }
        );
        const auto &recvlist_ref = recvlist_ref_all[myid_global];
        assert(equal_map_vector(sendlist_ref, proc2idlist.first));
        assert(equal_map_vector(recvlist_ref, proc2idlist.second));
    }

    // 2 atoms, both have 2 basis functions
    ab.set(std::vector<size_t>{2, 2});
    assert(ab.nb_total == m);
    {
        const auto proc2idlist = LIBRPA::utils::get_communicate_global_ids_list_ap_to_blacs(
            myid_global,
            {{0, {{0, 0}}}, {1, {{0, 1}}}, {2, {{1, 0}}}, {3, {{1, 1}}}},
            ab, ab, ad, true, false);
        // flattened indices with column major and sorted with row index going fastest
        const std::vector<std::unordered_map<int, std::vector<size_t>>> sendlist_ref_all({
            {},
            {},
            {},
            {},
        });
        const auto &sendlist_ref = sendlist_ref_all[myid_global];
        const std::vector<std::unordered_map<int, std::vector<size_t>>> recvlist_ref_all({
            {},
            {},
            {},
            {},
        });
        const auto &recvlist_ref = recvlist_ref_all[myid_global];
        assert(equal_map_vector(sendlist_ref, proc2idlist.first));
        assert(equal_map_vector(recvlist_ref, proc2idlist.second));
    }

    // 2 atoms, atom 0 with 3 and atom 1 with 1
    ab.set(std::vector<size_t>{3, 1});
    assert(ab.nb_total == m);
    {
        const auto proc2idlist = LIBRPA::utils::get_communicate_global_ids_list_ap_to_blacs(
            myid_global,
            {{0, {{0, 0}}}, {1, {{0, 1}}}, {2, {{1, 0}}}, {3, {{1, 1}}}},
            ab, ab, ad, true, false);
        // flattened indices with column major and sorted with row index going fastest
        const std::vector<std::unordered_map<int, std::vector<size_t>>> sendlist_ref_all({
            {{1, {8, 9}}, {2, {2, 6}}, {3, {10}}}, // 1:{(0,2) (1,2)} 2:{(2,0) (2,1)} 3:{2,2}
            {{3, {14}}}, // 3:{(2,3)}
            {{3, {11}}}, // 3:{(3,2)}
            {},
        });
        const auto &sendlist_ref = sendlist_ref_all[myid_global];
        const std::vector<std::unordered_map<int, std::vector<size_t>>> recvlist_ref_all({
            {},
            {{0, {8, 9}}},
            {{0, {2, 6}}},
            {{0, {10}}, {1, {14}}, {2, {11}}},
        });
        const auto &recvlist_ref = recvlist_ref_all[myid_global];
        assert(equal_map_vector(sendlist_ref, proc2idlist.first));
        assert(equal_map_vector(recvlist_ref, proc2idlist.second));
    }

    // 3 atoms, atom 0 and 2 with 1 and atom 1 with 2
    // | 0   0 | 0   1 |
    // | 0   3 | 3   1 |
    // |-------|-------|
    // | 0   3 | 3   1 |
    // | 2   2 | 2   3 |
    ab.set(std::vector<size_t>{1, 2, 1});
    assert(ab.nb_total == m);
    {
        const auto proc2idlist = LIBRPA::utils::get_communicate_global_ids_list_ap_to_blacs(
            myid_global,
            {
             {0, {{0, 0}, {1, 0}, {0, 1}}},
             {1, {{0, 2}, {1, 2}}},
             {2, {{2, 0}, {2, 1}}},
             {3, {{1, 1}, {2, 2}}}
            },
            ab, ab,
            ad, true, false);
        // flattened indices with column major and sorted with row index going fastest
        const std::vector<std::unordered_map<int, std::vector<size_t>>> sendlist_ref_all({
            {{1, {8}}, {2, {2}}}, // 1:{(0,2)} 2:{(2,0)}
            {{3, {14}}}, // 3:{(2,3)}
            {{3, {11}}}, // 3:{(3,2)}
            {{0, {5}}, {1, {9}}, {2, {6}}}, // 0:{(1,1)} 1:{(1,2)} 2:{(2,1)}
        });
        const auto &sendlist_ref = sendlist_ref_all[myid_global];
        const std::vector<std::unordered_map<int, std::vector<size_t>>> recvlist_ref_all({
            {{3, {5}}},
            {{0, {8}}, {3, {9}}},
            {{0, {2}}, {3, {6}}},
            {{1, {14}}, {2, {11}}},
        });
        const auto &recvlist_ref = recvlist_ref_all[myid_global];

        assert(equal_map_vector(sendlist_ref, proc2idlist.first));
        assert(equal_map_vector(recvlist_ref, proc2idlist.second));
    }

    // m = n = 6
    m = 6;
    n = m;
    ad.init_1b1p(m, n, 0, 0);
    assert(ad.initialized());
    assert(ad.mb() == 3);
    assert(ad.nb() == 3);
    ab.set(std::vector<size_t>{1, 2, 1, 2});
    assert(ab.nb_total == m);
    // 4 atoms, atom 0 and 2 with 1 and atom 1 and 3 with 2
    // | 0   1   1 | 2   3   3 |
    // | 0   1   1 | 2   3   3 |
    // | 0   1   1 | 2   3   3 |
    // |-----------|-----------|
    // | 0   1   1 | 2   3   3 |
    // | 0   1   1 | 2   3   3 |
    // | 0   1   1 | 2   3   3 |
    {
        const auto proc2idlist = LIBRPA::utils::get_communicate_global_ids_list_ap_to_blacs(
            myid_global,
            {
             {0, {{0, 0}, {1, 0}, {2, 0}, {3, 0}}},
             {1, {{0, 1}, {1, 1}, {2, 1}, {3, 1}}},
             {2, {{0, 2}, {1, 2}, {2, 2}, {3, 2}}},
             {3, {{0, 3}, {1, 3}, {2, 3}, {3, 3}}},
            },
            ab, ab,
            ad, true, false, false); // column major
        // flattened indices with column major and sorted with row index going fastest
        // const std::vector<std::unordered_map<int, std::vector<size_t>>> sendlist_ref_all({
        //     {{2, {3, 4, 5}}},
        //     {{0, {6, 12, 7, 8, 13, 14}}, {2, {9, 15, 10, 11, 16, 17}}},
        //     {{1, {18, 19, 20}}, {3, {21, 22, 23}}},
        //     {{1, {24, 30, 25, 26, 31, 32}}},
        // });
        // const auto &sendlist_ref = sendlist_ref_all[myid_global];
        // const std::vector<std::unordered_map<int, std::vector<size_t>>> recvlist_ref_all({
        //     {{1, {6, 7, 8, 12, 13, 14}}},
        //     {{2, {18, 19, 20,}}, {3, {24, 25, 26, 30, 31, 32}}},
        //     {{0, {3, 4, 5}}, {1, {9, 10, 11, 15, 16, 17,}}},
        //     {{2, {21, 22, 23}}},
        // });
        // const auto &recvlist_ref = recvlist_ref_all[myid_global];
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
        // assert(equal_map_vector(sendlist_ref, proc2idlist.first));
        // assert(equal_map_vector(recvlist_ref, proc2idlist.second));
    }

    blacs_ctxt_global_h.exit();
}


static void test_ap_to_blacs_local_indices_communicate()
{
    using namespace LIBRPA::envs;

    const size_t m = 4;
    const size_t n = m;
    blacs_ctxt_global_h.set_square_grid(true, LIBRPA::CTXT_LAYOUT::R);
    // Process grid:
    //    0  1
    //    2  3
    LIBRPA::Array_Desc ad(blacs_ctxt_global_h);
    ad.init_1b1p(m, n, 0, 0);
    assert(ad.initialized());
    assert(ad.mb() == 2);
    assert(ad.nb() == 2);

    LIBRPA::AtomicBasis ab;

    // 2 atoms, atom 0 with 1 basis, atom 1 with 3
    ab.set(std::vector<size_t>{1, 3});
    assert(ab.nb_total == m);
    {
        const auto proc2idlist = LIBRPA::utils::get_communicate_local_ids_list_ap_to_blacs(
            myid_global,
            {{0, {{0, 0}}}, {1, {{0, 1}}}, {2, {{1, 0}}}, {3, {{1, 1}}}},
            ab, ab, ad, true, false);
        // flattened indices with column major and sorted with row index going fastest
        const std::vector<std::unordered_map<int, std::vector<std::pair<atpair_t, std::vector<size_t>>>>> sendlist_ref_all(
            {
                {},
                {{0, {{{0, 1}, {0}}}}},
                {{0, {{{1, 0}, {0}}}}},
                {{0, {{{1, 1}, {0}}}}, {1, {{{1, 1}, {3, 6}}}},{2, {{{1, 1}, {1, 2}}}}},
            }
        );
        const auto &sendlist_ref = sendlist_ref_all[myid_global];
        const std::vector<std::unordered_map<int, std::vector<size_t>>> recvlist_ref_all(
            {
                {{1, {2}}, {2, {1}}, {3, {3}}},
                {{3, {1, 3}}},
                {{3, {2, 3}}},
                {}
            }
        );
        const auto &recvlist_ref = recvlist_ref_all[myid_global];
        assert(equal_map_vector_pv(sendlist_ref, proc2idlist.first));
        assert(equal_map_vector(recvlist_ref, proc2idlist.second));
    }

    // 3 atoms, atom 0 and 2 with 1 and atom 1 with 2
    // | 0   0 | 0   1 |
    // | 0   3 | 3   1 |
    // |-------|-------|
    // | 0   3 | 3   1 |
    // | 2   2 | 2   3 |
    ab.set(std::vector<size_t>{1, 2, 1});
    assert(ab.nb_total == m);
    {
        const auto proc2idlist = LIBRPA::utils::get_communicate_local_ids_list_ap_to_blacs(
            myid_global,
            {
             {0, {{0, 0}, {1, 0}, {0, 1}}},
             {1, {{0, 2}, {1, 2}}},
             {2, {{2, 0}, {2, 1}}},
             {3, {{1, 1}, {2, 2}}}
            },
            ab, ab,
            ad, true, false);
        // flattened indices with column major and sorted with row index going fastest
        // const std::vector<std::unordered_map<int, std::vector<size_t>>> sendlist_ref_all({
        //     {{1, {8}}, {2, {2}}}, // 1:{(0,2)} 2:{(2,0)}
        //     {{3, {14}}}, // 3:{(2,3)}
        //     {{3, {11}}}, // 3:{(3,2)}
        //     {{0, {5}}, {1, {9}}, {2, {6}}}, // 0:{(1,1)} 1:{(1,2)} 2:{(2,1)}
        // });
        // flattened indices with column major and sorted with row index going fastest
        const std::vector<std::unordered_map<int, std::vector<std::pair<atpair_t, std::vector<size_t>>>>> sendlist_ref_all(
            {
                {{1, {{{0, 1}, {1}}}}, {2, {{{1, 0}, {1}}}}},
                {{3, {{{1, 2}, {1}}}}},
                {{3, {{{2, 1}, {1}}}}},
                {{0, {{{1, 1}, {0}}}}, {1, {{{1, 1}, {2}}}},{2, {{{1, 1}, {1}}}}},
            }
        );
        const auto &sendlist_ref = sendlist_ref_all[myid_global];
        const std::vector<std::unordered_map<int, std::vector<size_t>>> recvlist_ref_all({
            {{3, {3}}},
            {{0, {0}}, {3, {1}}},
            {{0, {0}}, {3, {2}}},
            {{1, {2}}, {2, {1}}},
        });
        const auto &recvlist_ref = recvlist_ref_all[myid_global];
        assert(equal_map_vector_pv(sendlist_ref, proc2idlist.first));
        assert(equal_map_vector(recvlist_ref, proc2idlist.second));
    }
    blacs_ctxt_global_h.exit();
}


static void test_ap_to_blacs_global_indices_sy_communicate()
{
    using namespace LIBRPA::envs;

    const size_t m = 4;
    const size_t n = m;
    blacs_ctxt_global_h.set_square_grid(true, LIBRPA::CTXT_LAYOUT::R);
    // Process grid:
    //    0  1
    //    2  3
    LIBRPA::Array_Desc ad(blacs_ctxt_global_h);
    ad.init_1b1p(m, n, 0, 0);
    assert(ad.initialized());
    assert(ad.mb() == 2);
    assert(ad.nb() == 2);

    LIBRPA::AtomicBasis ab;

    // 2 atoms, atom 0 with 1 basis, atom 1 with 3
    // | 0   1 | 1   1 |
    // | 2   3 | 3   3 |
    // |-------|-------|
    // | 2   3 | 3   3 |
    // | 2   3 | 3   3 |
    ab.set(std::vector<size_t>{1, 3});
    assert(ab.nb_total == m);
    {
        const auto proc2idlist = LIBRPA::utils::get_communicate_global_ids_list_ap_to_blacs_sy(
            myid_global, 'u',
            {{0, {{0, 0}}}, {1, {{0, 1}}}, {3, {{1, 1}}}},
            ab, ad, true, false);
        // flattened indices with column major and sorted with row index going fastest
        const std::vector<std::unordered_map<int, std::vector<size_t>>> sendlist_ref_all(
            {
                {},
                {{0, {4, 4}}, {2, {8, 12}}},
                {},
                {{0, {5}}, {1, {9, 13}}, {2, {9, 13}}},
            }
        );
        const auto &sendlist_ref = sendlist_ref_all[myid_global];
        const std::vector<std::unordered_map<int, std::pair<std::vector<size_t>, std::vector<char>>>>
            recvlist_ref_all(
            {
                {{1, {{1, 4}, {'c', 'n'}}}, {3, {{5}, {'n'}}}},
                {{3, {{9, 13}, {'n', 'n'}}}},
                {{1, {{2, 3}, {'c', 'c'}}}, {3, {{6, 7}, {'c', 'c'}}}},
                {}
            }
        );
        const auto &recvlist_ref = recvlist_ref_all[myid_global];
        assert(equal_map_vector(sendlist_ref, proc2idlist.first));
        assert(equal_map_pair_vector(recvlist_ref, proc2idlist.second));
    }

    blacs_ctxt_global_h.exit();
}


static void test_ap_to_blacs_local_indices_sy_communicate()
{
    using namespace LIBRPA::envs;

    const size_t m = 4;
    const size_t n = m;
    blacs_ctxt_global_h.set_square_grid(true, LIBRPA::CTXT_LAYOUT::R);

    // Process grid:
    //    0  1
    //    2  3
    LIBRPA::Array_Desc ad(blacs_ctxt_global_h);
    ad.init_1b1p(m, n, 0, 0);
    assert(ad.initialized());
    assert(ad.mb() == 2);
    assert(ad.nb() == 2);

    LIBRPA::AtomicBasis ab;

    // 2 atoms, atom 0 with 1 basis, atom 1 with 3
    // | 0   1 | 1   1 |
    // | 2   3 | 3   3 |
    // |-------|-------|
    // | 2   3 | 3   3 |
    // | 2   3 | 3   3 |
    ab.set(std::vector<size_t>{1, 3});
    assert(ab.nb_total == m);
    {
        const auto proc2idlist = LIBRPA::utils::get_communicate_local_ids_list_ap_to_blacs_sy(
            myid_global, 'u',
            {{0, {{0, 0}}}, {1, {{0, 1}}}, {3, {{1, 1}}}},
            ab, ad, true, false);
        const std::vector<std::unordered_map<int,
                                   std::vector<std::pair<atpair_t, std::vector<size_t>>>>> sendlist_ref_all(
            // global
            // {
            //     {},
            //     {{0, {4, 4}}, {2, {8, 12}}},
            //     {},
            //     {{0, {5}}, {1, {9, 13}}, {2, {9, 13}}},
            // }
            {
                {},
                {{0, {{{0, 1}, {0, 0}}}}, {2, {{{0, 1}, {1, 2}}}}},
                {},
                {{0, {{{1, 1}, {0}}}}, {1, {{{1, 1}, {3, 6}}}}, {2, {{{1, 1}, {3, 6}}}}},
            }
        );
        const auto &sendlist_ref = sendlist_ref_all[myid_global];
        const std::vector<std::unordered_map<int, std::pair<std::vector<size_t>, std::vector<char>>>>
            recvlist_ref_all(
            // global
            // {
            //     {{1, {{1, 4}, {'c', 'n'}}}, {3, {{5}, {'n'}}}},
            //     {{3, {{9, 13}, {'n', 'n'}}}},
            //     {{1, {{2, 3}, {'c', 'c'}}}, {3, {{6, 7}, {'c', 'c'}}}},
            //     {}
            // }
            {
                {{1, {{1, 2}, {'c', 'n'}}}, {3, {{3}, {'n'}}}},
                {{3, {{1, 3}, {'n', 'n'}}}},
                {{1, {{0, 1}, {'c', 'c'}}}, {3, {{2, 3}, {'c', 'c'}}}},
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
        assert(equal_map_pair_vector(recvlist_ref, proc2idlist.second));
    }

    blacs_ctxt_global_h.exit();
}

static void test_blacs_to_ap_global_indices_communicate()
{
    using namespace LIBRPA::envs;

    const size_t m = 4;
    const size_t n = m;
    blacs_ctxt_global_h.set_square_grid(true, LIBRPA::CTXT_LAYOUT::R);
    // Process grid:
    //    0  1
    //    2  3
    LIBRPA::Array_Desc ad(blacs_ctxt_global_h);
    ad.init_1b1p(m, n, 0, 0);
    assert(ad.initialized());
    assert(ad.mb() == 2);
    assert(ad.nb() == 2);

    LIBRPA::AtomicBasis ab;

    // 2 atoms, atom 0 with 1 basis, atom 1 with 3
    // | 0   1 | 1   1 |
    // | 2   3 | 3   3 |
    // |-------|-------|
    // | 2   3 | 3   3 |
    // | 2   3 | 3   3 |
    ab.set(std::vector<size_t>{1, 3});
    assert(ab.nb_total == m);
    {
        const auto proc2idlist = LIBRPA::utils::get_communicate_global_ids_list_blacs_to_ap(
            myid_global,
            {{0, {{0, 0}}}, {1, {{0, 1}}}, {2, {{1, 0}}}, {3, {{1, 1}}}},
            ab, ab, ad, true, false);
        // flattened indices with column major and sorted with row index going fastest
        const std::vector<std::unordered_map<int, std::vector<size_t>>> sendlist_ref_all(
            {
                {{1, {4}}, {2, {1}}, {3, {5}}},
                {{3, {9, 13}}},
                {{3, {6, 7}}},
                {},
            }
        );
        const auto &sendlist_ref = sendlist_ref_all[myid_global];
        const std::vector<std::unordered_map<int, std::vector<size_t>>> recvlist_ref_all(
            {
                {},
                {{0, {4}}},
                {{0, {1}}},
                {{0, {5}}, {1, {9, 13}}, {2, {6, 7}}},
            }
        );
        const auto &recvlist_ref = recvlist_ref_all[myid_global];
        assert(equal_map_vector(sendlist_ref, proc2idlist.first));
        assert(equal_map_vector(recvlist_ref, proc2idlist.second));
    }

    // 2 atoms, both have 2 basis functions
    ab.set(std::vector<size_t>{2, 2});
    assert(ab.nb_total == m);
    {
        const auto proc2idlist = LIBRPA::utils::get_communicate_global_ids_list_blacs_to_ap(
            myid_global,
            {{0, {{0, 0}}}, {1, {{0, 1}}}, {2, {{1, 0}}}, {3, {{1, 1}}}},
            ab, ab, ad, true, false);
        const std::vector<std::unordered_map<int, std::vector<size_t>>> sendlist_ref_all({
            {},
            {},
            {},
            {},
        });
        const auto &sendlist_ref = sendlist_ref_all[myid_global];
        const std::vector<std::unordered_map<int, std::vector<size_t>>> recvlist_ref_all({
            {},
            {},
            {},
            {},
        });
        const auto &recvlist_ref = recvlist_ref_all[myid_global];
        assert(equal_map_vector(sendlist_ref, proc2idlist.first));
        assert(equal_map_vector(recvlist_ref, proc2idlist.second));
    }

    // 2 atoms, atom 0 with 3 and atom 1 with 1
    // | 0   0 | 0   1 |
    // | 0   0 | 0   1 |
    // |-------|-------|
    // | 0   0 | 0   1 |
    // | 2   2 | 2   3 |
    ab.set(std::vector<size_t>{3, 1});
    assert(ab.nb_total == m);
    {
        const auto proc2idlist = LIBRPA::utils::get_communicate_global_ids_list_blacs_to_ap(
            myid_global,
            {{0, {{0, 0}}}, {1, {{0, 1}}}, {2, {{1, 0}}}, {3, {{1, 1}}}},
            ab, ab, ad, true, false);
        // flattened indices with column major and sorted with row index going fastest
        const std::vector<std::unordered_map<int, std::vector<size_t>>> sendlist_ref_all({
            {},
            {{0, {8, 9}}},
            {{0, {2, 6}}},
            {{0, {10}}, {1, {14}}, {2, {11}}},
        });
        const auto &sendlist_ref = sendlist_ref_all[myid_global];
        const std::vector<std::unordered_map<int, std::vector<size_t>>> recvlist_ref_all({
            {{1, {8, 9}}, {2, {2, 6}}, {3, {10}}},
            {{3, {14}}},
            {{3, {11}}},
            {}
        });
        const auto &recvlist_ref = recvlist_ref_all[myid_global];
        assert(equal_map_vector(sendlist_ref, proc2idlist.first));
        assert(equal_map_vector(recvlist_ref, proc2idlist.second));
    }

    // 3 atoms, atom 0 and 2 with 1 and atom 1 with 2
    // | 0   0 | 0   1 |
    // | 0   3 | 3   1 |
    // |-------|-------|
    // | 0   3 | 3   1 |
    // | 2   2 | 2   3 |
    ab.set(std::vector<size_t>{1, 2, 1});
    assert(ab.nb_total == m);
    {
        const auto proc2idlist = LIBRPA::utils::get_communicate_global_ids_list_blacs_to_ap(
            myid_global,
            {
             {0, {{0, 0}, {1, 0}, {0, 1}}},
             {1, {{0, 2}, {1, 2}}},
             {2, {{2, 0}, {2, 1}}},
             {3, {{1, 1}, {2, 2}}}
            },
            ab, ab,
            ad, true, false);
        // flattened indices with column major and sorted with row index going fastest
        const std::vector<std::unordered_map<int, std::vector<size_t>>> sendlist_ref_all({
            {{3, {5}}},
            {{0, {8}}, {3, {9}}},
            {{0, {2}}, {3, {6}}},
            {{1, {14}}, {2, {11}}},
        });
        const auto &sendlist_ref = sendlist_ref_all[myid_global];
        const std::vector<std::unordered_map<int, std::vector<size_t>>> recvlist_ref_all({
            {{1, {8}}, {2, {2}}},
            {{3, {14}}},
            {{3, {11}}},
            {{0, {5}}, {1, {9}}, {2, {6}}},
        });
        const auto &recvlist_ref = recvlist_ref_all[myid_global];

        assert(equal_map_vector(sendlist_ref, proc2idlist.first));
        assert(equal_map_vector(recvlist_ref, proc2idlist.second));
    }
    // same basis, but a different IJ distribtution
    // | 0   1 | 1   2 |
    // | 3   0 | 0   1 |
    // |-------|-------|
    // | 3   0 | 0   1 |
    // | 2   3 | 3   0 |
    {
        const auto proc2idlist = LIBRPA::utils::get_communicate_global_ids_list_blacs_to_ap(
            myid_global,
            {
             {0, {{0, 0}, {1, 1}, {2, 2}}},
             {1, {{0, 1}, {1, 2}}},
             {2, {{0, 2}, {2, 0}}},
             {3, {{1, 0}, {2, 1}}}
            },
            ab, ab,
            ad, true, false);
        // flattened indices with column major and sorted with row index going fastest
        const std::vector<std::unordered_map<int, std::vector<size_t>>> sendlist_ref_all({
            {{1, {4}}, {3, {1}}},
            {{0, {9}}, {2, {12}}},
            {{0, {6}}, {3, {2, 7}}},
            {{1, {14}}, {0, {10, 15}}},
        });
        const auto &sendlist_ref = sendlist_ref_all[myid_global];
        const std::vector<std::unordered_map<int, std::vector<size_t>>> recvlist_ref_all({
            {{1, {9}}, {2, {6}}, {3, {10, 15}}},
            {{3, {14}}, {0, {4}}},
            {{1, {12}}},
            {{0, {1}}, {2, {2, 7}}},
        });
        const auto &recvlist_ref = recvlist_ref_all[myid_global];

        assert(equal_map_vector(sendlist_ref, proc2idlist.first));
        assert(equal_map_vector(recvlist_ref, proc2idlist.second));
    }
    blacs_ctxt_global_h.exit();
}

static void test_blacs_to_ap_local_indices_communicate()
{
    using namespace LIBRPA::envs;

    const size_t m = 4;
    const size_t n = m;
    blacs_ctxt_global_h.set_square_grid(true, LIBRPA::CTXT_LAYOUT::R);
    // Process grid:
    //    0  1
    //    2  3
    LIBRPA::Array_Desc ad(blacs_ctxt_global_h);
    ad.init_1b1p(m, n, 0, 0);
    assert(ad.initialized());
    assert(ad.mb() == 2);
    assert(ad.nb() == 2);

    LIBRPA::AtomicBasis ab;

    // 2 atoms, atom 0 with 1 basis, atom 1 with 3
    ab.set(std::vector<size_t>{1, 3});
    assert(ab.nb_total == m);
    {
        const auto proc2idlist = LIBRPA::utils::get_communicate_local_ids_list_blacs_to_ap(
            myid_global,
            {{0, {{0, 0}}}, {1, {{0, 1}}}, {2, {{1, 0}}}, {3, {{1, 1}}}},
            ab, ab, ad, true, false);
        // flattened indices with column major and sorted with row index going fastest
        const std::vector<std::unordered_map<int, std::vector<std::pair<atpair_t, std::vector<size_t>>>>> recvlist_ref_all(
            {
                {},
                {{0, {{{0, 1}, {0}}}}},
                {{0, {{{1, 0}, {0}}}}},
                {{0, {{{1, 1}, {0}}}}, {1, {{{1, 1}, {3, 6}}}},{2, {{{1, 1}, {1, 2}}}}},
            }
        );
        const std::vector<std::unordered_map<int, std::vector<size_t>>> sendlist_ref_all(
            {
                {{1, {2}}, {2, {1}}, {3, {3}}},
                {{3, {1, 3}}},
                {{3, {2, 3}}},
                {}
            }
        );
        const auto &recvlist_ref = recvlist_ref_all[myid_global];
        const auto &sendlist_ref = sendlist_ref_all[myid_global];
        assert(equal_map_vector(sendlist_ref, proc2idlist.first));
        assert(equal_map_vector_pv(recvlist_ref, proc2idlist.second));
    }

    // 3 atoms, atom 0 and 2 with 1 and atom 1 with 2
    // | 0   0 | 0   1 |
    // | 0   3 | 3   1 |
    // |-------|-------|
    // | 0   3 | 3   1 |
    // | 2   2 | 2   3 |
    ab.set(std::vector<size_t>{1, 2, 1});
    assert(ab.nb_total == m);
    {
        const auto proc2idlist = LIBRPA::utils::get_communicate_local_ids_list_blacs_to_ap(
            myid_global,
            {
             {0, {{0, 0}, {1, 0}, {0, 1}}},
             {1, {{0, 2}, {1, 2}}},
             {2, {{2, 0}, {2, 1}}},
             {3, {{1, 1}, {2, 2}}}
            },
            ab, ab,
            ad, true, false);
        const std::vector<std::unordered_map<int, std::vector<size_t>>> sendlist_ref_all({
            {{3, {3}}},
            {{0, {0}}, {3, {1}}},
            {{0, {0}}, {3, {2}}},
            {{1, {2}}, {2, {1}}},
        });
        const std::vector<std::unordered_map<int, std::vector<std::pair<atpair_t, std::vector<size_t>>>>>
            recvlist_ref_all({
                {{1, {{{0, 1}, {1}}}}, {2, {{{1, 0}, {1}}}}},
                {{3, {{{1, 2}, {1}}}}},
                {{3, {{{2, 1}, {1}}}}},
                {{0, {{{1, 1}, {0}}}}, {1, {{{1, 1}, {2}}}}, {2, {{{1, 1}, {1}}}}},
            });
        const auto &sendlist_ref = sendlist_ref_all[myid_global];
        const auto &recvlist_ref = recvlist_ref_all[myid_global];
        assert(equal_map_vector(sendlist_ref, proc2idlist.first));
        assert(equal_map_vector_pv(recvlist_ref, proc2idlist.second));
    }
    // same basis, but a different IJ distribtution
    // | 0   1 | 1   2 |
    // | 3   0 | 0   1 |
    // |-------|-------|
    // | 3   0 | 0   1 |
    // | 2   3 | 3   0 |
    {
        const auto proc2idlist = LIBRPA::utils::get_communicate_local_ids_list_blacs_to_ap(
            myid_global,
            {
             {0, {{0, 0}, {1, 1}, {2, 2}}},
             {1, {{0, 1}, {1, 2}}},
             {2, {{0, 2}, {2, 0}}},
             {3, {{1, 0}, {2, 1}}}
            },
            ab, ab,
            ad, true, false);
        const std::vector<std::unordered_map<int, std::vector<size_t>>> sendlist_ref_all({
            {{1, {2}}, {3, {1}}},
            {{0, {1}}, {2, {2}}},
            {{0, {2}}, {3, {0, 3}}},
            {{1, {2}}, {0, {0, 3}}},
        });
        const auto &sendlist_ref = sendlist_ref_all[myid_global];
        const std::vector<std::unordered_map<int, std::vector<std::pair<atpair_t, std::vector<size_t>>>>>
            recvlist_ref_all({
                {{1, {{{1, 1}, {2}}}}, {2, {{{1, 1}, {1}}}}, {3, {{{1, 1}, {3}}, {{2, 2}, {0}}}}},
                {{3, {{{1, 2}, {1}}}}, {0, {{{0, 1}, {0}}}}},
                {{1, {{{0, 2}, {0}}}}},
                {{0, {{{1, 0}, {0}}}}, {2, {{{1, 0}, {1}}, {{2, 1}, {0}}}}},
                // global:
                // {{1, {9}}, {2, {6}}, {3, {10, 15}}},
                // {{3, {14}}, {0, {4}}},
                // {{1, {12}}},
                // {{0, {1}}, {2, {2, 7}}},
        });
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
        assert(equal_map_vector_pv(recvlist_ref, proc2idlist.second));
    }

    blacs_ctxt_global_h.exit();
}

static void test_blacs_to_ap_global_indices_sy_communicate()
{
    using namespace LIBRPA::envs;

    const size_t m = 4;
    const size_t n = m;
    blacs_ctxt_global_h.set_square_grid(true, LIBRPA::CTXT_LAYOUT::R);
    // Process grid:
    //    0  1
    //    2  3
    LIBRPA::Array_Desc ad(blacs_ctxt_global_h);
    ad.init_1b1p(m, n, 0, 0);
    assert(ad.initialized());
    assert(ad.mb() == 2);
    assert(ad.nb() == 2);

    LIBRPA::AtomicBasis ab;

    // 2 atoms, atom 0 with 1 basis, atom 1 with 3
    // | 0   1 | 1   1 |
    // | 2   3 | 3   3 |
    // |-------|-------|
    // | 2   3 | 3   3 |
    // | 2   3 | 3   3 |
    ab.set(std::vector<size_t>{1, 3});
    assert(ab.nb_total == m);
    {
        // Using the function for general matrix is sufficient, just only specify either half of atom pairs.
        // This is because BLACS distribution generally contain the full matrix.
        const auto proc2idlist = LIBRPA::utils::get_communicate_global_ids_list_blacs_to_ap(
            myid_global,
            {{0, {{0, 0}}}, {1, {{0, 1}}}, {3, {{1, 1}}}},
            ab, ab, ad, true, false);
        // flattened indices with column major and sorted with row index going fastest
        const std::vector<std::unordered_map<int, std::vector<size_t>>> sendlist_ref_all(
            {
                {{1, {4}}, {3, {5}}},
                {{3, {9, 13}}},
                {{3, {6, 7}}},
                {},
            }
        );
        const auto &sendlist_ref = sendlist_ref_all[myid_global];
        const std::vector<std::unordered_map<int, std::vector<size_t>>>
            recvlist_ref_all(
            {
                {},
                {{0, {4}},},
                {},
                {{0, {5}}, {1, {9, 13}}, {2, {6, 7}}},
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
    }

    blacs_ctxt_global_h.exit();
}

static void test_blacs_to_ap_local_indices_sy_communicate()
{
    using namespace LIBRPA::envs;

    const size_t m = 4;
    const size_t n = m;
    blacs_ctxt_global_h.set_square_grid(true, LIBRPA::CTXT_LAYOUT::R);
    // Process grid:
    //    0  1
    //    2  3
    LIBRPA::Array_Desc ad(blacs_ctxt_global_h);
    ad.init_1b1p(m, n, 0, 0);
    assert(ad.initialized());
    assert(ad.mb() == 2);
    assert(ad.nb() == 2);

    LIBRPA::AtomicBasis ab;

    // 2 atoms, atom 0 with 1 basis, atom 1 with 3
    ab.set(std::vector<size_t>{1, 3});
    assert(ab.nb_total == m);
    {
        const auto proc2idlist = LIBRPA::utils::get_communicate_local_ids_list_blacs_to_ap(
            myid_global,
            {{0, {{0, 0}}}, {1, {{0, 1}}}, {3, {{1, 1}}}},
            ab, ab, ad, true, false);
        // flattened indices with column major and sorted with row index going fastest
        const std::vector<std::unordered_map<int, std::vector<std::pair<atpair_t, std::vector<size_t>>>>> recvlist_ref_all(
            {
                {},
                {{0, {{{0, 1}, {0}}}}},
                {},
                {{0, {{{1, 1}, {0}}}}, {1, {{{1, 1}, {3, 6}}}},{2, {{{1, 1}, {1, 2}}}}},
            }
        );
        const std::vector<std::unordered_map<int, std::vector<size_t>>> sendlist_ref_all(
            {
                {{1, {2}}, {3, {3}}},
                {{3, {1, 3}}},
                {{3, {2, 3}}},
                {}
            }
        );
        const auto &recvlist_ref = recvlist_ref_all[myid_global];
        const auto &sendlist_ref = sendlist_ref_all[myid_global];
        assert(equal_map_vector(sendlist_ref, proc2idlist.first));
        assert(equal_map_vector_pv(recvlist_ref, proc2idlist.second));
    }

    // 3 atoms, atom 0 and 2 with 1 and atom 1 with 2
    // | 0   0 | 0   1 |
    // | 0   3 | 3   1 |
    // |-------|-------|
    // | 0   3 | 3   1 |
    // | 2   2 | 2   3 |
    ab.set(std::vector<size_t>{1, 2, 1});
    assert(ab.nb_total == m);
    {
        const auto proc2idlist = LIBRPA::utils::get_communicate_local_ids_list_blacs_to_ap(
            myid_global,
            {
             {0, {{0, 0}, {0, 1}}},
             {1, {{0, 2}, {1, 2}}},
             {3, {{1, 1}, {2, 2}}}
            },
            ab, ab,
            ad, true, false);
        // flattened indices with column major and sorted with row index going fastest
        // const std::vector<std::unordered_map<int, std::vector<size_t>>> sendlist_ref_all({
        //     {{1, {8}}, {2, {2}}}, // 1:{(0,2)} 2:{(2,0)}
        //     {{3, {14}}}, // 3:{(2,3)}
        //     {{3, {11}}}, // 3:{(3,2)}
        //     {{0, {5}}, {1, {9}}, {2, {6}}}, // 0:{(1,1)} 1:{(1,2)} 2:{(2,1)}
        // });
        // flattened indices with column major and sorted with row index going fastest
        const std::vector<std::unordered_map<int, std::vector<std::pair<atpair_t, std::vector<size_t>>>>> recvlist_ref_all(
            {
                {{1, {{{0, 1}, {1}}}}, },
                {{3, {{{1, 2}, {1}}}}},
                {},
                {{0, {{{1, 1}, {0}}}}, {1, {{{1, 1}, {2}}}},{2, {{{1, 1}, {1}}}}},
            }
        );
        const std::vector<std::unordered_map<int, std::vector<size_t>>> sendlist_ref_all({
            {{3, {3}}},
            {{0, {0}}, {3, {1}}},
            {{3, {2}}},
            {{1, {2}},},
        });
        const auto &sendlist_ref = sendlist_ref_all[myid_global];
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
        assert(equal_map_vector_pv(recvlist_ref, proc2idlist.second));
    }
    blacs_ctxt_global_h.exit();
}

static void test_index_scheduler()
{
    using namespace LIBRPA::envs;

    blacs_ctxt_global_h.set_square_grid(true, LIBRPA::CTXT_LAYOUT::R);
    // Process grid:
    //    0  1
    //    2  3
    LIBRPA::Array_Desc ad(blacs_ctxt_global_h);
    LIBRPA::utils::IndexScheduler sched;
    LIBRPA::AtomicBasis ab;

    size_t m = 4;
    size_t n = m;
    ad.init_1b1p(m, n, 0, 0);
    assert(ad.initialized());
    assert(ad.mb() == 2);
    assert(ad.nb() == 2);

    ab.set(std::vector<size_t>{1, 3});
    assert(ab.nb_total == m);
    {
        const std::unordered_map<int, std::set<atpair_t>> map_proc_IJs
            {{0, {{0, 0}}}, {1, {{0, 1}}}, {2, {{1, 0}}}, {3, {{1, 1}}}};
        sched.init(map_proc_IJs, ab, ab, ad, false); // column major
        const std::vector<std::vector<size_t>> blacs_locid_ref_all({
            {2, 1, 3},
            {1, 3},
            {2, 3},
            {},
        });
        const std::vector<std::vector<size_t>> ap_ipair_ref_all({
            {},
            {0},
            {0},
            {0, 0, 0, 0, 0},
        }  // all zero because there is only one pair on each process in this case
        );
        const std::vector<std::vector<size_t>> ap_locid_ref_all({
            {},
            {0},
            {0},
            {0, 3, 6, 1, 2},
        });

        const auto &blacs_locid_ref = blacs_locid_ref_all[myid_global];
        const auto &ap_ipair_ref = ap_ipair_ref_all[myid_global];
        const auto &ap_locid_ref = ap_locid_ref_all[myid_global];

        assert(equal_vector(blacs_locid_ref, sched.ids_blacs_locid));
        assert(equal_vector(ap_ipair_ref, sched.ids_ap_ipair));
        assert(equal_vector(ap_locid_ref, sched.ids_ap_locid));
    }

    m = 6;
    n = m;
    ad.init_1b1p(m, n, 0, 0);
    assert(ad.initialized());
    assert(ad.mb() == 3);
    assert(ad.nb() == 3);
    ab.set(std::vector<size_t>{1, 2, 1, 2});
    assert(ab.nb_total == m);
    // 4 atoms, atom 0 and 2 with 1 and atom 1 and 3 with 2
    // | 0   1   1 | 2   3   3 |
    // | 0   1   1 | 2   3   3 |
    // | 0   1   1 | 2   3   3 |
    // |-----------|-----------|
    // | 0   1   1 | 2   3   3 |
    // | 0   1   1 | 2   3   3 |
    // | 0   1   1 | 2   3   3 |
    {
        const std::unordered_map<int, std::set<atpair_t>> map_proc_IJs{
            {0, {{0, 0}, {1, 0}, {2, 0}, {3, 0}}},
            {1, {{0, 1}, {1, 1}, {2, 1}, {3, 1}}},
            {2, {{0, 2}, {1, 2}, {2, 2}, {3, 2}}},
            {3, {{0, 3}, {1, 3}, {2, 3}, {3, 3}}},
        };
        sched.init(map_proc_IJs, ab, ab, ad, false); // column major
        const std::vector<std::vector<size_t>> blacs_locid_ref_all({
            {3, 4, 5, 6, 7, 8},
            {0, 1, 2, 3, 4, 5, 6, 7, 8},
            {0, 1, 2, 3, 4, 5, 6, 7, 8},
            {0, 1, 2},
        });
        const std::vector<std::vector<size_t>> ap_ipair_ref_all({
            {2, 3, 3},
            {0, 1, 1, 0, 1, 1, 2, 3, 3, 2, 3, 3},
            {0, 1, 1, 2, 3, 3},
            {0, 1, 1, 0, 1, 1},
        }
        );
        const std::vector<std::vector<size_t>> ap_locid_ref_all({
            {0, 0, 1},
            {0, 0, 1, 1, 2, 3, 0, 0, 1, 1, 2, 3},
            {0, 0, 1, 0, 0, 1},
            {0, 0, 1, 1, 2, 3},
        });

        const auto &blacs_locid_ref = blacs_locid_ref_all[myid_global];
        const auto &ap_ipair_ref = ap_ipair_ref_all[myid_global];
        const auto &ap_locid_ref = ap_locid_ref_all[myid_global];

        // for (int i = 0; i < 4; i++)
        // {
        //     blacs_ctxt_global_h.barrier();
        //     if (myid_global == i)
        //     {
        //         std::cout << "myid " << i << std::endl;
        //         std::cout << "BLACS Ref:" << std::endl << blacs_locid_ref << std::endl;
        //         std::cout << "BLACS Test:" << std::endl << sched.ids_blacs_locid << std::endl;
        //         std::cout << "AP ipair Ref:" << std::endl << ap_ipair_ref << std::endl;
        //         std::cout << "AP ipair Test:" << std::endl << sched.ids_ap_ipair << std::endl;
        //         std::cout << "AP locid Ref:" << std::endl << ap_locid_ref << std::endl;
        //         std::cout << "AP locid Test:" << std::endl << sched.ids_ap_locid << std::endl;
        //     }
        // }
        assert(equal_vector(blacs_locid_ref, sched.ids_blacs_locid));
        assert(equal_vector(ap_ipair_ref, sched.ids_ap_ipair));
        assert(equal_vector(ap_locid_ref, sched.ids_ap_locid));
    }
    blacs_ctxt_global_h.exit();
}

static void test_get_balanced_ap()
{
    using namespace LIBRPA::envs;

    blacs_ctxt_global_h.set_square_grid(true, LIBRPA::CTXT_LAYOUT::R);

    LIBRPA::Array_Desc ad(blacs_ctxt_global_h);
    LIBRPA::AtomicBasis ab;

    {
        ad.init_1b1p(6, 6, 0, 0);
        ab.set(std::vector<size_t>{1, 2, 1, 2});
        const auto map_proc_IJs = LIBRPA::utils::get_balanced_ap_distribution_for_consec_descriptor(ab, ab, ad);
        const std::unordered_map<int, std::vector<atpair_t>> map_proc_IJs_ref = 
        {
            {0, {{0, 0}, {0, 1}, {1, 0}, {1, 1}}},
            {1, {{0, 2}, {0, 3}, {1, 2}, {1, 3}}},
            {2, {{2, 0}, {2, 1}, {3, 0}, {3, 1}}},
            {3, {{2, 2}, {2, 3}, {3, 2}, {3, 3}}},
        };
        assert(map_proc_IJs.size() == map_proc_IJs_ref.size());
        for (const auto &[id, IJs]: map_proc_IJs)
        {
            const auto &IJs_ref = map_proc_IJs_ref.at(id);
            assert(IJs.size() == IJs_ref.size());
            for (const auto &IJ: IJs_ref)
            {
                assert(IJs.find(IJ) != IJs.end());
            }
        }
    }

    // 2-atom system with very different basis size and non-uniform local matrix size
    {
        ad.init_1b1p(11, 11, 0, 0);
        ab.set(std::vector<size_t>{1, 10});
        const auto map_proc_IJs = LIBRPA::utils::get_balanced_ap_distribution_for_consec_descriptor(ab, ab, ad);
        const std::unordered_map<int, std::vector<atpair_t>> map_proc_IJs_ref = 
        {
            {0, {{0, 0}}},
            {1, {{0, 1}}},
            {2, {{1, 0}}},
            {3, {{1, 1}}},
        };
        assert(map_proc_IJs.size() == map_proc_IJs_ref.size());
        for (const auto &[id, IJs]: map_proc_IJs)
        {
            const auto &IJs_ref = map_proc_IJs_ref.at(id);
            assert(IJs.size() == IJs_ref.size());
            for (const auto &IJ: IJs_ref)
            {
                assert(IJs.find(IJ) != IJs.end());
            }
        }
    }

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
    test_ap_to_blacs_global_indices_communicate();
    test_ap_to_blacs_local_indices_communicate();
    test_ap_to_blacs_global_indices_sy_communicate();
    test_ap_to_blacs_local_indices_sy_communicate();

    test_index_scheduler();

    test_blacs_to_ap_global_indices_communicate();
    test_blacs_to_ap_local_indices_communicate();
    test_blacs_to_ap_global_indices_sy_communicate();
    test_blacs_to_ap_local_indices_sy_communicate();

    test_get_balanced_ap();

    // test functions end

    finalize_blacs();
    finalize_mpi();
    MPI_Finalize();

    return 0;
}
