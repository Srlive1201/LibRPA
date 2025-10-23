#include "../envs_blacs.h"
#include "../envs_mpi.h"
#include "../utils_atomic_basis_blacs.h"

#include "../matrix_m.h"
#include "../utils_matrix_m_mpi.h"

#include "../stl_io_helper.h"
#include "testutils.h"

#include <stdexcept>
#include <cassert>

static void test_ap_to_blacs_global_indices_communicate()
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
        const auto proc2idlist = LIBRPA::utils::get_communicate_global_ids_list_ap_to_blacs(
            myid_global,
            {{0, {{0, 0}}}, {1, {{0, 1}}}, {2, {{1, 0}}}, {3, {{1, 1}}}},
            ab, ab, ad, true, false);
        // flattened indices with column major and sorted with row index going fastest
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
        const std::vector<std::map<int, std::vector<size_t>>> sendlist_ref_all({
            {},
            {},
            {},
            {},
        });
        const auto &sendlist_ref = sendlist_ref_all[myid_global];
        const std::vector<std::map<int, std::vector<size_t>>> recvlist_ref_all({
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
        const std::vector<std::map<int, std::vector<size_t>>> sendlist_ref_all({
            {{1, {8, 9}}, {2, {2, 6}}, {3, {10}}}, // 1:{(0,2) (1,2)} 2:{(2,0) (2,1)} 3:{2,2}
            {{3, {14}}}, // 3:{(2,3)}
            {{3, {11}}}, // 3:{(3,2)}
            {},
        });
        const auto &sendlist_ref = sendlist_ref_all[myid_global];
        const std::vector<std::map<int, std::vector<size_t>>> recvlist_ref_all({
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
        const std::vector<std::map<int, std::vector<size_t>>> sendlist_ref_all({
            {{1, {8}}, {2, {2}}}, // 1:{(0,2)} 2:{(2,0)}
            {{3, {14}}}, // 3:{(2,3)}
            {{3, {11}}}, // 3:{(3,2)}
            {{0, {5}}, {1, {9}}, {2, {6}}}, // 0:{(1,1)} 1:{(1,2)} 2:{(2,1)}
        });
        const auto &sendlist_ref = sendlist_ref_all[myid_global];
        const std::vector<std::map<int, std::vector<size_t>>> recvlist_ref_all({
            {{3, {5}}},
            {{0, {8}}, {3, {9}}},
            {{0, {2}}, {3, {6}}},
            {{1, {14}}, {2, {11}}},
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
        assert(equal_map_vector(recvlist_ref, proc2idlist.second));
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
        const std::vector<std::map<int, std::vector<std::pair<atpair_t, std::vector<size_t>>>>> sendlist_ref_all(
            {
                {},
                {{0, {{{0, 1}, {0}}}}},
                {{0, {{{1, 0}, {0}}}}},
                {{0, {{{1, 1}, {0}}}}, {1, {{{1, 1}, {3, 6}}}},{2, {{{1, 1}, {1, 2}}}}},
            }
        );
        const auto &sendlist_ref = sendlist_ref_all[myid_global];
        const std::vector<std::map<int, std::vector<size_t>>> recvlist_ref_all(
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
        // const std::vector<std::map<int, std::vector<size_t>>> sendlist_ref_all({
        //     {{1, {8}}, {2, {2}}}, // 1:{(0,2)} 2:{(2,0)}
        //     {{3, {14}}}, // 3:{(2,3)}
        //     {{3, {11}}}, // 3:{(3,2)}
        //     {{0, {5}}, {1, {9}}, {2, {6}}}, // 0:{(1,1)} 1:{(1,2)} 2:{(2,1)}
        // });
        // flattened indices with column major and sorted with row index going fastest
        const std::vector<std::map<int, std::vector<std::pair<atpair_t, std::vector<size_t>>>>> sendlist_ref_all(
            {
                {{1, {{{0, 1}, {1}}}}, {2, {{{1, 0}, {1}}}}},
                {{3, {{{1, 2}, {1}}}}},
                {{3, {{{2, 1}, {1}}}}},
                {{0, {{{1, 1}, {0}}}}, {1, {{{1, 1}, {2}}}},{2, {{{1, 1}, {1}}}}},
            }
        );
        const auto &sendlist_ref = sendlist_ref_all[myid_global];
        const std::vector<std::map<int, std::vector<size_t>>> recvlist_ref_all({
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
        const std::vector<std::map<int, std::vector<size_t>>> sendlist_ref_all(
            {
                {},
                {{0, {4, 4}}, {2, {8, 12}}},
                {},
                {{0, {5}}, {1, {9, 13}}, {2, {9, 13}}},
            }
        );
        const auto &sendlist_ref = sendlist_ref_all[myid_global];
        const std::vector<std::map<int, std::pair<std::vector<size_t>, std::vector<char>>>>
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
        const std::vector<std::map<int,
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
        const std::vector<std::map<int, std::pair<std::vector<size_t>, std::vector<char>>>>
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
        const std::vector<std::map<int, std::vector<size_t>>> sendlist_ref_all(
            {
                {{1, {4}}, {2, {1}}, {3, {5}}},
                {{3, {9, 13}}},
                {{3, {6, 7}}},
                {},
            }
        );
        const auto &sendlist_ref = sendlist_ref_all[myid_global];
        const std::vector<std::map<int, std::vector<size_t>>> recvlist_ref_all(
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
        const std::vector<std::map<int, std::vector<size_t>>> sendlist_ref_all({
            {},
            {},
            {},
            {},
        });
        const auto &sendlist_ref = sendlist_ref_all[myid_global];
        const std::vector<std::map<int, std::vector<size_t>>> recvlist_ref_all({
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
        const std::vector<std::map<int, std::vector<size_t>>> sendlist_ref_all({
            {},
            {{0, {8, 9}}},
            {{0, {2, 6}}},
            {{0, {10}}, {1, {14}}, {2, {11}}},
        });
        const auto &sendlist_ref = sendlist_ref_all[myid_global];
        const std::vector<std::map<int, std::vector<size_t>>> recvlist_ref_all({
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
        const std::vector<std::map<int, std::vector<size_t>>> sendlist_ref_all({
            {{3, {5}}},
            {{0, {8}}, {3, {9}}},
            {{0, {2}}, {3, {6}}},
            {{1, {14}}, {2, {11}}},
        });
        const auto &sendlist_ref = sendlist_ref_all[myid_global];
        const std::vector<std::map<int, std::vector<size_t>>> recvlist_ref_all({
            {{1, {8}}, {2, {2}}},
            {{3, {14}}},
            {{3, {11}}},
            {{0, {5}}, {1, {9}}, {2, {6}}},
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
        const std::vector<std::map<int, std::vector<std::pair<atpair_t, std::vector<size_t>>>>> recvlist_ref_all(
            {
                {},
                {{0, {{{0, 1}, {0}}}}},
                {{0, {{{1, 0}, {0}}}}},
                {{0, {{{1, 1}, {0}}}}, {1, {{{1, 1}, {3, 6}}}},{2, {{{1, 1}, {1, 2}}}}},
            }
        );
        const std::vector<std::map<int, std::vector<size_t>>> sendlist_ref_all(
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
        // flattened indices with column major and sorted with row index going fastest
        // const std::vector<std::map<int, std::vector<size_t>>> sendlist_ref_all({
        //     {{1, {8}}, {2, {2}}}, // 1:{(0,2)} 2:{(2,0)}
        //     {{3, {14}}}, // 3:{(2,3)}
        //     {{3, {11}}}, // 3:{(3,2)}
        //     {{0, {5}}, {1, {9}}, {2, {6}}}, // 0:{(1,1)} 1:{(1,2)} 2:{(2,1)}
        // });
        // flattened indices with column major and sorted with row index going fastest
        const std::vector<std::map<int, std::vector<std::pair<atpair_t, std::vector<size_t>>>>> recvlist_ref_all(
            {
                {{1, {{{0, 1}, {1}}}}, {2, {{{1, 0}, {1}}}}},
                {{3, {{{1, 2}, {1}}}}},
                {{3, {{{2, 1}, {1}}}}},
                {{0, {{{1, 1}, {0}}}}, {1, {{{1, 1}, {2}}}},{2, {{{1, 1}, {1}}}}},
            }
        );
        const std::vector<std::map<int, std::vector<size_t>>> sendlist_ref_all({
            {{3, {3}}},
            {{0, {0}}, {3, {1}}},
            {{0, {0}}, {3, {2}}},
            {{1, {2}}, {2, {1}}},
        });
        const auto &sendlist_ref = sendlist_ref_all[myid_global];
        const auto &recvlist_ref = recvlist_ref_all[myid_global];
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
        const std::vector<std::map<int, std::vector<size_t>>> sendlist_ref_all(
            {
                {{1, {4}}, {3, {5}}},
                {{3, {9, 13}}},
                {{3, {6, 7}}},
                {},
            }
        );
        const auto &sendlist_ref = sendlist_ref_all[myid_global];
        const std::vector<std::map<int, std::vector<size_t>>>
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
        const std::vector<std::map<int, std::vector<std::pair<atpair_t, std::vector<size_t>>>>> recvlist_ref_all(
            {
                {},
                {{0, {{{0, 1}, {0}}}}},
                {},
                {{0, {{{1, 1}, {0}}}}, {1, {{{1, 1}, {3, 6}}}},{2, {{{1, 1}, {1, 2}}}}},
            }
        );
        const std::vector<std::map<int, std::vector<size_t>>> sendlist_ref_all(
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
        // const std::vector<std::map<int, std::vector<size_t>>> sendlist_ref_all({
        //     {{1, {8}}, {2, {2}}}, // 1:{(0,2)} 2:{(2,0)}
        //     {{3, {14}}}, // 3:{(2,3)}
        //     {{3, {11}}}, // 3:{(3,2)}
        //     {{0, {5}}, {1, {9}}, {2, {6}}}, // 0:{(1,1)} 1:{(1,2)} 2:{(2,1)}
        // });
        // flattened indices with column major and sorted with row index going fastest
        const std::vector<std::map<int, std::vector<std::pair<atpair_t, std::vector<size_t>>>>> recvlist_ref_all(
            {
                {{1, {{{0, 1}, {1}}}}, },
                {{3, {{{1, 2}, {1}}}}},
                {},
                {{0, {{{1, 1}, {0}}}}, {1, {{{1, 1}, {2}}}},{2, {{{1, 1}, {1}}}}},
            }
        );
        const std::vector<std::map<int, std::vector<size_t>>> sendlist_ref_all({
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

    test_blacs_to_ap_global_indices_communicate();
    test_blacs_to_ap_local_indices_communicate();
    test_blacs_to_ap_global_indices_sy_communicate();
    test_blacs_to_ap_local_indices_sy_communicate();
    // test functions end

    finalize_blacs();
    finalize_mpi();
    MPI_Finalize();

    return 0;
}
