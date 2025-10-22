#include "../envs_blacs.h"
#include "../envs_mpi.h"
#include "../utils_atomic_basis_blacs.h"

#include "../matrix_m.h"
#include "../utils_matrix_m_mpi.h"

#include "../stl_io_helper.h"
#include "testutils.h"

#include <stdexcept>
#include <cassert>

static void test_ap_to_2d_global_indices_communicate()
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


static void test_ap_to_2d_local_indices_communicate()
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

// Example of performing transformation from atom-pair distribution to BLACS disctribution
static void test_ap_to_2d_matrix_m_communicate()
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

    // 3 atoms, atom 0 and 2 with 1 and atom 1 with 2
    // Process indices
    // | 0   0 | 0   1 |
    // | 0   3 | 3   1 |
    // |-------|-------|
    // | 0   3 | 3   1 |
    // | 2   2 | 2   3 |
    ab.set(std::vector<size_t>{1, 2, 1});
    assert(ab.nb_total == m);
    {
        // initialize data
        map<atpair_t, matrix_m<double>> data;
        // Global matrix
        // | -1  |  0   2  |  1 |
        // |-----|---------|----|
        // |  0  |  1   3  |  0 |
        // |  2  |  3  -2  | -1 |
        // |-----|---------|----|
        // |  1  |  0  -1  |  1 |
        if (myid_global == 0)
        {
             data[{0, 0}] = matrix_m<double>({{-1.0}}, MAJOR::COL);
             assert(data.at({0, 0}).nr() == 1 && data.at({0, 0}).nc() == 1);
             data[{1, 0}] = matrix_m<double>({{0.0}, {2.0}}, MAJOR::COL);
             assert(data.at({1, 0}).nr() == 2 && data.at({1, 0}).nc() == 1);
             data[{0, 1}] = matrix_m<double>({{0.0, 2.0}}, MAJOR::COL);
             assert(data.at({0, 1}).nr() == 1 && data.at({0, 1}).nc() == 2);
        }
        else if (myid_global == 1)
        {
             data[{0, 2}] = matrix_m<double>({{1.0}}, MAJOR::COL);
             assert(data.at({0, 2}).nr() == 1 && data.at({0, 2}).nc() == 1);
             data[{1, 2}] = matrix_m<double>({{0.0}, {-1.0}}, MAJOR::COL);
             assert(data.at({1, 2}).nr() == 2 && data.at({1, 2}).nc() == 1);
        }
        else if (myid_global == 2)
        {
             data[{2, 0}] = matrix_m<double>({{1.0}}, MAJOR::COL);
             assert(data.at({2, 0}).nr() == 1 && data.at({2, 0}).nc() == 1);
             data[{2, 1}] = matrix_m<double>({{0.0, -1.0}}, MAJOR::COL);
             assert(data.at({2, 1}).nr() == 1 && data.at({2, 1}).nc() == 2);
        }
        else if (myid_global == 3)
        {
             data[{1, 1}] = matrix_m<double>({{1.0, 3.0}, {3.0, -2.0}}, MAJOR::COL);
             assert(data.at({1, 1}).nr() == 2 && data.at({1, 1}).nc() == 2);
             data[{2, 2}] = matrix_m<double>({{1.0}}, MAJOR::COL);
             assert(data.at({2, 2}).nr() == 1 && data.at({2, 2}).nc() == 1);
        }
        // Reference BLACS local matrix
        // | -1   0  |  2   1 |
        // |  0   1  |  3   0 |
        // |---------|--------|
        // |  2   3  | -2  -1 |
        // |  1   0  | -1   1 |
        auto m_loc_ref = init_local_mat<double>(ad, MAJOR::COL);
        m_loc_ref = 0.0;
        assert(m_loc_ref.nr() == 2 && m_loc_ref.nc() == 2);
        if (myid_global == 0)
        {
            m_loc_ref(0, 0) = -1.0;
            m_loc_ref(1, 1) = 1.0;
        }
        else if (myid_global == 1)
        {
            m_loc_ref(0, 0) = 2.0;
            m_loc_ref(1, 0) = 3.0;
            m_loc_ref(0, 1) = 1.0;
        }
        else if (myid_global == 2)
        {
            m_loc_ref(0, 0) = 2.0;
            m_loc_ref(1, 0) = 1.0;
            m_loc_ref(0, 1) = 3.0;
        }
        else if (myid_global == 3)
        {
            m_loc_ref(0, 0) = -2.0;
            m_loc_ref(1, 0) = -1.0;
            m_loc_ref(0, 1) = -1.0;
            m_loc_ref(1, 1) = 1.0;
        }

        // Initialize BLACS local buffer
        auto m_loc = init_local_mat<double>(ad, MAJOR::COL);
        assert(m_loc.nr() == 2 && m_loc.nc() == 2);
        // Fill in data that is already available before communication
        int I, J, i, j;
        for (int ir = 0; ir < m_loc.nr(); ir++)
        {
            ab.get_local_index(ad.indx_l2g_r(ir), I, i);
            for (int ic = 0; ic < m_loc.nc(); ic++)
            {
                ab.get_local_index(ad.indx_l2g_c(ic), J, j);
                const atpair_t atpair{static_cast<atom_t>(I), static_cast<atom_t>(J)};
                if (data.count(atpair))
                {
                    m_loc(ir, ic) = data.at(atpair)(i, j);
                }
            }
        }

        // compute indices
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

        const auto &pid_ids_send = proc2idlist.first;
        const auto &pid_ids_recv = proc2idlist.second;

        // prepare recv buffer
        // computer recv buffer size and displacements
        map<int, pair<int, MPI_Count>> pid_recv_disp_count;
        std::vector<double> recvbuff(0);
        MPI_Count recvcount = 0;
        for (int pid = 0; pid < size_global; pid++)
        {
            if (pid_ids_recv.count(pid))
            {
                pid_recv_disp_count[pid] = {static_cast<int>(recvcount), pid_ids_recv.at(pid).size()};
                recvcount += pid_ids_recv.at(pid).size();
            }
        }
        if (recvcount > 0) recvbuff.resize(recvcount);

        // prepare send buffer
        // first run, computer send buffer size and displacements, allocate the whole buffer
        map<int, pair<int, MPI_Count>> pid_send_disp_count;
        std::vector<double> sendbuff(0);
        MPI_Count sendcount = 0;
        for (int pid = 0; pid < size_global; pid++)
        {
            if (pid_ids_send.count(pid))
            {
                const auto &ids_send = pid_ids_send.at(pid);
                MPI_Count sendcount_pid = 0;
                for (const auto &pair_ids: ids_send)
                {
                    sendcount_pid += pair_ids.second.size();
                }
                pid_send_disp_count[pid] = {static_cast<int>(sendcount), sendcount_pid};
                sendcount += sendcount_pid;
            }
        }
        if (sendcount > 0) sendbuff.resize(sendcount);
        // second run, fill in the data to be sent
        for (int pid = 0; pid < size_global; pid++)
        {
            if (pid_ids_send.count(pid))
            {
                const auto &ids_send = pid_ids_send.at(pid);
                MPI_Count disp = pid_send_disp_count.at(pid).first;
                for (const auto &pair_ids: ids_send)
                {
                    const auto &atpair = pair_ids.first;
                    const auto &atmat = data.at(atpair);
                    const auto &ids = pair_ids.second;
                    for (auto i = 0; i < ids.size(); i++)
                    {
                        const auto &id = ids[i];
                        sendbuff[disp+i] = atmat.ptr()[id];
                    }
                    disp += ids.size();
                }
            }
        }

        // Begin communication
        std::vector<MPI_Request> reqs;
        MPI_Request req;
        // First non-blocking receive
        for (const auto &pid_disp_count: pid_recv_disp_count)
        {
            const auto &pid = pid_disp_count.first;
            const auto &disp_count = pid_disp_count.second;
            const auto &disp = disp_count.first;
            const auto &count = disp_count.second;
            // ofs << "MPI_Irecv dist " << dist << " count " << count << " from pid " << pid << endl;
            MPI_Irecv(recvbuff.data() + disp, count, MPI_DOUBLE, pid, 0, MPI_COMM_WORLD, &req);
            reqs.push_back(req);
        }
        // Then non-blocking send
        for (const auto &pid_disp_count: pid_send_disp_count)
        {
            const auto &pid = pid_disp_count.first;
            const auto &disp_count = pid_disp_count.second;
            const auto &disp = disp_count.first;
            const auto &count = disp_count.second;
            // ofs << "MPI_Isend count " << count << " to pid " << pid << endl;
            MPI_Isend(sendbuff.data() + disp, count, MPI_DOUBLE, pid, 0, MPI_COMM_WORLD, &req);
            reqs.push_back(req);
        }

        // Wait all non-blocking communication to finish
        if (!reqs.empty()) MPI_Waitall((int)reqs.size(), reqs.data(), MPI_STATUSES_IGNORE);

        // // Intermediate check
        // for (int i = 0; i < 4; i++)
        // {
        //     blacs_ctxt_global_h.barrier();
        //     if (myid_global == i)
        //     {
        //         std::cout << "myid " << i << std::endl;
        //         std::cout << "Send buffer:" << sendbuff << std::endl;
        //         std::cout << "Recv buffer:" << recvbuff << std::endl;
        //     }
        // }

        // Map the received data to the correct location in the BLACS sub-block
        for (const auto &pid_disp_count: pid_recv_disp_count)
        {
            const auto &pid = pid_disp_count.first;
            const auto &ids = pid_ids_recv.at(pid);
            const auto &disp_count = pid_disp_count.second;
            const auto &disp = disp_count.first;
            const auto &count = disp_count.second;
            for (auto i = 0; i < count; i++)
            {
                m_loc.ptr()[ids[i]] = recvbuff[disp+i];
            }
        }

        // Final check
        for (int i = 0; i < 4; i++)
        {
            blacs_ctxt_global_h.barrier();
            if (myid_global == i)
            {
                std::cout << "myid " << i << std::endl;
                assert(fequal_array(m_loc_ref.size(), m_loc_ref.ptr(), m_loc.ptr(), true));
                // fequal_array(m_loc_ref.size(), m_loc_ref.ptr(), m_loc.ptr(), true);
                // std::cout << "Local matrix:" << std::endl << m_loc << std::endl;
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
    test_ap_to_2d_global_indices_communicate();
    test_ap_to_2d_local_indices_communicate();
    test_ap_to_2d_matrix_m_communicate();
    // test functions end

    finalize_blacs();
    finalize_mpi();
    MPI_Finalize();

    return 0;
}
