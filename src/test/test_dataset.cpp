#include "../api/dataset.h"
#include "../mpi/base_mpi.h"

#include "../mpi/global_mpi.h"
#include "../io/global_io.h"

static void test_set_comm_blacs_coul_single_blacs()
{
    using namespace librpa_int;

    std::vector<size_t> nbs{2, 3};
    Dataset ds(MPI_COMM_WORLD);
    ds.basis_aux.set(nbs);
    ds.pbc.set_kgrids_kvec(2, 2, 2,
                           {
                               0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.5, 0.0, 0.0, 0.5, 0.5,
                               0.5, 0.0, 0.0, 0.5, 0.0, 0.5, 0.5, 0.5, 0.0, 0.5, 0.5, 0.5,
                           });
    const int n_aux = ds.basis_aux.nb_total;
    const int m = n_aux, n = n_aux;
    const int mb = (m + 1) / 2, nb = (n + 1) / 2;
    const int nkpts = ds.pbc.klist.size();
    int lbrow, ubrow, lbcol, ubcol;

    // A column-major 2x2 process grid, all q-points on it
    lbrow = 0, ubrow = m;
    lbcol = 0, ubcol = n;
    if (ds.comm_h.myid % 2 == 0) ubrow = mb;
    else lbrow = mb;
    if (ds.comm_h.myid / 2 == 0) ubcol = nb;
    else lbcol = nb;
    ds.vq_lbrow = lbrow;
    ds.vq_ubrow = ubrow;
    ds.vq_lbcol = lbcol;
    ds.vq_ubcol = ubcol;
    for (const auto &k: ds.pbc.klist)
    {
        ds.vq_block_loc[k].create(ubrow - lbrow, ubcol - lbcol);
    }
    ds.initialize_comm_blacs_coul();
    assert(ds.blacs_coul_intra_q_h.layout == CTXT_LAYOUT::C);
    assert(ds.blacs_coul_intra_q_h.nprows == 2);
    assert(ds.blacs_coul_intra_q_h.npcols == 2);
    assert(ds.comm_coul_inter_q_h.nprocs == 1);
    ds.finalize_comm_blacs_coul();
    ds.vq_block_loc.clear();

    // A row-major 2x2 process grid, all q-points on it
    lbrow = 0, ubrow = m;
    lbcol = 0, ubcol = n;
    if (ds.comm_h.myid / 2 == 0) ubrow = mb;
    else lbrow = mb;
    if (ds.comm_h.myid % 2 == 0) ubcol = nb;
    else lbcol = nb;
    ds.vq_lbrow = lbrow;
    ds.vq_ubrow = ubrow;
    ds.vq_lbcol = lbcol;
    ds.vq_ubcol = ubcol;
    for (const auto &k: ds.pbc.klist)
    {
        ds.vq_block_loc[k].create(ubrow - lbrow, ubcol - lbcol);
    }
    ds.initialize_comm_blacs_coul();
    assert(ds.blacs_coul_intra_q_h.layout == CTXT_LAYOUT::R);
    assert(ds.blacs_coul_intra_q_h.nprows == 2);
    assert(ds.blacs_coul_intra_q_h.npcols == 2);
    assert(ds.comm_coul_inter_q_h.nprocs == 1);
    ds.finalize_comm_blacs_coul();
    ds.vq_block_loc.clear();

    // 1x1 process grid (major never mind), q-points distributed on 4 tasks
    lbrow = 0, ubrow = m;
    lbcol = 0, ubcol = n;
    ds.vq_lbrow = lbrow;
    ds.vq_ubrow = ubrow;
    ds.vq_lbcol = lbcol;
    ds.vq_ubcol = ubcol;
    for (int ik = 0; ik < nkpts; ik++)
    {
        if (ik % ds.comm_h.nprocs == ds.comm_h.myid)
        {
            const auto &k = ds.pbc.klist[ik];
            ds.vq_block_loc[k].create(ubrow - lbrow, ubcol - lbcol);
        }
    }
    ds.initialize_comm_blacs_coul();
    assert(ds.blacs_coul_intra_q_h.nprows == 1);
    assert(ds.blacs_coul_intra_q_h.npcols == 1);
    assert(ds.comm_coul_inter_q_h.nprocs == 4);
    ds.finalize_comm_blacs_coul();
    ds.vq_block_loc.clear();

    ds.free();
}

int main (int argc, char *argv[])
{
    using namespace librpa_int;
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

    int size_global = get_mpi_size(MPI_COMM_WORLD);
    if (size_global != 4) throw std::runtime_error("test imposes 4 MPI processes");

    global::init_global_mpi();
    global::init_global_io();

    test_set_comm_blacs_coul_single_blacs();

    global::finalize_global_io();
    global::finalize_global_mpi();

    MPI_Finalize();
    return 0;
}
