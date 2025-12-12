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
    const int n_aux = ds.basis_aux.nb_total;
    const int m = n_aux, n = n_aux;
    // A column-major 2x2 process grid, all q-points on it
    const int mb = (m + 1) / 2, nb = (n + 1) / 2;
    int lbrow = 0, ubrow = m;
    int lbcol = 0, ubcol = n;
    if (ds.comm_h.myid % 2 == 0) ubrow = mb;
    else lbrow = mb;
    if (ds.comm_h.myid / 2 == 0) ubcol = nb;
    else lbcol = nb;

    ds.vq_lbrow = lbrow;
    ds.vq_ubrow = ubrow;
    ds.vq_lbcol = lbcol;
    ds.vq_ubcol = ubcol;
    ds.pbc.set_kgrids_kvec(2, 2, 2,
                           {
                               0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.5, 0.0, 0.0, 0.5, 0.5,
                               0.5, 0.0, 0.0, 0.5, 0.0, 0.5, 0.5, 0.5, 0.0, 0.5, 0.5, 0.5,
                           });
    for (const auto &k: ds.pbc.klist)
    {
        ds.vq_block_loc[k].create(ubrow - lbrow, ubcol - lbcol);
    }
    ds.initialize_comm_blacs_coul();
    ds.finalize_comm_blacs_coul();
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
