#include "../src/matrix_m_parallel_utils.h"
#include "testutils.h"

using namespace LIBRPA;

template <typename T>
void test_pgemm(const T &m1_lb, const T &m1_ub)
{
    typedef T type;
    blacs_ctxt_world_h.set_square_grid();
    assert(blacs_ctxt_world_h.nprocs == 4);
    assert(blacs_ctxt_world_h.nprows == 2);
    assert(blacs_ctxt_world_h.npcols == 2);

    const int m = 2, n = 6, k = 4;
    const int mb = 1, nb = 3, kb = 2;
    matrix_m<T> m1(0, 0, MAJOR::COL), m2(0, 0, MAJOR::COL), prod(0, 0, MAJOR::COL), prod_lapack(0, 0, MAJOR::COL);
    int pid_src, irsrc, icsrc;
    pid_src = 2;
    blacs_ctxt_world_h.get_pcoord(pid_src, irsrc, icsrc);
    printf("pid_src irsrc icsrc: %d %d %d\n", pid_src, irsrc, icsrc);

    // prepare the real data in the source process
    if (blacs_ctxt_world_h.myid == pid_src)
    {
        m1.resize(m, k);
        m2.resize(k, n);
        prod.resize(m, n);
        m1.randomize(m1_lb, m1_ub);
        m2.randomize();
        printf("Mat1:\n%s", str(m1).c_str());
        printf("Mat2:\n%s", str(m2).c_str());
        prod_lapack = m1 * m2;
        printf("Product:\n%s", str(prod_lapack).c_str());
    }
    blacs_ctxt_world_h.barrier();

    auto pair_desc_m1 = prepare_array_desc_mr2d_src_and_all(blacs_ctxt_world_h, m, k, mb, kb, irsrc, icsrc);
    auto pair_desc_m2 = prepare_array_desc_mr2d_src_and_all(blacs_ctxt_world_h, k, n, kb, nb, irsrc, icsrc);
    auto m1_local = init_local_mat<type>(pair_desc_m1.second, MAJOR::COL);
    auto m2_local = init_local_mat<type>(pair_desc_m2.second, MAJOR::COL);
    // distribute data
    ScalapackConnector::pgemr2d_f(m, k,
                                  m1.c, 1, 1, pair_desc_m1.first.desc,
                                  m1_local.c, 1, 1, pair_desc_m1.second.desc,
                                  blacs_ctxt_world_h.ictxt);
    ScalapackConnector::pgemr2d_f(k, n,
                                  m2.c, 1, 1, pair_desc_m2.first.desc,
                                  m2_local.c, 1, 1, pair_desc_m2.second.desc,
                                  blacs_ctxt_world_h.ictxt);

    // initialize product matrix
    auto pair_desc_prod = prepare_array_desc_mr2d_src_and_all(blacs_ctxt_world_h, m, n, mb, nb, irsrc, icsrc);
    auto prod_local = init_local_mat<type>(pair_desc_prod.second, MAJOR::COL);

    // carry out distributed multiplication
    ScalapackConnector::pgemm_f('N', 'N', m, n, k, 1.0,
                                m1_local.c, 1, 1, pair_desc_m1.second.desc,
                                m2_local.c, 1, 1, pair_desc_m2.second.desc,
                                0.0,
                                prod_local.c, 1, 1, pair_desc_prod.second.desc);
    printf("prod_local on proc %d:\n%s", blacs_ctxt_world_h.myid, str(prod_local).c_str());
    // collect data back to source
    ScalapackConnector::pgemr2d_f(m, n,
                                  prod_local.c, 1, 1, pair_desc_prod.second.desc,
                                  prod.c, 1, 1, pair_desc_prod.first.desc,
                                  blacs_ctxt_world_h.ictxt);

    if (blacs_ctxt_world_h.myid == pid_src)
    {
        T thres = 1e-14;
        // small threshold when float is used
        if (std::is_same<typename to_real<type>::type, float>::value)
            thres = 1e-6;
        assert(fequal_array(m*n, prod_lapack.c, prod.c, false, thres));
    }
    blacs_ctxt_world_h.exit();
}

template <typename T>
void test_power_hemat_blacs(const T &m_lb, const T &m_ub)
{
    typedef T type;
    typedef typename to_real<type>::type real_type;
    const T thres =
        std::is_same<typename to_real<type>::type, float>::value ? 1e-5 : 1e-14;
    blacs_ctxt_world_h.set_square_grid();
    assert(blacs_ctxt_world_h.nprocs == 4);
    assert(blacs_ctxt_world_h.nprows == 2);
    assert(blacs_ctxt_world_h.npcols == 2);

    const int n = 6, nb = n / 2;
    matrix_m<T> mat(0, 0, MAJOR::COL), mat_gather(0, 0, MAJOR::COL);
    int pid_src, irsrc, icsrc;
    pid_src = 2;
    blacs_ctxt_world_h.get_pcoord(pid_src, irsrc, icsrc);
    if (blacs_ctxt_world_h.myid == pid_src)
    {
        vector<real_type> evs {0.1, 0.2, 0.5, 1.0, 2.0, 4.0};
        mat = random_he_selected_ev(n, evs, MAJOR::COL);
        mat_gather = mat;
    }
    auto pair_desc_m = prepare_array_desc_mr2d_src_and_all(blacs_ctxt_world_h, n, n, nb, nb, irsrc, icsrc);
    auto mat_loc = init_local_mat<T>(pair_desc_m.second, MAJOR::COL);
    auto eig_loc = init_local_mat<T>(pair_desc_m.second, MAJOR::COL);
    // printf("eig_loc\n%s", str(eig_loc).c_str());
    ScalapackConnector::pgemr2d_f(n, n,
                                  mat.c, 1, 1, pair_desc_m.first.desc,
                                  mat_loc.c, 1, 1, pair_desc_m.second.desc,
                                  blacs_ctxt_world_h.ictxt);
    // printf("mat_loc of PID %d before\n%s", blacs_ctxt_world_h.myid, str(mat_loc).c_str());
    auto mat_loc_back = mat_loc;

    real_type *W = new real_type [n];

    size_t n_filtered;

    power_hemat_blacs<real_type>(mat_loc, pair_desc_m.second,
                                 eig_loc, pair_desc_m.second, n_filtered, W, 1.0/3.0, -1.0e5);
    assert(n_filtered == 0);
    if (blacs_ctxt_world_h.myid == pid_src)
    {
        printf("Eigenvalues: ");
        for (int i = 0; i < n; i++)
            printf("%f ", W[i]);
        printf("\n");
    }
    blacs_ctxt_world_h.barrier();
    // printf("mat_loc of PID %d middle\n%s", blacs_ctxt_world_h.myid, str(mat_loc).c_str());
    power_hemat_blacs<real_type>(mat_loc, pair_desc_m.second,
                                 eig_loc, pair_desc_m.second, n_filtered, W, 3.0, -1.0e5);
    assert(n_filtered == 0);
    // printf("eig_loc\n%s", str(eig_loc).c_str());
    // printf("mat_loc of PID %d after\n%s", blacs_ctxt_world_h.myid, str(mat_loc).c_str());

    ScalapackConnector::pgemr2d_f(n, n,
                                  mat_loc.c, 1, 1, pair_desc_m.second.desc,
                                  mat_gather.c, 1, 1, pair_desc_m.first.desc,
                                  blacs_ctxt_world_h.ictxt);

    if (blacs_ctxt_world_h.myid == pid_src)
    {
        printf("mat global at pid_src %d\n%s", pid_src, str(mat).c_str());
        printf("mat gathered at pid_src %d\n%s", pid_src, str(mat_gather).c_str());
        assert(fequal_array(n*n, mat.c, mat_gather.c, false, thres));
    }

    delete [] W;

    blacs_ctxt_world_h.exit();
}

int main (int argc, char *argv[])
{
    MPI_Wrapper::init(argc, argv);
    if ( MPI_Wrapper::nprocs_world != 4 )
        throw invalid_argument("test imposes 4 MPI processes");
    mpi_comm_world_h.init();
    blacs_ctxt_world_h.init();

    test_pgemm<double>(-2, 1);
    test_pgemm<complex<double>>({-2, -1}, {1, 0});
    test_pgemm<float>(-2, 1);
    test_pgemm<complex<float>>({-2, -1}, {1, 0});

    test_power_hemat_blacs<complex<double>>(0.0, {1.0, 1.0});
    test_power_hemat_blacs<complex<float>>(0.0, {1.0, 2.0});

    MPI_Wrapper::finalize();

    return 0;
}

