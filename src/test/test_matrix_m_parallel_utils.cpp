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
    // size 1 matrices to avoid nullptr in pgemr2d
    matrix_m<T> m1(1, 1, MAJOR::COL), m2(1, 1, MAJOR::COL), prod(1, 1, MAJOR::COL), prod_lapack(1, 1, MAJOR::COL);
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
    matrix_m<T> m1_local = init_local_mat<T>(pair_desc_m1.second, MAJOR::COL);
    matrix_m<T> m2_local = init_local_mat<T>(pair_desc_m2.second, MAJOR::COL);
    printf("mat1_local addr %p\n%s", m1_local.ptr(), str(m1_local).c_str());
    printf("mat2_local addr %p\n%s", m2_local.ptr(), str(m2_local).c_str());
    // distribute data
    ScalapackConnector::pgemr2d_f(m, k,
                                  m1.ptr(), 1, 1, pair_desc_m1.first.desc,
                                  m1_local.ptr(), 1, 1, pair_desc_m1.second.desc,
                                  blacs_ctxt_world_h.ictxt);
    ScalapackConnector::pgemr2d_f(k, n,
                                  m2.ptr(), 1, 1, pair_desc_m2.first.desc,
                                  m2_local.ptr(), 1, 1, pair_desc_m2.second.desc,
                                  blacs_ctxt_world_h.ictxt);

    // initialize product matrix
    auto pair_desc_prod = prepare_array_desc_mr2d_src_and_all(blacs_ctxt_world_h, m, n, mb, nb, irsrc, icsrc);
    auto prod_local = init_local_mat<type>(pair_desc_prod.second, MAJOR::COL);

    // carry out distributed multiplication
    ScalapackConnector::pgemm_f('N', 'N', m, n, k, 1.0,
                                m1_local.ptr(), 1, 1, pair_desc_m1.second.desc,
                                m2_local.ptr(), 1, 1, pair_desc_m2.second.desc,
                                0.0,
                                prod_local.ptr(), 1, 1, pair_desc_prod.second.desc);
    printf("prod_local on proc %d:\n%s", blacs_ctxt_world_h.myid, str(prod_local).c_str());
    // collect data back to source
    ScalapackConnector::pgemr2d_f(m, n,
                                  prod_local.ptr(), 1, 1, pair_desc_prod.second.desc,
                                  prod.ptr(), 1, 1, pair_desc_prod.first.desc,
                                  blacs_ctxt_world_h.ictxt);

    if (blacs_ctxt_world_h.myid == pid_src)
    {
        T thres = 1e-14;
        // small threshold when float is used
        if (std::is_same<typename to_real<type>::type, float>::value)
            thres = 1e-6;
        assert(fequal_array(m*n, prod_lapack.ptr(), prod.ptr(), false, thres));
    }

    blacs_ctxt_world_h.exit();
}

template <typename T>
void test_invert_scalapack()
{
    typedef T type;
    typedef typename to_real<type>::type real_type;
    // a looser threshold than gemm
    const T thres = std::is_same<real_type, float>::value ? 1e-4 : 1e-12;

    blacs_ctxt_world_h.set_square_grid();
    assert(blacs_ctxt_world_h.nprocs == 4);
    assert(blacs_ctxt_world_h.nprows == 2);
    assert(blacs_ctxt_world_h.npcols == 2);

    const int n = 8;
    matrix_m<T> mat(1, 1, MAJOR::COL);

    int pid_src, irsrc, icsrc;
    pid_src = 2;
    blacs_ctxt_world_h.get_pcoord(pid_src, irsrc, icsrc);
    // initialize a non-singular matrix at source process
    if (blacs_ctxt_world_h.myid == pid_src)
    {
        bool is_singular = true;
        while (is_singular)
        {
            mat = random<T>(n, n, -1.0, 1.0, MAJOR::COL);
            is_singular = fequal(get_determinant(mat), T(0.0), thres);
        }
    }

    Array_Desc desc_mat(blacs_ctxt_world_h), desc_mat_fb_src(blacs_ctxt_world_h);
    desc_mat.init_1b1p(n, n, irsrc, icsrc);
    desc_mat_fb_src.init(n, n, n, n, irsrc, icsrc);
    // create local matrix and distribute the source to the process grid
    matrix_m<T> mat_loc = init_local_mat<T>(desc_mat, MAJOR::COL);
    ScalapackConnector::pgemr2d_f(n, n,
                                  mat.ptr(), 1, 1, desc_mat_fb_src.desc,
                                  mat_loc.ptr(), 1, 1, desc_mat.desc,
                                  blacs_ctxt_world_h.ictxt);
    auto mat_loc_orig = mat_loc.copy();
    invert_scalapack(mat_loc, desc_mat);
    auto mat_times_invmat = multiply_scalapack(mat_loc, desc_mat, mat_loc_orig, desc_mat, desc_mat);
    ScalapackConnector::pgemr2d_f(n, n,
                                  mat_times_invmat.ptr(), 1, 1, desc_mat.desc,
                                  mat.ptr(), 1, 1, desc_mat_fb_src.desc,
                                  blacs_ctxt_world_h.ictxt);
    printf("mat * invmat on proc %d\n%s", blacs_ctxt_world_h.myid, str(mat_times_invmat).c_str());
    if (blacs_ctxt_world_h.myid == pid_src)
    {
        matrix_m<T> identity(n, n, MAJOR::COL);
        identity.zero_out();
        identity.set_diag(1.0);
        assert(fequal_array(n*n, mat.ptr(), identity.ptr(), false, thres));
    }

    blacs_ctxt_world_h.exit();
}

template <typename T>
void test_power_hemat_blacs_square_grid(const T &m_lb, const T &m_ub)
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
    matrix_m<T> mat(1, 1, MAJOR::COL), mat_gather(1, 1, MAJOR::COL);
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
                                  mat.ptr(), 1, 1, pair_desc_m.first.desc,
                                  mat_loc.ptr(), 1, 1, pair_desc_m.second.desc,
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
                                  mat_loc.ptr(), 1, 1, pair_desc_m.second.desc,
                                  mat_gather.ptr(), 1, 1, pair_desc_m.first.desc,
                                  blacs_ctxt_world_h.ictxt);

    if (blacs_ctxt_world_h.myid == pid_src)
    {
        printf("mat global at pid_src %d\n%s", pid_src, str(mat).c_str());
        printf("mat gathered at pid_src %d\n%s", pid_src, str(mat_gather).c_str());
        assert(fequal_array(n*n, mat.ptr(), mat_gather.ptr(), false, thres));
    }

    delete [] W;

    blacs_ctxt_world_h.exit();
}

template<typename T>
void test_collect_block_from_IJ_storage()
{
    // typedef T type;
    // typedef typename to_real<type>::type real_type;
    blacs_ctxt_world_h.set_square_grid();
    // assert(blacs_ctxt_world_h.nprocs == 4);
    // assert(blacs_ctxt_world_h.nprows == 2);
    // assert(blacs_ctxt_world_h.npcols == 2);
    map<int, map<int, matrix_m<T>>> IJmap;

    AtomicBasis ab(vector<size_t>{1, 2, 2, 1});

    for (int i = 0; i < ab.n_atoms; i++)
    {
        IJmap[i][i].resize(ab.get_atom_nb(i), ab.get_atom_nb(i), MAJOR::COL).randomize(0, 1, !is_complex<T>(), is_complex<T>());
        for (int j = i + 1; j < ab.n_atoms; j++)
        {
            IJmap[i][j].resize(ab.get_atom_nb(i), ab.get_atom_nb(j), MAJOR::COL).randomize(-1, 0, false, false);
            IJmap[j][i] = transpose(IJmap[i][j], is_complex<T>());
        }
    }

    if (blacs_ctxt_world_h.myid == 0)
    {
        for (int I = 0; I < ab.n_atoms; I++)
            for (int J = 0; J < ab.n_atoms; J++)
                cout << I << " " << J << endl << str(IJmap[I][J]);
    }
    blacs_ctxt_world_h.barrier();
    T alpha = 1.0;

    Array_Desc desc_fb(blacs_ctxt_world_h);
    desc_fb.init(ab.nb_total, ab.nb_total, ab.nb_total, ab.nb_total, 0, 0);
    auto mat_loc = init_local_mat<T>(desc_fb, IJmap[0][0].major());

    for (int I = 0; I < ab.n_atoms; I++)
        for (int J = 0; J < ab.n_atoms; J++)
            collect_block_from_IJ_storage(mat_loc, desc_fb, ab, ab, I, J, alpha, IJmap[I][J].ptr(), MAJOR::COL);

    // if (blacs_ctxt_world_h.myid == 0)
    // {
    //     cout << "Use collect_block_from_IJ_storage" << endl << str(mat_loc);
    // }
    blacs_ctxt_world_h.barrier();

    auto mat_loc_syhe = init_local_mat<T>(desc_fb, IJmap[0][0].major());
    for (int I = 0; I < ab.n_atoms; I++)
        for (int J = I; J < ab.n_atoms; J++)
        {
            // printf("myid %d I %d J %d\n", blacs_ctxt_world_h.myid, I, J);
            collect_block_from_IJ_storage_syhe(mat_loc_syhe, desc_fb, ab, I, J, is_complex<T>(), alpha, IJmap[I][J].ptr(), MAJOR::COL);
        }
    // if (blacs_ctxt_world_h.myid == 0)
    // {
    //     cout << "Use collect_block_from_IJ_storage_syhe" << endl << str(mat_loc_syhe);
    // }

    map<int, map<int, matrix_m<T>>> IJmap_collect;
    if (blacs_ctxt_world_h.myid == 0)
    {
        for (int I = 0; I < ab.n_atoms; I++)
            for (int J = 0; J < ab.n_atoms; J++)
                cout << I << " " << J << endl << str(IJmap_collect[I][J]);
    }
    // cout << "mapping the 2D block to IJ" << endl;
    map_block_to_IJ_storage(IJmap_collect, ab, ab, mat_loc_syhe, desc_fb, MAJOR::COL);
    // cout << "done map" << endl;
    if (blacs_ctxt_world_h.myid == 0)
    {
        for (int I = 0; I < ab.n_atoms; I++)
            for (int J = 0; J < ab.n_atoms; J++)
                cout << I << " " << J << endl << str(IJmap_collect[I][J]);
    }

    blacs_ctxt_world_h.exit();
}

int main (int argc, char *argv[])
{
    MPI_Wrapper::init(argc, argv);
    // if ( MPI_Wrapper::nprocs_world != 4 )
    //     throw invalid_argument("test imposes 4 MPI processes");
    mpi_comm_world_h.init();
    blacs_ctxt_world_h.init();

    // test_pgemm<float>(-2, 1);
    // test_pgemm<double>(-2, 1);
    // test_pgemm<complex<float>>({-2, -1}, {1, 0});
    // test_pgemm<complex<double>>({-2, -1}, {1, 0});
    //
    // test_power_hemat_blacs_square_grid<complex<double>>(0.0, {1.0, 1.0});
    // test_power_hemat_blacs_square_grid<complex<float>>(0.0, {1.0, 2.0});

    test_invert_scalapack<float>();
    test_invert_scalapack<double>();
    test_invert_scalapack<complex<float>>();
    test_invert_scalapack<complex<double>>();

    // test_collect_block_from_IJ_storage<double>();
    // test_collect_block_from_IJ_storage<complex<double>>();

    MPI_Wrapper::finalize();

    return 0;
}

