#pragma once
#include "matrix_m.h"
#include "constants.h"
#include "parallel_mpi.h"
#include "scalapack_connector.h"

template <typename T>
matrix_m<T> init_local_mat(const LIBRPA::Array_Desc &ad, MAJOR major)
{
    matrix_m<T> mat_lo(ad.m_loc(), ad.n_loc(), major);
    if (mat_lo.size() == 0)
        mat_lo.resize(1, 1);
    return mat_lo;
}

template <typename T>
matrix_m<T> get_local_mat(const matrix_m<T> &mat_go, const LIBRPA::Array_Desc &ad, MAJOR major)
{
    // assert the shape of global matrix conforms with the array descriptor
    assert(mat_go.nr() == ad.m() && mat_go.nc() == ad.n());

    matrix_m<T> mat_lo(ad.m_loc(), ad.n_loc(), major);
    for (int i = 0; i != mat_go.nr(); i++)
    {
        auto i_lo = ad.indx_g2l_r(i);
        if (i_lo < 0) continue;
        for (int j = 0; j != mat_go.nc(); j++)
        {
            auto j_lo = ad.indx_g2l_c(j);
            if (j_lo < 0) continue;
            mat_lo(i_lo, j_lo) = mat_go(i, j);
        }
    }
    if (mat_lo.size() == 0)
        mat_lo.resize(1, 1);
    return mat_lo;
}

template <typename T>
matrix_m<T> get_local_mat(const T *pv, MAJOR major_pv, const LIBRPA::Array_Desc &ad, MAJOR major)
{
    const int nr = ad.m(), nc = ad.n();
    matrix_m<T> mat_lo(ad.m_loc(), ad.n_loc(), major);
    const auto picker = Indx_pickers_2d[major_pv];
    for (int i = 0; i != nr; i++)
    {
        auto i_lo = ad.indx_g2l_r(i);
        if (i_lo < 0) continue;
        for (int j = 0; j != nc; j++)
        {
            auto j_lo = ad.indx_g2l_c(j);
            if (j_lo < 0) continue;
            mat_lo(i_lo, j_lo) = pv[picker(nr, i, nc, j)];
        }
    }
    if (mat_lo.size() == 0)
        mat_lo.resize(1, 1);
    return mat_lo;
}

template <typename Tdst, typename Tsrc>
void collect_block_from_IJ_storage(
    matrix_m<Tdst> &mat_lo,
    const LIBRPA::Array_Desc &ad,
    const LIBRPA::AtomicBasis &atbasis_row,
    const LIBRPA::AtomicBasis &atbasis_col,
    const int &I, const int &J,
    Tdst alpha, const Tsrc *pvIJ, MAJOR major_pv)
{
    // assert(mat_lo.nr() == ad.m_loc() && mat_lo.nc() == ad.n_loc());
    assert(ad.m() == atbasis_row.nb_total && ad.n() == atbasis_col.nb_total );
    const int row_start_id = atbasis_row.get_part_range()[I];
    const int col_start_id = atbasis_col.get_part_range()[J];
    const int row_nb = atbasis_row.get_atom_nb(I);
    const int col_nb = atbasis_col.get_atom_nb(J);
    const auto picker = Indx_pickers_2d[major_pv];
    for (int iI = 0; iI != row_nb; iI++)
    {
        for (int jJ = 0; jJ != col_nb; jJ++)
        {
            int ilo = ad.indx_g2l_r(row_start_id + iI);
            int jlo = ad.indx_g2l_c(col_start_id + jJ);
            if (ilo < 0 || jlo < 0) continue;
            mat_lo(ilo, jlo) += alpha * pvIJ[picker(row_nb, iI, col_nb, jJ)];
        }
    }
}

//! collect 2D block from a IJ-pair storage exploiting its symmetric/Hermitian property
template <typename Tdst, typename Tsrc>
void collect_block_from_IJ_storage_syhe(
    matrix_m<Tdst> &mat_lo,
    const LIBRPA::Array_Desc &ad,
    const LIBRPA::AtomicBasis &atbasis,
    const int &I, const int &J, bool conjugate,
    Tdst alpha, const Tsrc *pvIJ, MAJOR major_pv)
{
    // assert(mat_lo.nr() == ad.m_loc() && mat_lo.nc() == ad.n_loc());
    assert(ad.m() == atbasis.nb_total && ad.n() == atbasis.nb_total);
    // blocks on diagonal is trivial
    if (I == J)
    {
        collect_block_from_IJ_storage(mat_lo, ad, atbasis, atbasis, I, J, alpha, pvIJ, major_pv);
        return;
    }
    const std::function<Tsrc(Tsrc)> filter = conjugate ? [](Tsrc a) { return get_conj(a); }
                                                       : [](Tsrc a) { return a; };
    const int I_start_id = atbasis.get_part_range()[I];
    const int J_start_id = atbasis.get_part_range()[J];
    const int I_nb = atbasis.get_atom_nb(I);
    const int J_nb = atbasis.get_atom_nb(J);
    const auto picker = Indx_pickers_2d[major_pv];
    for (int ilo = 0; ilo != ad.m_loc(); ilo++)
        for (int jlo = 0; jlo != ad.n_loc(); jlo++)
        {
            int i_gl = ad.indx_l2g_r(ilo);
            int j_gl = ad.indx_l2g_c(jlo);
            Tsrc temp;
            if ((i_gl > j_gl && I > J) || (i_gl < j_gl && I < J))
                temp = pvIJ[picker(I_nb, i_gl - I_start_id, J_nb, j_gl - J_start_id)];
            else
                temp = filter(pvIJ[picker(J_nb, i_gl - J_start_id, I_nb, j_gl - I_start_id)]);
            mat_lo(ilo, jlo) += alpha * temp;
        }
}

template <typename T>
matrix_m<std::complex<T>> power_hemat_blacs(matrix_m<std::complex<T>> &A_local,
                                            const LIBRPA::Array_Desc &ad_A,
                                            matrix_m<std::complex<T>> &Z_local,
                                            const LIBRPA::Array_Desc &ad_Z,
                                            size_t &n_filtered, T *W, T power,
                                            const T &threshold = -1.e5)
{
    assert (A_local.is_col_major() && Z_local.is_col_major());
    const bool is_int_power = fabs(power - int(power)) < 1e-4;
    const int n = ad_A.m();
    const char jobz = 'V';
    const char uplo = 'U';

    int lwork = -1, lrwork = -1, info = 0;
    std::complex<T>  *work;
    T *rwork;
    {
        work  = new std::complex<T>[1];
        rwork = new T[1];
        // query the optimal lwork and lrwork
        T *Wquery = new T[1];
        ScalapackConnector::pheev_f(jobz, uplo,
                n, A_local.c, 1, 1, ad_A.desc,
                Wquery, Z_local.c, 1, 1, ad_Z.desc, work, lwork, rwork, lrwork, info);
        lwork = int(work[0].real());
        lrwork = int(rwork[0]);
        delete [] work, Wquery, rwork;
    }
    // printf("lwork %d lrwork %d\n", lwork, lrwork); // debug

    work = new std::complex<T> [lwork];
    rwork = new T [lrwork];
    ScalapackConnector::pheev_f(jobz, uplo,
            n, A_local.c, 1, 1, ad_A.desc,
            W, Z_local.c, 1, 1, ad_Z.desc, work, lwork, rwork, lrwork, info);
    delete [] work, rwork;

    // check the number of non-singular eigenvalues,
    // using the fact that W is in ascending order
    n_filtered = n;
    for (int i = 0; i != n; i++)
        if (W[i] >= threshold)
        {
            n_filtered = i;
            break;
        }

    // filter and scale the eigenvalues, store in a temp array
    T W_temp[n];
    for (int i = 0; i != n_filtered; i++)
        W_temp[i] = 0.0;
    for (int i = n_filtered; i != n; i++)
    {
        if (W[i] < 0 && !is_int_power)
            printf("Warning! unfiltered negative eigenvalue with non-integer power: # %d ev = %f , pow = %f\n", i, W[i], power);
        if (fabs(W[i]) < 1e-10 && power < 0)
            printf("Warning! unfiltered nearly-singular eigenvalue with negative power: # %d ev = %f , pow = %f\n", i, W[i], power);
        W_temp[i] = std::pow(W[i], power);
    }

    // create scaled eigenvectors
    matrix_m<std::complex<T>> scaled_Z_local(Z_local);
    for (int i = 0; i != n; i++)
        ScalapackConnector::pscal_f(n, W_temp[i], scaled_Z_local.c, 1, 1+i, ad_Z.desc, 1);
    ScalapackConnector::pgemm_f('N', 'C', n, n, n, 1.0, Z_local.c, 1, 1, ad_Z.desc, scaled_Z_local.c, 1, 1, ad_Z.desc, 0.0, A_local.c, 1, 1, ad_A.desc);
    return scaled_Z_local;
}

template <typename T>
matrix_m<T> multiply_scalpack(const matrix_m<T> &m1_loc, const LIBRPA::Array_Desc &desc_m1,
                              const matrix_m<T> &m2_loc, const LIBRPA::Array_Desc &desc_m2,
                              const LIBRPA::Array_Desc &desc_prod)
{
    if (m1_loc.major() != m2_loc.major())
        throw std::invalid_argument("m1 and m2 in different major");
    if (desc_m1.n() != desc_m2.m())
        throw std::invalid_argument("m1.col != m2.row");
    const int m = desc_m1.m();
    const int n = desc_m2.n();
    const int k = desc_m1.n();
    matrix_m<T> prod_loc = init_local_mat<T>(desc_prod, m1_loc.major());
    if (m1_loc.is_col_major())
    {
        ScalapackConnector::pgemm_f('N', 'N', m, n, k, 1.0,
                m1_loc.c, 1, 1, desc_m1.desc,
                m2_loc.c, 1, 1, desc_m2.desc,
                0.0,
                prod_loc.c, 1, 1, desc_prod.desc);
    }
    else
    {
        throw std::logic_error("row-major parallel multiply is not implemented");
    }
    return prod_loc;
}

template <typename T>
void invert_scalapack(matrix_m<T> &m_loc, const LIBRPA::Array_Desc &desc_m)
{
    assert(m_loc.is_col_major());
    int info = 0;
    if (m_loc.is_col_major())
    {
        // NOTE: use local leading dimension desc_m.lld() will lead to invalid pointer error
        int *ipiv = new int [desc_m.m()];
        ScalapackConnector::pgetrf_f(desc_m.m(), desc_m.n(), m_loc.c, 1, 1, desc_m.desc, ipiv, info);
        // get optimized work size
        int lwork = -1, liwork = -1;
        {
            T *work = new T [1];
            int *iwork = new int [1];
            ScalapackConnector::pgetri_f(desc_m.m(), m_loc.c, 1, 1, desc_m.desc, ipiv, work, lwork, iwork, liwork, info);
            lwork = int(get_real(work[0]));
            liwork = iwork[0];
            // printf("lwork %d liwork %d\n", lwork, liwork);
            delete [] work;
            delete [] iwork;
        }

        T *work = new T[lwork];
        int *iwork = new int[liwork];
        ScalapackConnector::pgetri_f(desc_m.m(), m_loc.c, 1, 1, desc_m.desc, ipiv, work, lwork, iwork, liwork, info);
        delete [] iwork;
        delete [] work;
        delete [] ipiv;
        // printf("done free\n");
    }
    else
    {
        throw std::logic_error("row-major invert is not implemented");
    }
}
