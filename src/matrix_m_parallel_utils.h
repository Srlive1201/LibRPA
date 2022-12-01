#pragma once
#include "matrix_m.h"
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
