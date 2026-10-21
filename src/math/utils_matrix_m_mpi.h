#ifndef UTILS_MATRIX_M_MPI
#define UTILS_MATRIX_M_MPI

#pragma once
#include <omp.h>

#include <algorithm>
#include <cassert>
#include <functional>
#include <map>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "../gpu/la_connector.h"
#if defined(LIBRPA_USE_CUDA) || defined(LIBRPA_USE_HIP)
#include <ddla/ddla_connector.h>
#endif

#include "../mpi/base_blacs.h"
#include "../utils/profiler.h"
#include "matrix_m.h"
#include "scalapack_connector.h"

namespace librpa_int
{

template <typename T>
matrix_m<T> init_local_mat(const ArrayDesc &ad, MAJOR major)
{
    matrix_m<T> mat_lo(ad.m_loc(), ad.n_loc(), major);
    return mat_lo;
}

template <typename T>
matrix_m<T> get_local_mat(const matrix_m<T> &mat_go, const ArrayDesc &ad,
                          MAJOR major = MAJOR::AUTO)
{
    // assert the shape of global matrix conforms with the array descriptor
    assert(mat_go.nr() == ad.m() && mat_go.nc() == ad.n());

    if (major == MAJOR::AUTO) major = mat_go.major();
    const auto major_src = mat_go.major();
    matrix_m<T> mat_lo(ad.m_loc(), ad.n_loc(), major);
    const auto nr = mat_lo.nr();
    const auto nc = mat_lo.nc();

    if (major_src == MAJOR::ROW)
    {
        for (int i = 0; i != nr; i++)
        {
            const auto i_go = ad.indx_l2g_r(i);
            for (int j = 0; j != nc; j++)
            {
                const auto j_go = ad.indx_l2g_c(j);
                mat_lo(i, j) = mat_go(i_go, j_go);
            }
        }
    }
    else
    {
        for (int j = 0; j != nc; j++)
        {
            const auto j_go = ad.indx_l2g_c(j);
            for (int i = 0; i != nr; i++)
            {
                const auto i_go = ad.indx_l2g_r(i);
                mat_lo(i, j) = mat_go(i_go, j_go);
            }
        }
    }
    return mat_lo;
}

template <typename T>
matrix_m<T> get_local_mat(const T *pv, MAJOR major_pv, const ArrayDesc &ad, MAJOR major)
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

/*!
 * @brief Compute power of Hermitian matrix using BLACS
 *
 * @param  [in,out]  A_local       Process-local part of global matrix A to be powered
 * @param  [in]      ad_A          Array descriptor of A
 * @param  [out]     Z_local       Process-local part of the eigenvector matrix (Z) of A
 * @param  [in]      ad_Z          Array descriptor of Z
 * @param  [out]     n_filtered    Array descriptor of Z
 * @param  [out]     W             Eigenvalues of A, including those smaller than threshold
 * @param  [in]      power         Power to perform
 * @param  [in]      threshold     The threshold to filter the eigenvalues
 *
 * @retval           scale_Z       Eigenvectors scaled by the power of eigenvalues, using ad_Z
 */
template <typename T>
matrix_m<std::complex<T>> power_hemat_blacs(matrix_m<std::complex<T>> &A_local,
                                            const ArrayDesc &ad_A,
                                            matrix_m<std::complex<T>> &Z_local,
                                            const ArrayDesc &ad_Z,
                                            size_t &n_filtered, T *W, T power,
                                            const T &threshold = -1.e5)
{
    using global::ofs_myid;
    using global::profiler;

    profiler.start(__FUNCTION__, LIBRPA_VERBOSE_DEBUG);

    assert (A_local.is_col_major() && Z_local.is_col_major());
    const bool is_int_power = fabs(power - int(power)) < 1e-4;
    const size_t n = ad_A.m();
    const char jobz = 'V';
    const char uplo = 'U';

    // temporary array for heev with optmized block size
    const int blocksize_row_opt = std::min(ad_A.mb(), 128);
    const int blocksize_col_opt = std::min(ad_A.nb(), 128);

    // initialize descriptor of array A for optimized block size
    ArrayDesc ad_A_opt(ad_A.ictxt());
    ad_A_opt.init(n, n, blocksize_row_opt, blocksize_col_opt, 0, 0);
    auto A_local_opt = init_local_mat<std::complex<T>>(ad_A_opt, MAJOR::COL);
    // NOTE: imply A and Z should be in the same context
    ScalapackConnector::pgemr2d_f(n, n, A_local.ptr(), 1, 1, ad_A.desc,
                                  A_local_opt.ptr(), 1, 1, ad_A_opt.desc, ad_A.ictxt());

    // initialize descriptor of array Z for optimized block size
    if (ad_A.ictxt() != ad_Z.ictxt())
    {
        ofs_myid << "Warning(power_hemat_blacs): input contexts of A and Z are different!" << std::endl;
    }
    ArrayDesc ad_Z_opt(ad_Z.ictxt());
    ad_Z_opt.init(n, n, blocksize_row_opt, blocksize_col_opt, 0, 0);
    auto Z_local_opt = init_local_mat<std::complex<T>>(ad_Z_opt, MAJOR::COL);
    // printf("Z_local_opt size: %d\n", Z_local_opt.size());

    if (n == 1)
    {
        // The scalar eigensystem is exact. Some PZHEEV implementations return
        // (1,1) instead of the normalized eigenvector (1,0) in this quick case.
        W[0] = ad_A_opt.is_src() ? A_local_opt(0, 0).real() : T(0);
        if constexpr (std::is_same<T, double>::value)
            Cdgsum2d(ad_A_opt.ictxt(), "A", " ", 1, 1, W, 1, -1, -1);
        else
            Csgsum2d(ad_A_opt.ictxt(), "A", " ", 1, 1, W, 1, -1, -1);
        if (ad_Z_opt.is_src()) Z_local_opt(0, 0) = std::complex<T>(1, 0);
    }
    else
    {
        profiler.start("power_hemat_blacs_1", LIBRPA_VERBOSE_DEBUG);
        int lwork = -1, lrwork = -1, info = 0;
        std::complex<T> *work;
        T *rwork;
        {
            work = new std::complex<T>[1];
            rwork = new T[1];
            // query the optimal lwork and lrwork
            T *Wquery = new T[1];
            // librpa_int::global::lib_printf("power_hemat_blacs descA %s\n",
            // ad_A.info_desc().c_str());
            ScalapackConnector::pheev_f(jobz, uplo, n, A_local_opt.ptr(), 1, 1, ad_A_opt.desc,
                                        Wquery, Z_local_opt.ptr(), 1, 1, ad_A_opt.desc, work, lwork,
                                        rwork, lrwork, info);
            lwork = std::max(int(work[0].real()), 1);
            lrwork = std::max(int(rwork[0]), 1);
            delete[] work;
            delete[] Wquery;
            delete[] rwork;
        }
        profiler.stop("power_hemat_blacs_1");

        profiler.start("power_hemat_blacs_2", LIBRPA_VERBOSE_DEBUG);
        work = new std::complex<T>[lwork];
        rwork = new T[lrwork];
        ScalapackConnector::pheev_f(jobz, uplo, n, A_local_opt.ptr(), 1, 1, ad_A_opt.desc, W,
                                    Z_local_opt.ptr(), 1, 1, ad_Z_opt.desc, work, lwork, rwork,
                                    lrwork, info);
        delete[] work;
        delete[] rwork;
        profiler.stop("power_hemat_blacs_2");
    }
    // Optimized A no longer used
    profiler.start("power_hemat_blacs_3", LIBRPA_VERBOSE_DEBUG);
    A_local_opt.clear();
    // send back the eigenvector matrix
    ScalapackConnector::pgemr2d_f(n, n, Z_local_opt.ptr(), 1, 1, ad_Z_opt.desc, Z_local.ptr(), 1, 1,
                                  ad_Z.desc, ad_Z.ictxt());
    profiler.stop("power_hemat_blacs_3");

    // check the number of non-singular eigenvalues,
    // using the fact that W is in ascending order
    n_filtered = n;
    for (size_t i = 0; i != n; i++)
        if (W[i] >= threshold)
        {
            n_filtered = i;
            break;
        }

    // filter and scale the eigenvalues, store in a temp array
    profiler.start("power_hemat_blacs_4", LIBRPA_VERBOSE_DEBUG);
    std::vector<T> W_temp(n);
    for (size_t i = 0; i < n_filtered; i++) W_temp[i] = 0.0;
    for (size_t i = n_filtered; i != n; i++)
    {
        if (W[i] < 0 && !is_int_power)
        {
            global::lib_printf(
                LIBRPA_VERBOSE_WARN,
                "Warning! unfiltered negative eigenvalue with non-integer power: # %d ev = %f , "
                "pow = %f\n",
                i, W[i], power);
        }
        if (fabs(W[i]) < 1e-10 && power < 0)
        {
            global::lib_printf(
                LIBRPA_VERBOSE_WARN,
                "Warning! unfiltered nearly-singular eigenvalue with negative power: # %d ev = %f "
                ", pow = %f\n",
                i, W[i], power);
        }
        W_temp[i] = std::pow(W[i], power);
    }
    profiler.stop("power_hemat_blacs_4");
    // debug print
    // for (int i = 0; i != n; i++)
    // {
    //     librpa_int::global::lib_printf("%d %f %f\n", i, W[i], W_temp[i]);
    // }

    profiler.start("power_hemat_blacs_5", LIBRPA_VERBOSE_DEBUG);
    // create scaled eigenvectors
    auto scaled_opt = Z_local_opt.copy();
    for (size_t i = 0; i != n; i++)
    {
        ScalapackConnector::pscal_f(n, W_temp[i], scaled_opt.ptr(), 1, 1+i, ad_Z_opt.desc, 1);
    }
    ScalapackConnector::pgemm_f('N', 'C', n, n, n, 1.0, Z_local_opt.ptr(), 1, 1, ad_Z_opt.desc,
                                scaled_opt.ptr(), 1, 1, ad_Z_opt.desc, 0.0, A_local.ptr(), 1, 1,
                                ad_A.desc);
    auto scaled = Z_local.copy();
    // send back the scaled eigenvector matrix with descriptor using optimized block size to that
    // with input descriptor
    ScalapackConnector::pgemr2d_f(n, n, scaled_opt.ptr(), 1, 1, ad_Z_opt.desc, scaled.ptr(), 1, 1,
                                  ad_Z.desc, ad_Z.ictxt());
    profiler.stop("power_hemat_blacs_5");
    profiler.stop(__FUNCTION__);

    return scaled;
}

template <typename T>
matrix_m<std::complex<T>> power_hemat_blacs_desc(matrix_m<std::complex<T>> &A_local,
                                                 const ArrayDesc &ad_A,
                                                 matrix_m<std::complex<T>> &Z_local,
                                                 const ArrayDesc &ad_Z, size_t &n_filtered, T *W,
                                                 T power, const T &threshold = -1.e5)
{
    using global::ofs_myid;
    using global::profiler;

    profiler.start(__FUNCTION__, LIBRPA_VERBOSE_DEBUG);

    A_local *= -1.0;
    assert(A_local.is_col_major() && Z_local.is_col_major());
    const bool is_int_power = fabs(power - int(power)) < 1e-4;
    const int n = ad_A.m();
    const char jobz = 'V';
    const char uplo = 'U';

    // temporary array for heev with optmized block size
    const int blocksize_row_opt = std::min(ad_A.mb(), 128);
    const int blocksize_col_opt = std::min(ad_A.nb(), 128);

    // initialize descriptor of array A for optimized block size
    ArrayDesc ad_A_opt(ad_A.ictxt());
    ad_A_opt.init(n, n, blocksize_row_opt, blocksize_col_opt, 0, 0);
    auto A_local_opt = init_local_mat<std::complex<T>>(ad_A_opt, MAJOR::COL);
    // NOTE: imply A and Z should be in the same context
    ScalapackConnector::pgemr2d_f(n, n, A_local.ptr(), 1, 1, ad_A.desc, A_local_opt.ptr(), 1, 1,
                                  ad_A_opt.desc, ad_A.ictxt());

    // initialize descriptor of array Z for optimized block size
    if (ad_A.ictxt() != ad_Z.ictxt())
    {
        ofs_myid << "Warning(power_hemat_blacs): input contexts of A and Z are different!\n";
    }
    ArrayDesc ad_Z_opt(ad_Z.ictxt());
    ad_Z_opt.init(n, n, blocksize_row_opt, blocksize_col_opt, 0, 0);
    auto Z_local_opt = init_local_mat<std::complex<T>>(ad_Z_opt, MAJOR::COL);
    // printf("Z_local_opt size: %d\n", Z_local_opt.size());

    profiler.start("power_hemat_blacs_1", LIBRPA_VERBOSE_DEBUG);
    int lwork = -1, lrwork = -1, info = 0;
    std::complex<T> *work;
    T *rwork;
    {
        work = new std::complex<T>[1];
        rwork = new T[1];
        // query the optimal lwork and lrwork
        T *Wquery = new T[1];
        // LIBRPA::utils::lib_printf("power_hemat_blacs descA %s\n", ad_A.info_desc().c_str());
        ScalapackConnector::pheev_f(jobz, uplo, n, A_local_opt.ptr(), 1, 1, ad_A_opt.desc, Wquery,
                                    Z_local_opt.ptr(), 1, 1, ad_A_opt.desc, work, lwork, rwork,
                                    lrwork, info);
        lwork = std::max(int(work[0].real()), 1);
        lrwork = std::max(int(rwork[0]), 1);
        delete[] work;
        delete[] Wquery;
        delete[] rwork;
    }
    profiler.stop("power_hemat_blacs_1");

    profiler.start("power_hemat_blacs_2", LIBRPA_VERBOSE_DEBUG);
    work = new std::complex<T>[lwork];
    rwork = new T[lrwork];
    ScalapackConnector::pheev_f(jobz, uplo, n, A_local_opt.ptr(), 1, 1, ad_A_opt.desc, W,
                                Z_local_opt.ptr(), 1, 1, ad_Z_opt.desc, work, lwork, rwork, lrwork,
                                info);
    delete[] work;
    delete[] rwork;
    for (int i = 0; i != n; i++)
    {
        W[i] *= -1.0;
    }
    // Optimized A no longer used
    A_local_opt.clear();

    // change to descending order
    /*auto Z_tmp = Z_local_opt.copy();
    for (int i = 0; i != n; i++)
    {
        for (int j = 0; j != n; j++)
            Z_local_opt.dataobj[i * n + j] = Z_tmp.dataobj[(n - i - 1) * n + j];
    }*/

    ScalapackConnector::pgemr2d_f(n, n, Z_local_opt.ptr(), 1, 1, ad_Z_opt.desc, Z_local.ptr(), 1, 1,
                                  ad_Z.desc, ad_Z.ictxt());
    profiler.stop("power_hemat_blacs_2");

    // check the number of non-singular eigenvalues,
    // using the fact that W is in descending order
    n_filtered = as_size(n);
    for (int i = 0; i != n; i++)
        if (W[n - 1 - i] >= threshold)
        {
            n_filtered = as_size(i);
            break;
        }

    // filter and scale the eigenvalues, store in a temp array
    profiler.start("power_hemat_blacs_3", LIBRPA_VERBOSE_DEBUG);
    std::vector<T> W_temp(as_size(n), T(0));
    const auto n_kept = as_size(n) - n_filtered;
    for (std::size_t i = 0; i != n_kept; i++)
    {
        if (W[i] < 0 && !is_int_power)
        {
            global::lib_printf(
                LIBRPA_VERBOSE_WARN,
                "Warning! unfiltered negative eigenvalue with non-integer power: # %zu ev = %f , "
                "pow = %f\n",
                i, W[i], power);
        }
        if (fabs(W[i]) < 1e-10 && power < 0)
        {
            global::lib_printf(
                LIBRPA_VERBOSE_WARN,
                "Warning! unfiltered nearly-singular eigenvalue with negative power: # %zu ev = %f "
                ", pow = %f\n",
                i, W[i], power);
        }
        W_temp[i] = std::pow(W[i], power);
    }
    profiler.stop("power_hemat_blacs_3");
    // debug print
    // for (int i = 0; i != n; i++)
    // {
    //     LIBRPA::utils::lib_printf("%d %f %f\n", i, W[i], W_temp[i]);
    // }

    profiler.start("power_hemat_blacs_4", LIBRPA_VERBOSE_DEBUG);

    // create scaled eigenvectors
    auto scaled_opt = Z_local_opt.copy();
    for (int i = 0; i != n; i++)
    {
        ScalapackConnector::pscal_f(n, W_temp[i], scaled_opt.ptr(), 1, 1 + i, ad_Z_opt.desc, 1);
    }
    ScalapackConnector::pgemm_f('N', 'C', n, n, n, 1.0, Z_local_opt.ptr(), 1, 1, ad_Z_opt.desc,
                                scaled_opt.ptr(), 1, 1, ad_Z_opt.desc, 0.0, A_local.ptr(), 1, 1,
                                ad_A.desc);
    auto scaled = Z_local.copy();
    // send back the scaled eigenvector matrix with descriptor using optimized block size to that
    // with input descriptor
    ScalapackConnector::pgemr2d_f(n, n, scaled_opt.ptr(), 1, 1, ad_Z_opt.desc, scaled.ptr(), 1, 1,
                                  ad_Z.desc, ad_Z.ictxt());
    profiler.stop("power_hemat_blacs_4");
    profiler.stop(__FUNCTION__);

    return scaled;
}

template <typename T>
matrix_m<std::complex<T>> power_hemat_blacs_real(matrix_m<std::complex<T>> &A_local,
                                                 const ArrayDesc &ad_A,
                                                 matrix_m<std::complex<T>> &Z_local,
                                                 const ArrayDesc &ad_Z, size_t &n_filtered,
                                                 T *W, T power, const T &threshold = -1.e5)
{
    using global::ofs_myid;
    using global::profiler;

    profiler.start(__FUNCTION__, LIBRPA_VERBOSE_DEBUG);
    // Step 1: Extract real part of the complex matrix
    auto A_local_real = A_local.get_real();
    A_local_real *= -1.0;
    assert(A_local.is_col_major() && Z_local.is_col_major());
    const bool is_int_power = fabs(power - int(power)) < 1e-4;
    const int n = ad_A.m();
    const char jobz = 'V';
    const char uplo = 'U';

    // temporary array for syev with optimized block size
    const int blocksize_row_opt = std::min(ad_A.mb(), 128);
    const int blocksize_col_opt = std::min(ad_A.nb(), 128);

    // initialize descriptor of array A for optimized block size
    ArrayDesc ad_A_opt(ad_A.ictxt());
    ad_A_opt.init(n, n, blocksize_row_opt, blocksize_col_opt, 0, 0);
    // Initialize as real matrix for diagonalization
    auto A_local_opt = init_local_mat<T>(ad_A_opt, MAJOR::COL);
    ScalapackConnector::pgemr2d_f(n, n, A_local_real.ptr(), 1, 1, ad_A.desc, A_local_opt.ptr(), 1,
                                  1, ad_A_opt.desc, ad_A.ictxt());

    // initialize descriptor of array Z for optimized block size
    if (ad_A.ictxt() != ad_Z.ictxt())
    {
        ofs_myid << "Warning(power_hemat_blacs): input contexts of A and Z are different!\n";
    }
    ArrayDesc ad_Z_opt(ad_Z.ictxt());
    ad_Z_opt.init(n, n, blocksize_row_opt, blocksize_col_opt, 0, 0);
    // Initialize Z as real matrix initially
    auto Z_local_opt = init_local_mat<T>(ad_Z_opt, MAJOR::COL);

    profiler.start("power_hemat_blacs_1", LIBRPA_VERBOSE_DEBUG);
    int lwork = -1, lrwork = -1, info = 0;

    // Query optimal workspace using the provided psyev_f interface
    T work_query, rwork_query, Wquery;
    ScalapackConnector::psyev_f(jobz, uplo, n, A_local_opt.ptr(), 1, 1, ad_A_opt.desc, &Wquery,
                                Z_local_opt.ptr(), 1, 1, ad_A_opt.desc, &work_query, lwork, &rwork_query,
                                lrwork, info);
    lwork = std::max(as_int(work_query), 1);
    lrwork = std::max(as_int(rwork_query), 1);
    profiler.stop("power_hemat_blacs_1");

    std::vector<T> work(lwork);
    std::vector<T> rwork(lrwork);
    profiler.start("power_hemat_blacs_2", LIBRPA_VERBOSE_DEBUG);
    // Perform real symmetric diagonalization using the provided interface
    ofs_myid << lwork << " " << lrwork << " " << n << std::endl;
    ScalapackConnector::psyev_f(jobz, uplo, n, A_local_opt.ptr(), 1, 1, ad_A_opt.desc, W,
                                Z_local_opt.ptr(), 1, 1, ad_Z_opt.desc, work.data(), lwork, rwork.data(), lrwork,
                                info);
    work.clear();
    rwork.clear();
    for (int i = 0; i < n; i++)
    {
        W[i] *= -1.0;
    }
    A_local_opt.clear();

    // Convert real eigenvectors to complex (with zero imaginary part)
    auto Z_local_opt_complex = Z_local_opt.to_complex();
    // Transfer using complex matrix descriptor
    ScalapackConnector::pgemr2d_f(n, n, Z_local_opt_complex.ptr(), 1, 1, ad_Z_opt.desc,
                                  Z_local.ptr(), 1, 1, ad_Z.desc, ad_Z.ictxt());
    profiler.stop("power_hemat_blacs_2");

    // Check number of non-singular eigenvalues
    n_filtered = as_size(n);
    for (int i = 0; i < n; i++)
    {
        if (W[n - 1 - i] >= threshold)
        {
            n_filtered = as_size(i);
            break;
        }
    }

    // Filter and scale eigenvalues
    profiler.start("power_hemat_blacs_3", LIBRPA_VERBOSE_DEBUG);
    std::vector<T> W_temp(as_size(n), T(0));
    const std::size_t n_kept = as_size(n) - n_filtered;
    for (std::size_t i = 0; i < n_kept; i++)
    {
        if (W[i] < 0 && !is_int_power)
        {
            global::lib_printf(
                LIBRPA_VERBOSE_WARN,
                "Warning! unfiltered negative eigenvalue with non-integer power: # %zu ev = %f , "
                "pow = %f\n",
                i, W[i], power);
        }
        if (fabs(W[i]) < 1e-10 && power < 0)
        {
            global::lib_printf(
                LIBRPA_VERBOSE_WARN,
                "Warning! unfiltered nearly-singular eigenvalue with negative power: # %zu ev = %f "
                ", pow = %f\n",
                i, W[i], power);
        }
        W_temp[i] = std::pow(W[i], power);
    }
    profiler.stop("power_hemat_blacs_3");

    profiler.start("power_hemat_blacs_4", LIBRPA_VERBOSE_DEBUG);
    // Scale eigenvectors using complex matrix operations
    auto scaled_opt = Z_local_opt_complex.copy();
    for (int i = 0; i < n; i++)
    {
        // Use ScaLAPACK scaling function with complex matrix
        ScalapackConnector::pscal_f(n, W_temp[i], scaled_opt.ptr(), 1, 1 + i, ad_Z_opt.desc, 1);
    }

    // Compute Z * diag(W_temp) * Z^H
    ScalapackConnector::pgemm_f('N', 'C', n, n, n, 1.0, Z_local_opt_complex.ptr(), 1, 1,
                                ad_Z_opt.desc, scaled_opt.ptr(), 1, 1, ad_Z_opt.desc, 0.0,
                                A_local.ptr(), 1, 1, ad_A.desc);

    auto scaled = Z_local.copy();
    // Transfer back to original complex matrix
    ScalapackConnector::pgemr2d_f(n, n, scaled_opt.ptr(), 1, 1, ad_Z_opt.desc, scaled.ptr(), 1, 1,
                                  ad_Z.desc, ad_Z.ictxt());
    profiler.stop("power_hemat_blacs_4");
    profiler.stop(__FUNCTION__);

    return scaled;
}

template <typename T, typename Treal = typename to_real<T>::type>
void print_matrix_mm_parallel(std::ostream &os, const matrix_m<T> &mat_loc, const ArrayDesc &ad,
                              const std::string &comment = "", Treal threshold = 1e-15, bool row_first = true)
{
    ArrayDesc ad_fb(ad.ictxt());
    const int nr = ad.m(), nc = ad.n();
    const int irsrc = ad.irsrc(), icsrc = ad.icsrc();
    ad_fb.init(nr, nc, nr, nc, irsrc, icsrc);
    matrix_m<T> mat_glo = init_local_mat<T>(ad_fb, mat_loc.major());

    ScalapackConnector::pgemr2d_f(nr, nc, mat_loc.ptr(), 1, 1, ad.desc, mat_glo.ptr(), 1, 1,
                                  ad_fb.desc, ad.ictxt());
    if (ad_fb.is_src() && os.good())
    {
        print_matrix_mm(mat_glo, os, comment, threshold, row_first);
    }
}

template <typename T, typename Treal = typename to_real<T>::type>
void print_matrix_mm_file_parallel(const string &fn, const matrix_m<T> &mat_loc,
                                   const ArrayDesc &ad, const std::string &comment = "",
                                   Treal threshold = 1e-15, bool row_first = true)
{
    ofstream fs;
    if (ad.is_src())
        fs.open(fn);
    print_matrix_mm_parallel(fs, mat_loc, ad, comment, threshold, row_first);
    if (ad.is_src())
        fs.close();
}

template <typename T, typename Treal = typename to_real<T>::type>
void write_matrix_elsi_csc_parallel(const string &fn, const matrix_m<T> &mat_loc, const ArrayDesc &ad, Treal threshold = 1e-15)
{
    ArrayDesc ad_fb(ad.ictxt());
    const int nr = ad.m(), nc = ad.n();
    const int irsrc = ad.irsrc(), icsrc = ad.icsrc();
    ad_fb.init(nr, nc, nr, nc, irsrc, icsrc);
    matrix_m<T> mat_glo = init_local_mat<T>(ad_fb, mat_loc.major());

    ScalapackConnector::pgemr2d_f(nr, nc, mat_loc.ptr(), 1, 1, ad.desc, mat_glo.ptr(), 1, 1, ad_fb.desc, ad.ictxt());
    if (ad_fb.is_src()) write_matrix_elsi_csc(mat_glo, fn, threshold);
    ad.barrier();
}

template <typename T>
matrix_m<T> multiply_scalapack(const matrix_m<T> &m1_loc, const ArrayDesc &desc_m1,
                               const matrix_m<T> &m2_loc, const ArrayDesc &desc_m2,
                               const ArrayDesc &desc_prod)
{
    if (m1_loc.major() != m2_loc.major())
    {
        throw std::invalid_argument("m1 and m2 in different major");
    }

    const int m = desc_m1.m();
    const int n = desc_m2.n();
    const int k = desc_m1.n();

    const int m2 = desc_m2.m();
    if (k != m2)
    {
        throw std::invalid_argument("m1.col != m2.row: " + std::to_string(k) + " != " + std::to_string(m2));
    }

    matrix_m<T> prod_loc = init_local_mat<T>(desc_prod, m1_loc.major());
    if (m1_loc.is_col_major())
    {
        ScalapackConnector::pgemm_f('N', 'N', m, n, k, 1.0,
                m1_loc.ptr(), 1, 1, desc_m1.desc,
                m2_loc.ptr(), 1, 1, desc_m2.desc,
                0.0,
                prod_loc.ptr(), 1, 1, desc_prod.desc);
    }
    else
    {
        throw std::logic_error("row-major parallel multiply is not implemented");
    }
    return prod_loc;
}

template <typename T>
void invert_scalapack(matrix_m<T> &m_loc, const ArrayDesc &desc_m)
{
    assert(m_loc.is_col_major());
    int info = 0;
    if (m_loc.is_col_major())
    {
        // NOTE: use local leading dimension desc_m.lld() will lead to invalid pointer error
        int *ipiv = new int [desc_m.m()];
        // librpa_int::global::lib_printf("invert_scalpack desc %s\n", desc_m.info_desc().c_str());
        ScalapackConnector::pgetrf_f(desc_m.m(), desc_m.n(), m_loc.ptr(), 1, 1, desc_m.desc, ipiv,
                                     info);
        // get optimized work size
        int lwork = -1, liwork = -1;
        {
            T *work = new T [1];
            int *iwork = new int [1];
            ScalapackConnector::pgetri_f(desc_m.m(), m_loc.ptr(), 1, 1, desc_m.desc, ipiv, work,
                                         lwork, iwork, liwork, info);
            lwork = std::max(as_int(get_real(work[0])), 1);
            liwork = std::max(iwork[0], 1);
            // librpa_int::global::lib_printf("lwork %d liwork %d\n", lwork, liwork);
            delete [] work;
            delete [] iwork;
        }

        T *work = new T[lwork];
        int *iwork = new int[liwork];
        ScalapackConnector::pgetri_f(desc_m.m(), m_loc.ptr(), 1, 1, desc_m.desc, ipiv, work, lwork,
                                     iwork, liwork, info);
        delete [] iwork;
        delete [] work;
        delete [] ipiv;
        // librpa_int::global::lib_printf("done free\n");
    }
    else
    {
        throw std::logic_error("row-major invert is not implemented");
    }
}

} /* end of namespace librpa_int */


#endif // UTILS_MATRIX_M_MPI_H
