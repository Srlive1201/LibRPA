#pragma once
#include <omp.h>

#include <functional>
#include <map>

#include "atomic_basis.h"
#include "base_blacs.h"
#include "matrix_m.h"
#include "profiler.h"
#include "scalapack_connector.h"
#include "constants.h"
#include "utils_io.h"
#ifdef LIBRPA_USE_LIBRI
#include <RI/global/Tensor.h>
#else
#include "libri_stub.h"
#endif

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
    // cout << str(mat_lo);
    for (int iI = 0; iI != row_nb; iI++)
    {
        int ilo = ad.indx_g2l_r(row_start_id + iI);
        // LIBRPA::utils::lib_printf("row igo %d ilo %d\n", row_start_id + iI, ilo);
        if (ilo < 0) continue;
        for (int jJ = 0; jJ != col_nb; jJ++)
        {
            int jlo = ad.indx_g2l_c(col_start_id + jJ);
            // LIBRPA::utils::lib_printf("col jgo %d jlo %d\n", col_start_id + jJ, jlo);
            if (jlo < 0) continue;
            // cout << "pickerID " << picker(row_nb, iI, col_nb, jJ) << " " << pvIJ[picker(row_nb, iI, col_nb, jJ)] << " matloc befor " << mat_lo(ilo, jlo);
            mat_lo(ilo, jlo) += alpha * pvIJ[picker(row_nb, iI, col_nb, jJ)];
            // cout << " after " << mat_lo(ilo, jlo) << endl;
        }
    }
}

template <typename Tdst, typename Tsrc, typename TA, typename TC, typename TAC = std::pair<TA, TC>>
void collect_block_from_IJ_storage_tensor(
    matrix_m<Tdst> &mat_lo,
    const LIBRPA::Array_Desc &ad,
    const LIBRPA::AtomicBasis &atbasis_row,
    const LIBRPA::AtomicBasis &atbasis_col,
    const TC &cell,
    Tdst alpha, const std::map<TA,std::map<TAC, RI::Tensor<Tsrc>>> &TMAP)
{
    // assert(mat_lo.nr() == ad.m_loc() && mat_lo.nc() == ad.n_loc());
    assert(ad.m() == atbasis_row.nb_total && ad.n() == atbasis_col.nb_total );

    matrix_m<Tdst> tmp_loc(mat_lo.nr(),mat_lo.nc(), MAJOR::ROW);
    size_t cp_size= ad.n_loc()*sizeof(Tdst);

    omp_lock_t mat_lock;
    omp_init_lock(&mat_lock);

    #pragma omp parallel for
    for (int ilo = 0; ilo != ad.m_loc(); ilo++)
    {
        int I_loc, J_loc, i_ab, j_ab;
        int i_gl = ad.indx_l2g_r(ilo);
        atbasis_row.get_local_index(i_gl, I_loc, i_ab);
        vector<Tdst> tmp_loc_row(ad.n_loc());
        for (int jlo = 0; jlo != ad.n_loc(); jlo++)
        {
            int j_gl = ad.indx_l2g_c(jlo);
            atbasis_col.get_local_index(j_gl, J_loc, j_ab);
            //Tdst temp;
            // LIBRPA::utils::lib_printf("i_gl I_loc i_ab %d %d %d j_gl J_loc j_ab %d %d %d\n", i_gl, I_loc, i_ab, j_gl, J_loc, j_ab);
            tmp_loc_row[jlo] = TMAP.at(I_loc).at({J_loc, cell})(i_ab, j_ab);
        }
        Tdst *row_ptr=tmp_loc.ptr() + ilo* ad.n_loc();
        omp_set_lock(&mat_lock);
        memcpy(row_ptr,&tmp_loc_row[0],cp_size);
        omp_unset_lock(&mat_lock);
    }
    #pragma omp barrier
    omp_destroy_lock(&mat_lock);
    if(mat_lo.is_col_major())
    {
        tmp_loc.swap_to_col_major();
    }
    mat_lo = tmp_loc;
}

template <typename Tdst, typename Tsrc, typename TA, typename TAC>
void collect_block_from_IJ_storage_tensor_transform(
    matrix_m<Tdst> &mat_lo,
    const LIBRPA::Array_Desc &ad,
    const LIBRPA::AtomicBasis &atbasis_row,
    const LIBRPA::AtomicBasis &atbasis_col,
    const std::function<Tdst(const TA &, const TAC &)> &transform,
    const std::map<TA,std::map<TAC, RI::Tensor<Tsrc>>> &TMAP)
{
    // assert(mat_lo.nr() == ad.m_loc() && mat_lo.nc() == ad.n_loc());
    assert(ad.m() == atbasis_row.nb_total && ad.n() == atbasis_col.nb_total );

    matrix_m<Tdst> tmp_loc(mat_lo.nr(),mat_lo.nc(), MAJOR::ROW);
    size_t cp_size= ad.n_loc()*sizeof(Tdst);

    omp_lock_t mat_lock;
    omp_init_lock(&mat_lock);

    #pragma omp parallel for
    for (int ilo = 0; ilo != ad.m_loc(); ilo++)
    {
        int I_loc, J_loc, i_ab, j_ab;
        int i_gl = ad.indx_l2g_r(ilo);
        atbasis_row.get_local_index(i_gl, I_loc, i_ab);
        vector<Tdst> tmp_loc_row(ad.n_loc(), 0);
        for (int jlo = 0; jlo != ad.n_loc(); jlo++)
        {
            int j_gl = ad.indx_l2g_c(jlo);
            atbasis_col.get_local_index(j_gl, J_loc, j_ab);
            if (TMAP.count(I_loc) > 0)
            {
                for (const auto &acell_mat: TMAP.at(I_loc))
                {
                    const auto &acell = acell_mat.first;
                    const auto &mat = acell_mat.second;
                    if (acell.first != J_loc)
                    {
                        continue;
                    }
                    tmp_loc_row[jlo] += mat(i_ab, j_ab) * transform(I_loc, acell);
                }
            }
        }
        Tdst *row_ptr=tmp_loc.ptr() + ilo* ad.n_loc();
        omp_set_lock(&mat_lock);
        memcpy(row_ptr,&tmp_loc_row[0],cp_size);
        omp_unset_lock(&mat_lock);
    }
    #pragma omp barrier
    omp_destroy_lock(&mat_lock);
    if(mat_lo.is_col_major())
    {
        tmp_loc.swap_to_col_major();
    }
    mat_lo = tmp_loc;
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
    const std::function<Tdst(Tdst)> filter = conjugate ? [](Tdst a) { return get_conj(a); }
                                                 : [](Tdst a) { return a; };
    const int I_nb = atbasis.get_atom_nb(I);
    const int J_nb = atbasis.get_atom_nb(J);
    const auto picker = Indx_pickers_2d[major_pv];
    int I_loc, J_loc, i_ab, j_ab;
    for (int ilo = 0; ilo != ad.m_loc(); ilo++)
        for (int jlo = 0; jlo != ad.n_loc(); jlo++)
        {
            int i_gl = ad.indx_l2g_r(ilo);
            int j_gl = ad.indx_l2g_c(jlo);
            atbasis.get_local_index(i_gl, I_loc, i_ab);
            atbasis.get_local_index(j_gl, J_loc, j_ab);
            Tdst temp;
            // LIBRPA::utils::lib_printf("i_gl I_loc i_ab %d %d %d j_gl J_loc j_ab %d %d %d\n", i_gl, I_loc, i_ab, j_gl, J_loc, j_ab);
            if ((I_loc == I) && (J_loc == J))
            {
                // in the same block
                temp = pvIJ[picker(I_nb, i_ab, J_nb, j_ab)];
                // cout << "same block pickerID " << picker(I_nb, i_ab, J_nb, j_ab) << " " << temp << " matloc befor " << mat_lo(ilo, jlo);
            }
            else if ((I_loc == J) && (J_loc == I))
            {
                // in the opposite block
                temp = filter(pvIJ[picker(I_nb, j_ab, J_nb, i_ab)]);
                // cout << "opposite block pickerID " << picker(I_nb, j_ab, J_nb, i_ab) << " " << temp << " matloc befor " << mat_lo(ilo, jlo);
            }
            else
                continue;
            mat_lo(ilo, jlo) += alpha * temp;
            // cout << " after " << mat_lo(ilo, jlo) << endl;
        }
}

template <typename Tdst, typename TA, typename TC, typename TAC = std::pair<TA, TC>>
void collect_block_from_ALL_IJ_Tensor(
    matrix_m<Tdst> &mat_lo,
    const LIBRPA::Array_Desc &ad,
    const LIBRPA::AtomicBasis &atbasis,
    const TC &cell,
    bool conjugate,
    Tdst alpha, const std::map<TA,std::map<TAC, RI::Tensor<Tdst>>> &TMAP, MAJOR major_pv)
{
    // assert(mat_lo.nr() == ad.m_loc() && mat_lo.nc() == ad.n_loc());
    assert(ad.m() == atbasis.nb_total && ad.n() == atbasis.nb_total);
    matrix_m<Tdst> tmp_loc(mat_lo.nr(),mat_lo.nc(), MAJOR::ROW);
    size_t cp_size= ad.n_loc()*sizeof(Tdst);

#ifdef LIBRPA_DEBUG
    LIBRPA::envs::ofs_myid << "cp_size: " << cp_size << endl;
    LIBRPA::envs::ofs_myid << "Available TMAP keys: ";
    print_keys(LIBRPA::envs::ofs_myid, TMAP);
    LIBRPA::envs::ofs_myid << endl;
    LIBRPA::envs::ofs_myid << "mat_lo dims: " << mat_lo.nr() << " " << mat_lo.nc() << endl;
    LIBRPA::envs::ofs_myid << "ad_loc dims: " << ad.m_loc() << " " << ad.n_loc() << endl;
#endif

    omp_lock_t mat_lock;
    omp_init_lock(&mat_lock);

    #pragma omp parallel for
    for (int ilo = 0; ilo != ad.m_loc(); ilo++)
    {
        int I_loc, J_loc, i_ab, j_ab;
        int i_gl = ad.indx_l2g_r(ilo);
        atbasis.get_local_index(i_gl, I_loc, i_ab);
        vector<Tdst> tmp_loc_row(ad.n_loc());
        for (int jlo = 0; jlo != ad.n_loc(); jlo++)
        {

            int j_gl = ad.indx_l2g_c(jlo);

            atbasis.get_local_index(j_gl, J_loc, j_ab);
            //Tdst temp;
            // LIBRPA::utils::lib_printf("i_gl I_loc i_ab %d %d %d j_gl J_loc j_ab %d %d %d\n", i_gl, I_loc, i_ab, j_gl, J_loc, j_ab);
            if(I_loc<=J_loc)
            {
                tmp_loc_row[jlo]= TMAP.at(I_loc).at({J_loc, cell})(i_ab,j_ab);
            }
            else
            {
                Tdst tmp_ele;
                tmp_ele= TMAP.at(J_loc).at({I_loc, cell})(j_ab,i_ab);
                if(conjugate)
                    tmp_ele=get_conj(tmp_ele);
                tmp_loc_row[jlo]=tmp_ele;
            }
        }
        Tdst *row_ptr=tmp_loc.ptr()+ilo* ad.n_loc();
        omp_set_lock(&mat_lock);
        memcpy(row_ptr,&tmp_loc_row[0],cp_size);
        omp_unset_lock(&mat_lock);
    }
    #pragma omp barrier
    omp_destroy_lock(&mat_lock);
    if(mat_lo.is_col_major())
    {
        tmp_loc.swap_to_col_major();
    }
    mat_lo=tmp_loc;
}


//! collect a IJ-pair storage from a 2D block local matrix
template <typename T>
void map_block_to_IJ_storage(map<int, map<int, matrix_m<T>>> &IJmap,
                             const LIBRPA::AtomicBasis &atbasis_row,
                             const LIBRPA::AtomicBasis &atbasis_col,
                             const matrix_m<T> &mat_lo,
                             const LIBRPA::Array_Desc &desc, MAJOR major_map)
{
    assert(desc.m() == atbasis_row.nb_total && desc.n() == atbasis_col.nb_total);
    int I, J, iI, jJ;
    for (int i_lo = 0; i_lo != desc.m_loc(); i_lo++)
    {
        int i_glo = desc.indx_l2g_r(i_lo);
        atbasis_row.get_local_index(i_glo, I, iI);
        for (int j_lo = 0; j_lo != desc.n_loc(); j_lo++)
        {
            int j_glo = desc.indx_l2g_c(j_lo);
            atbasis_col.get_local_index(j_glo, J, jJ);
            if (IJmap.count(I) == 0 || IJmap.at(I).count(J) == 0 || IJmap.at(I).at(J).size() == 0)
                IJmap[I][J] = matrix_m<T>{static_cast<int>(atbasis_row.get_atom_nb(I)), static_cast<int>(atbasis_col.get_atom_nb(J)), major_map};
            IJmap[I][J](iI, jJ) = mat_lo(i_lo, j_lo);
        }
    }
}


template <typename T>
void map_block_to_IJ_storage_new(map<int, map<int, matrix_m<T>>> &IJmap,
                                 const LIBRPA::AtomicBasis &atbasis,
                                 const map<int, vector<int>> &map_lor_I_is,
                                 const map<int, vector<int>> &map_loc_J_js,
                                 const matrix_m<T> &mat_lo, const LIBRPA::Array_Desc &desc, MAJOR major_map)
{
    // map<int, map<int, matrix_m<T>>> IJmap_local;
    int max_threads = omp_get_max_threads(); // Get the total number of available threads

    std::vector<int> Is;
    std::vector<int> Js;
    for (const auto &I_is: map_lor_I_is)
    {
        Is.push_back(I_is.first);
    }
    for (const auto &J_js: map_loc_J_js)
    {
        Js.push_back(J_js.first);
    }
    const auto n_pairs_total = Is.size() * Js.size();

    // Create necessay atom-pair blocks
#pragma omp parallel for collapse(2) schedule(dynamic)
    for (int i_I = 0; i_I != Is.size(); i_I++)
    {
        for (int j_J = 0; j_J != Js.size(); j_J++)
        {
            const auto &I = Is[i_I];
            const auto &J = Js[j_J];
            const auto n_i = static_cast<int>(atbasis.get_atom_nb(I));
            const auto n_j = static_cast<int>(atbasis.get_atom_nb(J));
            auto mat = matrix_m<T>(n_i, n_j, major_map);
#pragma omp critical
            {
                // IJmap_local[I][J] = std::move(mat);
                IJmap[I][J] = std::move(mat);
            }
        }
    }

    int subdiv = max_threads / n_pairs_total;
    if (subdiv * n_pairs_total < max_threads)
    {
        subdiv += 1;
    }

    // Compute starting and ending j indices
    std::vector<std::array<int, 4>> blocks;
    for (const auto &I: Is)
    {
        for (const auto &J: Js)
        {
            const int n_j_loc = map_loc_J_js.at(J).size();
            const auto step = n_j_loc / subdiv;
            const auto resi = n_j_loc % subdiv;
            for (int j_subdiv = 0; j_subdiv < subdiv; j_subdiv++)
            {
                int st = j_subdiv * step + j_subdiv * int(j_subdiv < resi) + resi * int(j_subdiv >= resi);
                int ed = std::min(n_j_loc, st + step + int(j_subdiv < resi));
                blocks.push_back({I, J, st, ed});
            }
        }
    }

#pragma omp parallel for schedule(dynamic)
    for (int i_block = 0; i_block != blocks.size(); i_block++)
    {
        const auto &index = blocks[i_block];
        const auto &I = index[0];
        const auto &J = index[1];
        const auto &i_indices = map_lor_I_is.at(I);
        const auto &j_indices = map_loc_J_js.at(J);

        // auto mat = IJmap_local[I][J];
        auto &mat = IJmap[I][J];

        const auto &j_st = index[2];
        const auto &j_ed = index[3];
        for (int i_idx = 0; i_idx < i_indices.size(); i_idx++)
        {
            const int i = i_indices[i_idx];
            const auto i_loc = desc.indx_g2l_r(atbasis.get_global_index(I, i));
            for (int j_idx = j_st; j_idx < j_ed; j_idx++)
            {
                const int j = j_indices[j_idx];
                const auto l_loc = desc.indx_g2l_c(atbasis.get_global_index(J, j));
                mat(i, j) = mat_lo(i_loc, l_loc);
            }
        }
    }

    // Merge local map into the global map
    // IJmap.merge(IJmap_local);
}

/*!
 * @brief Compute power of Hermitian matrix using BLACS
 *
 * @param  [in,out]  A_local       Process-local part of global matrix A to be powered
 * @param  [in]      ad_A          Array descriptor of A
 * @param  [out]     A_local       Process-local part of the eigenvector matrix (Z) of A
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
                                            const LIBRPA::Array_Desc &ad_A,
                                            matrix_m<std::complex<T>> &Z_local,
                                            const LIBRPA::Array_Desc &ad_Z,
                                            size_t &n_filtered, T *W, T power,
                                            const T &threshold = -1.e5)
{
    using LIBRPA::envs::ofs_myid;

    assert (A_local.is_col_major() && Z_local.is_col_major());
    const bool is_int_power = fabs(power - int(power)) < 1e-4;
    const int n = ad_A.m();
    const char jobz = 'V';
    const char uplo = 'U';

    // temporary array for heev with optmized block size
    const int blocksize_row_opt = min(ad_A.mb(), 128);
    const int blocksize_col_opt = min(ad_A.nb(), 128);

    // initialize descriptor of array A for optimized block size
    LIBRPA::Array_Desc ad_A_opt(ad_A.ictxt());
    ad_A_opt.init(n, n, blocksize_row_opt, blocksize_col_opt, 0, 0);
    auto A_local_opt = init_local_mat<std::complex<T>>(ad_A_opt, MAJOR::COL);
    // NOTE: imply A and Z should be in the same context
    ScalapackConnector::pgemr2d_f(n, n, A_local.ptr(), 1, 1, ad_A.desc,
                                  A_local_opt.ptr(), 1, 1, ad_A_opt.desc, ad_A.ictxt());

    // initialize descriptor of array Z for optimized block size
    if (ad_A.ictxt() != ad_Z.ictxt())
    {
        ofs_myid << "Warning(power_hemat_blacs): input contexts of A and Z are different!" << endl;
    }
    LIBRPA::Array_Desc ad_Z_opt(ad_Z.ictxt());
    ad_Z_opt.init(n, n, blocksize_row_opt, blocksize_col_opt, 0, 0);
    auto Z_local_opt = init_local_mat<std::complex<T>>(ad_Z_opt, MAJOR::COL);
    // printf("Z_local_opt size: %d\n", Z_local_opt.size());

    Profiler::start("power_hemat_blacs_1");
    int lwork = -1, lrwork = -1, info = 0;
    std::complex<T>  *work;
    T *rwork;
    {
        work  = new std::complex<T>[1];
        rwork = new T[1];
        // query the optimal lwork and lrwork
        T *Wquery = new T[1];
        // LIBRPA::utils::lib_printf("power_hemat_blacs descA %s\n", ad_A.info_desc().c_str());
        ScalapackConnector::pheev_f(jobz, uplo,
                n, A_local_opt.ptr(), 1, 1, ad_A_opt.desc,
                Wquery, Z_local_opt.ptr(), 1, 1, ad_A_opt.desc, work, lwork, rwork, lrwork, info);
        lwork = int(work[0].real());
        lrwork = int(rwork[0]);
        delete [] work;
        delete [] Wquery;
        delete [] rwork;
    }
    Profiler::stop("power_hemat_blacs_1");

    Profiler::start("power_hemat_blacs_2");
    work = new std::complex<T> [lwork];
    rwork = new T [lrwork];
    ScalapackConnector::pheev_f(jobz, uplo,
            n, A_local_opt.ptr(), 1, 1, ad_A_opt.desc,
            W, Z_local_opt.ptr(), 1, 1, ad_Z_opt.desc, work, lwork, rwork, lrwork, info);
    delete [] work;
    delete [] rwork;
    Profiler::stop("power_hemat_blacs_2");
    // Optimized A no longer used
    Profiler::start("power_hemat_blacs_3");
    A_local_opt.clear();
    // send back the eigenvector matrix
    ScalapackConnector::pgemr2d_f(n, n, Z_local_opt.ptr(), 1, 1, ad_Z_opt.desc,
                                  Z_local.ptr(), 1, 1, ad_Z.desc, ad_Z.ictxt());
    Profiler::stop("power_hemat_blacs_3");

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
    Profiler::start("power_hemat_blacs_4");
    T W_temp[n];
    for (int i = 0; i != n_filtered; i++)
        W_temp[i] = 0.0;
    for (int i = n_filtered; i != n; i++)
    {
        if (W[i] < 0 && !is_int_power)
        {
            LIBRPA::utils::lib_printf("Warning! unfiltered negative eigenvalue with non-integer power: # %d ev = %f , pow = %f\n", i, W[i], power);
        }
        if (fabs(W[i]) < 1e-10 && power < 0)
        {
            LIBRPA::utils::lib_printf("Warning! unfiltered nearly-singular eigenvalue with negative power: # %d ev = %f , pow = %f\n", i, W[i], power);
        }
        W_temp[i] = std::pow(W[i], power);
    }
    Profiler::stop("power_hemat_blacs_4");
    // debug print
    // for (int i = 0; i != n; i++)
    // {
    //     LIBRPA::utils::lib_printf("%d %f %f\n", i, W[i], W_temp[i]);
    // }

    Profiler::start("power_hemat_blacs_5");
    // create scaled eigenvectors
    auto scaled_opt = Z_local_opt.copy();
    for (int i = 0; i != n; i++)
    {
        ScalapackConnector::pscal_f(n, W_temp[i], scaled_opt.ptr(), 1, 1+i, ad_Z_opt.desc, 1);
    }
    ScalapackConnector::pgemm_f('N', 'C', n, n, n, 1.0, Z_local_opt.ptr(), 1, 1, ad_Z_opt.desc, scaled_opt.ptr(), 1, 1, ad_Z_opt.desc, 0.0, A_local.ptr(), 1, 1, ad_A.desc);
    auto scaled = Z_local.copy();
    // send back the scaled eigenvector matrix with descriptor using optimized block size to that with input descriptor
    ScalapackConnector::pgemr2d_f(n, n, scaled_opt.ptr(), 1, 1, ad_Z_opt.desc,
                                  scaled.ptr(), 1, 1, ad_Z.desc, ad_Z.ictxt());
    Profiler::stop("power_hemat_blacs_5");

    return scaled;
}

template <typename T, typename Treal = typename to_real<T>::type>
void print_matrix_mm_parallel(ostream &os, const matrix_m<T> &mat_loc, const LIBRPA::Array_Desc &ad, Treal threshold = 1e-15, bool row_first = true)
{
    LIBRPA::Array_Desc ad_fb(ad.ictxt());
    const int nr = ad.m(), nc = ad.n();
    const int irsrc = ad.irsrc(), icsrc = ad.icsrc();
    ad_fb.init(nr, nc, nr, nc, irsrc, icsrc);
    matrix_m<T> mat_glo = init_local_mat<T>(ad_fb, mat_loc.major());

    ScalapackConnector::pgemr2d_f(nr, nc, mat_loc.ptr(), 1, 1, ad.desc, mat_glo.ptr(), 1, 1, ad_fb.desc, ad.ictxt());
    if (ad_fb.is_src() && os.good())
    {
        print_matrix_mm(mat_glo, os, threshold, row_first);
    }
}

template <typename T, typename Treal = typename to_real<T>::type>
void print_matrix_mm_file_parallel(const char *fn, const matrix_m<T> &mat_loc, const LIBRPA::Array_Desc &ad, Treal threshold = 1e-15, bool row_first = true)
{
    ofstream fs;
    if (ad.is_src())
        fs.open(fn);
    print_matrix_mm_parallel(fs, mat_loc, ad, threshold, row_first);
}

template <typename T>
matrix_m<T> multiply_scalapack(const matrix_m<T> &m1_loc, const LIBRPA::Array_Desc &desc_m1,
                               const matrix_m<T> &m2_loc, const LIBRPA::Array_Desc &desc_m2,
                               const LIBRPA::Array_Desc &desc_prod)
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
void invert_scalapack(matrix_m<T> &m_loc, const LIBRPA::Array_Desc &desc_m)
{
    assert(m_loc.is_col_major());
    int info = 0;
    if (m_loc.is_col_major())
    {
        // NOTE: use local leading dimension desc_m.lld() will lead to invalid pointer error
        int *ipiv = new int [desc_m.m()];
        // LIBRPA::utils::lib_printf("invert_scalpack desc %s\n", desc_m.info_desc().c_str());
        ScalapackConnector::pgetrf_f(desc_m.m(), desc_m.n(), m_loc.ptr(), 1, 1, desc_m.desc, ipiv, info);
        // get optimized work size
        int lwork = -1, liwork = -1;
        {
            T *work = new T [1];
            int *iwork = new int [1];
            ScalapackConnector::pgetri_f(desc_m.m(), m_loc.ptr(), 1, 1, desc_m.desc, ipiv, work, lwork, iwork, liwork, info);
            lwork = int(get_real(work[0]));
            liwork = iwork[0];
            // LIBRPA::utils::lib_printf("lwork %d liwork %d\n", lwork, liwork);
            delete [] work;
            delete [] iwork;
        }

        T *work = new T[lwork];
        int *iwork = new int[liwork];
        ScalapackConnector::pgetri_f(desc_m.m(), m_loc.ptr(), 1, 1, desc_m.desc, ipiv, work, lwork, iwork, liwork, info);
        delete [] iwork;
        delete [] work;
        delete [] ipiv;
        // LIBRPA::utils::lib_printf("done free\n");
    }
    else
    {
        throw std::logic_error("row-major invert is not implemented");
    }
}
