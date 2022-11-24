#pragma once
#include "matrix_m.h"
#include "parallel_mpi.h"

template <typename T>
inline matrix_m<T> init_local_mat(const LIBRPA::Array_Desc &ad, MAJOR major)
{
    matrix_m<T> mat_lo(ad.m_loc(), ad.n_loc(), major);
    if (mat_lo.size() == 0)
        mat_lo.resize(1, 1);
    return mat_lo;
}

template <typename T>
inline matrix_m<T> get_local_mat(const matrix_m<T> &mat_go, const LIBRPA::Array_Desc &ad, MAJOR major)
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
inline matrix_m<T> get_local_mat(const T *pv, MAJOR major_pv, const LIBRPA::Array_Desc &ad, MAJOR major)
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
inline void collect_block_from_IJ_storage(
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

