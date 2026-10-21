#pragma once

#include "matrix_m.h"

namespace librpa_int
{

template <typename T>
void fill_block_from_global(matrix_m<T> &m_ap, const matrix_m<T> &m_global, const int ir_st,
                            const int ic_st)
{
    const auto nI = m_ap.nr();
    const auto nJ = m_ap.nc();
    assert(ir_st >= 0 && ir_st + nI <= m_global.nr());
    assert(ic_st >= 0 && ic_st + nJ <= m_global.nc());
    const auto major = m_global.major();
    if (major == MAJOR::ROW)
    {
        for (int i = 0; i < nI; i++)
        {
            for (int j = 0; j < nJ; j++)
            {
                m_ap(i, j) = m_global(ir_st + i, ic_st + j);
            }
        }
    }
    else
    {
        for (int j = 0; j < nJ; j++)
        {
            for (int i = 0; i < nI; i++)
            {
                m_ap(i, j) = m_global(ir_st + i, ic_st + j);
            }
        }
    }
}

} /* end of namespace librpa_int */
