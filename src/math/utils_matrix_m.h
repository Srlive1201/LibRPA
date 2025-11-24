#pragma once

#include "matrix_m.h"

#include "../core/atomic_basis.h"

namespace LIBRPA
{

namespace utils
{

template <typename T>
void fill_ap_block_from_global(matrix_m<T> &m_ap, const matrix_m<T> &m_global,
                               const atpair_t &IJ,
                               const AtomicBasis &atbasis_r, const AtomicBasis &atbasis_c)
{
    const auto nI = m_ap.nr();
    const auto nJ = m_ap.nc();
    // assert (m_ap.major() == m_global.major());
    assert (as_size(nI) == atbasis_r.get_atom_nb(as_int(IJ.first)));
    assert (as_size(nJ) == atbasis_c.get_atom_nb(as_int(IJ.second)));
    const auto ir_st = atbasis_r.get_part_range()[IJ.first];
    const auto ic_st = atbasis_c.get_part_range()[IJ.second];
    const auto major = m_global.major();
    if (major == MAJOR::ROW)
    {
        for (int i = 0; i < nI; i++)
        {
            for (int j = 0; j < nJ; j++)
            {
                m_ap(i, j) = m_global(ir_st+i, ic_st+j);
            }
        }
    }
    else
    {
        for (int j = 0; j < nJ; j++)
        {
            for (int i = 0; i < nI; i++)
            {
                m_ap(i, j) = m_global(ir_st+i, ic_st+j);
            }
        }
    }
}

template <typename T>
matrix_m<T> get_ap_block_from_global(const matrix_m<T> &m_global,
                                     const atpair_t &IJ,
                                     const AtomicBasis &atbasis_r, const AtomicBasis &atbasis_c,
                                     MAJOR major = MAJOR::AUTO)
{
    const auto &nI = atbasis_r.get_atom_nb(as_int(IJ.first));
    const auto &nJ = atbasis_c.get_atom_nb(as_int(IJ.second));
    if (major == MAJOR::AUTO)
    {
        major = m_global.major();
    }
    matrix_m<T> m_ap(as_int(nI), as_int(nJ), major);
    m_ap.zero_out();
    fill_ap_block_from_global(m_ap, m_global, IJ, atbasis_r, atbasis_c);
    return m_ap;
}

} /* end of namespace utils */

} /* end of namespace LIBRPA */

