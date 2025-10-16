#pragma once

#include "base_blacs.h"
#include "atomic_basis.h"

#include <set>

namespace LIBRPA
{

namespace utils
{

//! prepare array descriptors for distributing(collecting) submatrices
//! from(to) a full matrix on source process with p?gemr2d
std::pair<Array_Desc, Array_Desc> prepare_array_desc_mr2d_src_and_all(
    const BLACS_CTXT_handler &ctxt_h, const int &m, const int &n, const int &mb,
    const int &nb, const int &irsrc, const int &icsrc);

//! obtain the necessary atom pair of atomic basis to build the block-cyclic submatrix
std::set<std::pair<int, int>> get_necessary_IJ_from_block_2D(const AtomicBasis &atbasis_row, const AtomicBasis &atbasis_col, const Array_Desc& arrdesc);

std::set<std::pair<int, int>> get_necessary_IJ_from_block_2D_sy(const char &uplo, const AtomicBasis &atbasis, const Array_Desc& arrdesc);

} /* end of namespace utils */

} /* end of namespace LIBRPA */
