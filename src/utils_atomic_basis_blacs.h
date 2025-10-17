#pragma once

#include <set>

#include "atomic_basis.h"
#include "base_blacs.h"

namespace LIBRPA
{

namespace utils
{

//! obtain the necessary atom pair of atomic basis to build the block-cyclic submatrix
std::set<std::pair<int, int>> get_necessary_IJ_from_block_2D(const AtomicBasis &atbasis_row, const AtomicBasis &atbasis_col, const Array_Desc& arrdesc);

std::set<std::pair<int, int>> get_necessary_IJ_from_block_2D_sy(const char &uplo, const AtomicBasis &atbasis, const Array_Desc& arrdesc);

} /* end of namespace utils */

} /* end of namespace LIBRPA */
