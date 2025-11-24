#pragma once

#include "base_blacs.h"

namespace LIBRPA
{

namespace utils
{

//! prepare array descriptors for distributing(collecting) submatrices
//! from(to) a full matrix on source process with p?gemr2d
std::pair<ArrayDesc, ArrayDesc> prepare_array_desc_mr2d_src_and_all(
    const BlacsCtxtHandler &ctxt_h, const int &m, const int &n, const int &mb,
    const int &nb, const int &irsrc, const int &icsrc);

} /* end of namespace utils */

} /* end of namespace LIBRPA */
