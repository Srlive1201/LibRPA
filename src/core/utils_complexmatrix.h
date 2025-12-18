#pragma once

#include "../math/complexmatrix.h"
#include "atomic_basis.h"

namespace librpa_int
{

ComplexMatrix get_ap_block_from_global(const ComplexMatrix &cm_global,
                                       const atpair_t &IJ,
                                       const AtomicBasis &atbasis_r, const AtomicBasis &atbasis_c);

}
