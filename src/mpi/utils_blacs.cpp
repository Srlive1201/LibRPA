#include "utils_blacs.h"

namespace librpa_int
{

namespace utils
{

std::pair<ArrayDesc, ArrayDesc> prepare_array_desc_mr2d_src_and_all(const BlacsCtxtHandler &ctxt_h, const int &m, const int &n, const int &mb, const int &nb, const int &irsrc, const int &icsrc)
{
    ArrayDesc desc_fullblk_on_src(ctxt_h), desc_all_procs(ctxt_h);
    desc_fullblk_on_src.init(m, n, m, n, irsrc, icsrc);
    desc_all_procs.init(m, n, mb, nb, irsrc, icsrc);
    return {desc_fullblk_on_src, desc_all_procs};
}

} /* end of namespace utils */

} /* end of namespace librpa_int */
