#include "utils_blacs.h"

namespace LIBRPA
{

namespace utils
{

std::pair<Array_Desc, Array_Desc> prepare_array_desc_mr2d_src_and_all(const BLACS_CTXT_handler &ctxt_h, const int &m, const int &n, const int &mb, const int &nb, const int &irsrc, const int &icsrc)
{
    Array_Desc desc_fullblk_on_src(ctxt_h), desc_all_procs(ctxt_h);
    desc_fullblk_on_src.init(m, n, m, n, irsrc, icsrc);
    desc_all_procs.init(m, n, mb, nb, irsrc, icsrc);
    return {desc_fullblk_on_src, desc_all_procs};
}

} /* end of namespace utils */

} /* end of namespace LIBRPA */
