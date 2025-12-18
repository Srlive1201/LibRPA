#include "utils_complexmatrix.h"

namespace librpa_int
{

ComplexMatrix get_ap_block_from_global(const ComplexMatrix &cm_global,
                                       const atpair_t &IJ,
                                       const AtomicBasis &atbasis_r, const AtomicBasis &atbasis_c)
{
    const auto I = IJ.first;
    const auto J = IJ.first;
    const auto I_num = atbasis_r.get_atom_nb(I);
    const auto J_num = atbasis_c.get_atom_nb(J);
    ComplexMatrix cm_IJ(I_num, J_num);
    for (size_t i = 0; i != I_num; i++)
    {
        size_t i_glo = atbasis_r.get_global_index(I, i);
        for (size_t j = 0; j != J_num; j++)
        {
            size_t j_glo = atbasis_c.get_global_index(J, j);
            cm_IJ(i, j) = cm_global(i_glo, j_glo);
        }
    }
    return cm_IJ;
}

}
