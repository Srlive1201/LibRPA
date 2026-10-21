#pragma once

#include "../io/output_matrix.h"
#include "atomic_basis.h"
#include "pbc.h"

namespace librpa_int
{

inline std::vector<RSpaceMatrixBlock> make_rspace_matrix_blocks(
    const ap_p_map<std::map<Vector3_Order<int>, Matz>> &blocks, const AtomicBasis &basis,
    const PeriodicBoundaryData &pbc)
{
    std::vector<RSpaceMatrixBlock> output;
    for (const auto &[atoms, translations] : blocks)
    {
        const auto nr = basis.get_atom_nb(atoms.first);
        const auto nc = basis.get_atom_nb(atoms.second);
        for (const auto &[translation, values] : translations)
        {
            assert(values.nr() == as_int(nr) && values.nc() == as_int(nc));
            output.push_back(
                {as_size(pbc.get_R_index(translation)), atoms.first, atoms.second, values});
        }
    }
    return output;
}

}  // namespace librpa_int
