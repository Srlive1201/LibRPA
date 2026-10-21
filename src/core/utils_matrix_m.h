#pragma once

#include "../math/utils_matrix_m.h"
#include "atomic_basis.h"

namespace librpa_int
{

template <typename T>
void fill_ap_block_from_global(matrix_m<T> &block, const matrix_m<T> &global, const atpair_t &atoms,
                               const AtomicBasis &row_basis, const AtomicBasis &column_basis)
{
    assert(as_size(block.nr()) == row_basis.get_atom_nb(as_int(atoms.first)));
    assert(as_size(block.nc()) == column_basis.get_atom_nb(as_int(atoms.second)));
    fill_block_from_global(block, global, row_basis.get_part_range()[atoms.first],
                           column_basis.get_part_range()[atoms.second]);
}

template <typename T>
matrix_m<T> get_ap_block_from_global(const matrix_m<T> &global, const atpair_t &atoms,
                                     const AtomicBasis &row_basis, const AtomicBasis &column_basis,
                                     MAJOR major = MAJOR::AUTO)
{
    if (major == MAJOR::AUTO) major = global.major();
    matrix_m<T> block(as_int(row_basis.get_atom_nb(as_int(atoms.first))),
                      as_int(column_basis.get_atom_nb(as_int(atoms.second))), major);
    block.zero_out();
    fill_ap_block_from_global(block, global, atoms, row_basis, column_basis);
    return block;
}

}  // namespace librpa_int
