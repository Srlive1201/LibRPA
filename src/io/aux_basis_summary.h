#pragma once

#include <cstddef>
#include <map>
#include <string>
#include <vector>

#include "../core/atom.h"

namespace librpa_int
{

std::string format_aux_basis_compression_summary(const std::vector<std::size_t> &large_nbs,
                                                 const std::vector<std::size_t> &small_nbs,
                                                 const std::map<atom_t, int> &atom_to_type);

}
