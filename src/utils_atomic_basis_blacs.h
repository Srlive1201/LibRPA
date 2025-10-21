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

/*!
 * @brief Get the list of global matrix indices (1D) to communicate during
 *        conversion between atom-pair and BLACS block-cyclic distributions
 *
 * @retval    map_send, map_receive
 */
std::pair<std::map<int, std::vector<size_t>>, std::map<int, std::vector<size_t>>>
get_communicate_global_ids_list_ap_to_blacs(const int &myid,
                                            const std::map<int, std::vector<atpair_t>> &map_proc_IJs_avail,
                                            const AtomicBasis &atbasis_r,
                                            const AtomicBasis &atbasis_c,
                                            const Array_Desc &ad, bool row_fast, bool row_major);

/*!
 * @brief Similar to get_communicate_global_ids_list_ap_to_blacs,
 *        but the indices are local instead of global.
 *
 * @retval    map_send, map_recv
 *            map_send: process id -> [ { <I,J>, [ local indices in atom-pair sub-matrix ] } , ...]
 *            map_recv: process id -> [ local indices in BLACS sub-matrix ]
 */
std::pair<std::map<int, std::vector<std::pair<atpair_t, std::vector<size_t>>>>,
          std::map<int, std::vector<size_t>>>
get_communicate_local_ids_list_ap_to_blacs(const int &myid,
                                           const std::map<int, std::vector<atpair_t>> &map_proc_IJs_avail,
                                           const AtomicBasis &atbasis_r,
                                           const AtomicBasis &atbasis_c,
                                           const Array_Desc &ad, bool row_fast, bool row_major);

} /* end of namespace utils */

} /* end of namespace LIBRPA */
