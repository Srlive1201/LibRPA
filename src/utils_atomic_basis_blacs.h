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
 * @param  [in]  myid                  Process rank to inquire
 * @param  [in]  map_proc_IJs_avail    Available atom-pair blocks on all processes, with process rank as map key
 * @param  [in]  atbasis_r             AtomicBasis object for row
 * @param  [in]  atbasis_c             AtomicBasis object for column
 * @param  [in]  ad                    Array descriptor for BLACS distribution
 * @param  [in]  row_fast              Flag to set the row basis index faster
 * @param  [in]  row_major             flag to compute the 1D indices from 2D in row-major (column major if false)
 *
 * @retval    map_send, map_recv
 *
 *            map_send: process id -> [ global indices (in atom-pair context) ]
 *
 *            map_send: process id -> [ global indices (in BLACS context) ]
 */
std::pair<std::map<int, std::vector<size_t>>, std::map<int, std::vector<size_t>>>
get_communicate_global_ids_list_ap_to_blacs(const int &myid,
                                            const std::map<int, std::vector<atpair_t>> &map_proc_IJs_avail,
                                            const AtomicBasis &atbasis_r,
                                            const AtomicBasis &atbasis_c,
                                            const Array_Desc &ad, bool row_fast, bool row_major);

/*!
 * @brief Similar to get_communicate_global_ids_list_ap_to_blacs, but the atom-pair data are only
 *        available with unordered atom pairs, i.e. either upper or lower half.
 *        This is intended for symmetric or Hermitian matrices.
 *
 * @param  [in]  myid                  Process rank to inquire
 * @param  [in]  uplo                  whether the atom pairs data are available in upper ('U', 'u') or lower ('L', 'l') half
 * @param  [in]  map_proc_IJs_avail    Available atom-pair blocks on all processes, with process rank as map key
 *                                     The user should ensure that either lower or upper half atom pairs are in the map,
 *                                     consistent with uplo.
 * @param  [in]  atbasis               AtomicBasis object for row and column
 * @param  [in]  ad                    Array descriptor for BLACS distribution
 * @param  [in]  row_fast              Flag to set the row basis index faster
 * @param  [in]  row_major             flag to compute the 1D indices from 2D in row-major (column major if false)
 *
 * @retval    map_send, map_recv
 *
 *            map_send: process id -> [ global indices (in atom-pair context) ]
 *
 *            map_recv: process id -> { [ global indices (in BLACS context) ], [ if conjugate should be taken when applicable ] }
 */
std::pair<std::map<int, std::vector<size_t>>, std::map<int, std::pair<std::vector<size_t>, std::vector<bool>>>>
get_communicate_global_ids_list_ap_to_blacs_sy(const int &myid,
                                               const char &uplo,
                                               const std::map<int, std::vector<atpair_t>> &map_proc_IJs_avail,
                                               const AtomicBasis &atbasis,
                                               const Array_Desc &ad, bool row_fast, bool row_major);

/*!
 * @brief Similar to get_communicate_global_ids_list_ap_to_blacs, but the indices are local instead of global.
 *
 * @param  [in]  myid                  Process rank to inquire
 * @param  [in]  map_proc_IJs_avail    Available atom-pair blocks on all processes, with process rank as map key
 * @param  [in]  atbasis_r             AtomicBasis object for row
 * @param  [in]  atbasis_c             AtomicBasis object for column
 * @param  [in]  ad                    Array descriptor for BLACS distribution
 * @param  [in]  row_fast              Flag to set the row basis index faster
 * @param  [in]  row_major             flag to compute the 1D indices from 2D in row-major (column major if false)
 *
 * @retval    map_send, map_recv.
 *
 *            map_send: process id -> [ { <I,J>, [ local indices in atom-pair sub-matrix ] } , ...]
 *
 *            map_recv: process id -> [ local indices in BLACS sub-matrix ]
 */
std::pair<std::map<int, std::vector<std::pair<atpair_t, std::vector<size_t>>>>,
          std::map<int, std::vector<size_t>>>
get_communicate_local_ids_list_ap_to_blacs(const int &myid,
                                           const std::map<int, std::vector<atpair_t>> &map_proc_IJs_avail,
                                           const AtomicBasis &atbasis_r,
                                           const AtomicBasis &atbasis_c,
                                           const Array_Desc &ad, bool row_fast, bool row_major);

/*!
 * @brief Similar to get_communicate_local_ids_list_ap_to_blacs, but the atom-pair data are only
 *        available with unordered atom pairs, i.e. either upper or lower half.
 *        This is intended for symmetric or Hermitian matrices.
 *
 * @param  [in]  myid                  Process rank to inquire
 * @param  [in]  uplo                  whether the atom pairs data are available in upper ('U', 'u') or lower ('L', 'l') half
 * @param  [in]  map_proc_IJs_avail    Available atom-pair blocks on all processes, with process rank as map key
 *                                     The user should ensure that either lower or upper half atom pairs are in the map,
 *                                     consistent with uplo.
 * @param  [in]  atbasis               AtomicBasis object for row and column
 * @param  [in]  ad                    Array descriptor for BLACS distribution
 * @param  [in]  row_fast              Flag to set the row basis index faster
 * @param  [in]  row_major             flag to compute the 1D indices from 2D in row-major (column major if false)
 *
 * @retval    map_send, map_recv.
 *
 *            map_send: process id -> [ { <I,J>, [ local indices in atom-pair sub-matrix ] } , ...]
 *
 *            map_recv: process id -> { [ local indices in BLACS sub-matrix ], [ if conjugate should be taken when applicable ] }
 */
std::pair<std::map<int, std::vector<std::pair<atpair_t, std::vector<size_t>>>>,
          std::map<int, std::pair<std::vector<size_t>, std::vector<bool>>>>
get_communicate_local_ids_list_ap_to_blacs_sy(const int &myid,
                                              const char &uplo,
                                              const std::map<int, std::vector<atpair_t>> &map_proc_IJs_avail,
                                              const AtomicBasis &atbasis,
                                              const Array_Desc &ad, bool row_fast, bool row_major);
} /* end of namespace utils */

} /* end of namespace LIBRPA */
