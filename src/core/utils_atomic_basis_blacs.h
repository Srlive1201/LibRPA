#pragma once

#include <set>
#include <unordered_map>

#include "../mpi/base_blacs.h"
#include "atomic_basis.h"

namespace librpa_int
{

// aliases for mapping indices from process rank
typedef std::unordered_map<int, std::vector<size_t>> myid_blacsid_map_t;
typedef std::unordered_map<int, std::vector<size_t>> myid_ap_glid_map_t;
typedef std::unordered_map<int, std::vector<std::pair<atpair_t, std::vector<size_t>>>> myid_ap_loid_map_t;
typedef std::unordered_map<int, std::pair<std::vector<size_t>, std::vector<char>>> myid_blacs_idop_map_t;

//! obtain the necessary atom pair of atomic basis to build the block-cyclic submatrix
std::set<std::pair<int, int>> get_necessary_IJ_from_block_2D(const AtomicBasis &atbasis_row, const AtomicBasis &atbasis_col, const ArrayDesc& arrdesc);

std::set<std::pair<int, int>> get_necessary_IJ_from_block_2D_sy(const char &uplo, const AtomicBasis &atbasis, const ArrayDesc& arrdesc);

// Scheduler for index and communication for conversion between atom-pair and BLACS distribution
class IndexScheduler
{
private:
    bool initialized_;
public:
    // atom pairs that this process should handle
    std::vector<atpair_t> atpairs;

    // Local 1D indices in atomic basis context to send when converting from AP to BLACS
    std::vector<size_t> ids_ap_ipair;
    std::vector<size_t> ids_ap_locid;
    // Local 1D indices in BLACS context to send when converting from BLACS to AP
    std::vector<size_t> ids_blacs_locid;

    // Displacements and counts to send for AP->BLACS
    MPI_Count total_count_ap;
    std::vector<MPI_Count> disp_ap;
    std::vector<MPI_Count> counts_ap;
    // Displacements and counts to send for BLACS->AP
    MPI_Count total_count_blacs;
    std::vector<MPI_Count> disp_blacs;
    std::vector<MPI_Count> counts_blacs;

    // default constructor
    IndexScheduler()
        : initialized_(false),
          atpairs(),
          ids_ap_ipair(),
          ids_ap_locid(),
          ids_blacs_locid(),
          total_count_ap(0),
          disp_ap(),
          counts_ap(),
          total_count_blacs(0),
          disp_blacs(),
          counts_blacs() {};
    // destructor
    ~IndexScheduler() {};
    // copy & move constructor
    IndexScheduler(const IndexScheduler& s) = delete;
    IndexScheduler(IndexScheduler&& s) = delete;
    // copy & move assignment operator
    IndexScheduler& operator=(const IndexScheduler& s) = delete;
    IndexScheduler& operator=(IndexScheduler&& s) = delete;

    inline bool initialized() const noexcept { return initialized_; }

    void init(const std::map<int, std::set<atpair_t>> &map_proc_IJs,
              const AtomicBasis &atbasis_r, const AtomicBasis &atbasis_c,
              const ArrayDesc &ad, const bool row_major) noexcept;

    void init(const std::unordered_map<int, std::set<atpair_t>> &map_proc_IJs,
              const AtomicBasis &atbasis_r, const AtomicBasis &atbasis_c,
              const ArrayDesc &ad, const bool row_major) noexcept;

    void reset();

private:
    void init_impl_(const ap_p_map<int> &map_IJ_proc,
                    const AtomicBasis &atbasis_r, const AtomicBasis &atbasis_c,
                    const ArrayDesc &ad, const bool row_major) noexcept;
};

std::unordered_map<int, std::set<atpair_t>>
get_balanced_ap_distribution_for_consec_descriptor(const AtomicBasis &atbasis_r,
                                                   const AtomicBasis &atbasis_c,
                                                   const ArrayDesc &ad, bool create_empty = false);

/*!
 * @brief Get the list of global matrix indices (1D) to communicate for
 *        conversion from atom-pair to BLACS block-cyclic distributions
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
 *            map_recv: process id -> [ global indices (in BLACS context) ]
 */
std::pair<myid_ap_glid_map_t, myid_blacsid_map_t>
get_communicate_global_ids_list_ap_to_blacs(const int &myid,
                                            const std::unordered_map<int, std::vector<atpair_t>> &map_proc_IJs_avail,
                                            const AtomicBasis &atbasis_r,
                                            const AtomicBasis &atbasis_c,
                                            const ArrayDesc &ad,
                                            bool row_fast, bool row_major, bool include_self = false);

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
 *            map_recv: process id -> { [ global indices (in BLACS context) ], [ conjugate label ] }.
 *                      conjugate label is 'c' is necessary conjugate should be performed.
 */
std::pair<myid_ap_glid_map_t, myid_blacs_idop_map_t>
get_communicate_global_ids_list_ap_to_blacs_sy(const int &myid,
                                               const char &uplo,
                                               const std::unordered_map<int, std::vector<atpair_t>> &map_proc_IJs_avail,
                                               const AtomicBasis &atbasis,
                                               const ArrayDesc &ad,
                                               bool row_fast, bool row_major, bool include_self = false);

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
std::pair<myid_ap_loid_map_t, myid_blacsid_map_t>
get_communicate_local_ids_list_ap_to_blacs(const int &myid,
                                           const std::unordered_map<int, std::vector<atpair_t>> &map_proc_IJs_avail,
                                           const AtomicBasis &atbasis_r,
                                           const AtomicBasis &atbasis_c,
                                           const ArrayDesc &ad,
                                           bool row_fast, bool row_major, bool include_self = false);

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
 *            map_recv: process id -> { [ local indices in BLACS sub-matrix ], [ conjugate label ] }.
 *                      conjugate label is 'c' is necessary conjugate should be performed.
 */
std::pair<myid_ap_loid_map_t, myid_blacs_idop_map_t>
get_communicate_local_ids_list_ap_to_blacs_sy(const int &myid,
                                              const char &uplo,
                                              const std::unordered_map<int, std::vector<atpair_t>> &map_proc_IJs_avail,
                                              const AtomicBasis &atbasis,
                                              const ArrayDesc &ad,
                                              bool row_fast, bool row_major, bool include_self = false);

/*!
 * @brief Get the list of global matrix indices (1D) to communicate for
 *        conversion from BLACS block-cyclic to atom-pair distributions
 *
 * @param  [in]  myid                  Process rank to inquire
 * @param  [in]  map_proc_IJs_require  Required atom-pair blocks on all processes, with process rank as map key
 * @param  [in]  atbasis_r             AtomicBasis object for row
 * @param  [in]  atbasis_c             AtomicBasis object for column
 * @param  [in]  ad                    Array descriptor for BLACS distribution
 * @param  [in]  row_fast              Flag to set the row basis index faster
 * @param  [in]  row_major             flag to compute the 1D indices from 2D in row-major (column major if false)
 *
 * @retval    map_send, map_recv
 *
 *            map_send: process id -> [ global indices (in BLACS context) ]
 *
 *            map_recv: process id -> [ global indices (in atom-pair context) ]
 */
std::pair<myid_blacsid_map_t, myid_ap_glid_map_t>
get_communicate_global_ids_list_blacs_to_ap(const int &myid,
                                            const std::unordered_map<int, std::vector<atpair_t>> &map_proc_IJs_require,
                                            const AtomicBasis &atbasis_r,
                                            const AtomicBasis &atbasis_c,
                                            const ArrayDesc &ad,
                                            bool row_fast, bool row_major, bool include_self = false);

/*!
 * @brief Similar to get_communicate_global_ids_list_blacs_to_ap, but the indices are local instead of global.
 *
 * @param  [in]  myid                  Process rank to inquire
 * @param  [in]  map_proc_IJs_require  Required atom-pair blocks on all processes, with process rank as map key
 * @param  [in]  atbasis_r             AtomicBasis object for row
 * @param  [in]  atbasis_c             AtomicBasis object for column
 * @param  [in]  ad                    Array descriptor for BLACS distribution
 * @param  [in]  row_fast              Flag to set the row basis index faster
 * @param  [in]  row_major             flag to compute the 1D indices from 2D in row-major (column major if false)
 *
 * @retval    map_send, map_recv
 *
 *            map_send: process id -> [ local indices in BLACS sub-matrix ]
 *
 *            map_recv: process id -> [ { <I,J>, [ local indices in atom-pair sub-matrix ] } , ...]
 */
std::pair<myid_blacsid_map_t, myid_ap_loid_map_t>
get_communicate_local_ids_list_blacs_to_ap(const int &myid,
                                           const std::unordered_map<int, std::vector<atpair_t>> &map_proc_IJs_require,
                                           const AtomicBasis &atbasis_r,
                                           const AtomicBasis &atbasis_c,
                                           const ArrayDesc &ad,
                                           bool row_fast, bool row_major, bool include_self = false);

} /* end of namespace librpa_int */
