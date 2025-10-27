#include "utils_atomic_basis_blacs.h"

#include <algorithm>
#include <cassert>
// #include <iostream>
#include <map>
#include <stdexcept>

#include "base_utility.h"
#include "scalapack_connector.h"
#include "profiler.h"

// #include "stl_io_helper.h"

namespace LIBRPA
{

namespace utils
{

std::set<std::pair<int, int>> get_necessary_IJ_from_block_2D(const AtomicBasis &atbasis_row, const AtomicBasis &atbasis_col, const Array_Desc& arrdesc)
{
    std::set<std::pair<int, int>> IJs;
    if (as_size(arrdesc.m()) != atbasis_row.nb_total)
        throw std::invalid_argument("row basis and array descriptor inconsistent");
    if (as_size(arrdesc.n()) != atbasis_col.nb_total)
        throw std::invalid_argument("col basis and array descriptor inconsistent");

    const auto mlo = arrdesc.m_loc();
    const auto nlo = arrdesc.n_loc();
    size_t glo;
    int I, J;
    for (int ilo = 0; ilo != mlo; ilo++)
        for (int jlo = 0; jlo != nlo; jlo++)
        {
            glo = arrdesc.indx_l2g_r(ilo);
            I = atbasis_row.get_i_atom(glo);
            glo = arrdesc.indx_l2g_c(jlo);
            J = atbasis_col.get_i_atom(glo);
            IJs.insert({I, J});
        }
    return IJs;
}

std::set<std::pair<int, int>> get_necessary_IJ_from_block_2D_sy(const char &uplo, const AtomicBasis &atbasis, const Array_Desc& arrdesc)
{
    std::set<std::pair<int, int>> IJs;
    if (uplo != 'U' && uplo != 'L')
        throw std::invalid_argument("uplo should be U or L");
    if (as_size(arrdesc.m()) != atbasis.nb_total)
        throw std::invalid_argument("row basis and array descriptor inconsistent");
    if (as_size(arrdesc.n()) != atbasis.nb_total)
        throw std::invalid_argument("col basis and array descriptor inconsistent");

    const auto mlo = arrdesc.m_loc();
    const auto nlo = arrdesc.n_loc();
    size_t glo;
    int I, J;
    for (int ilo = 0; ilo != mlo; ilo++)
        for (int jlo = 0; jlo != nlo; jlo++)
        {
            glo = arrdesc.indx_l2g_r(ilo);
            I = atbasis.get_i_atom(glo);
            glo = arrdesc.indx_l2g_c(jlo);
            J = atbasis.get_i_atom(glo);
            if (uplo == 'L')
                I > J? IJs.insert({I, J}): IJs.insert({J, I});
            else
                I > J? IJs.insert({J, I}): IJs.insert({I, J});
        }
    return IJs;
}


// indices computation backend for AP<->BLACS conversion of general matrix
static 
std::pair<std::unordered_map<int, std::vector<std::pair<atpair_t, std::vector<std::pair<size_t, size_t>>>>>,
          std::unordered_map<int, std::vector<std::pair<size_t, size_t>>>>
_get_communicate_ids_list_ap_and_blacs(const int &myid,
                                       const std::unordered_map<int, std::vector<atpair_t>> &map_proc_IJs,
                                       const AtomicBasis &atbasis_r,
                                       const AtomicBasis &atbasis_c,
                                       const Array_Desc &ad, const bool row_fast,
                                       const bool return_local, const bool include_self = false)
{
    assert(atbasis_r.nb_total == as_size(ad.m()));
    assert(ad.m() > 0 && atbasis_r.nb_total == as_size(ad.m()));
    assert(ad.n() > 0 && atbasis_c.nb_total == as_size(ad.n()));

    // Objects to return
    // process ID -> list of basis indices (AP)
    std::unordered_map<int, std::vector<std::pair<atpair_t, std::vector<std::pair<size_t, size_t>>>>> map_proc_id2d_ap;
    // process ID -> list of basis indices (BLACS)
    std::unordered_map<int, std::vector<std::pair<size_t, size_t>>> map_proc_id2d_blacs;

    int myprow, mypcol;
    int dummy = 0;

    // basis indices availble at myid in the atom-pair distribution
    const std::vector<atpair_t> &IJs =
        map_proc_IJs.count(myid) ? map_proc_IJs.at(myid) : std::vector<atpair_t>();
    // if (myid == 0) std::cout << "IJs_avail " << IJs_avail << std::endl;
    Profiler::start("get_2d_mat_indices_atpair");
    const auto ids_ap = get_2d_mat_indices_atpair(atbasis_r, atbasis_c, IJs, row_fast);
    Profiler::stop("get_2d_mat_indices_atpair");

    // global basis indices required by BLACS
    Profiler::start("get_2d_mat_indices_blacs");
    const auto ids_blacs = get_2d_mat_indices_blacs(ad, myid, row_fast);
    Profiler::stop("get_2d_mat_indices_blacs");

    Profiler::start("ap_and_blacs_map_proc_map_id2d_ap");
    std::unordered_map<int, atom_mapping<std::vector<std::pair<size_t, size_t>>>::pair_t> map_proc_map_id2d_ap;
    std::unordered_map<int, std::vector<atpair_t>> map_proc_IJs_save_order;
    for (const auto &id: ids_ap)
    {
        const auto &i_gl = as_int(id.first);
        const auto &j_gl = as_int(id.second);
        // filter out the required BLACS indices that already owns by myid
        myprow = ScalapackConnector::indxg2p(i_gl, ad.mb(), dummy, ad.irsrc(), ad.nprows());
        mypcol = ScalapackConnector::indxg2p(j_gl, ad.nb(), dummy, ad.icsrc(), ad.npcols());
        const auto target = Cblacs_pnum(ad.ictxt(), myprow, mypcol);
        // if (myid == 0) std::cout << "id/target " << id << " " << target << std::endl;
        if (myid != target)
        {
            if (return_local)
            {
                int I, J;
                int i_lo, j_lo;
                atbasis_r.get_local_index(i_gl, I, i_lo);
                atbasis_r.get_local_index(j_gl, J, j_lo);
                auto &IJs_save_order = map_proc_IJs_save_order[target];
                const atpair_t IJ = {static_cast<atom_t>(I), static_cast<atom_t>(J)};
                if (std::find(IJs_save_order.cbegin(), IJs_save_order.cend(), IJ) == IJs_save_order.cend())
                {
                    IJs_save_order.emplace_back(IJ);
                }
                const atpair_t ij_lo = {static_cast<atom_t>(i_lo), static_cast<atom_t>(j_lo)};
                map_proc_map_id2d_ap[target][IJ].emplace_back(ij_lo);
            }
            else
            {
                // NOTE: global index requested, save all to the (0,0) dummy pair index
                const atpair_t IJ{0, 0};
                auto &IJs_save_order = map_proc_IJs_save_order[target];
                if (std::find(IJs_save_order.cbegin(), IJs_save_order.cend(), IJ) == IJs_save_order.cend())
                {
                    IJs_save_order.emplace_back(IJ);
                }
                map_proc_map_id2d_ap[target][IJ].emplace_back(id);
            }
        }
    }
    for (const auto &target_map: map_proc_map_id2d_ap)
    {
        const auto &target = target_map.first;
        const auto &IJs_save_order = map_proc_IJs_save_order.at(target);
        for (const auto &IJ: IJs_save_order)
        {
            map_proc_id2d_ap[target].push_back({IJ, target_map.second.at(IJ)});
        }
    }
    Profiler::stop("ap_and_blacs_map_proc_map_id2d_ap");

    Profiler::start("reverse_map_proc_IJs_avail");
    // reverse map of map_proc_IJs_avail
    atom_mapping<int>::pair_t map_IJs_proc;
    for (const auto &proc_IJs: map_proc_IJs)
    {
        const auto &proc = proc_IJs.first;
        for (const auto &IJ: proc_IJs.second)
        {
            map_IJs_proc[IJ] = proc;
        }
    }
    Profiler::stop("reverse_map_proc_IJs_avail");

    Profiler::start("map_proc_id2d_blacs");
    Cblacs_pcoord(ad.ictxt(), myid, &myprow, &mypcol);
    for (const auto &id: ids_blacs)
    {
        const auto I = static_cast<atom_t>(atbasis_r.get_i_atom(id.first));
        const auto J = static_cast<atom_t>(atbasis_c.get_i_atom(id.second));
        const atpair_t IJ({I, J});
        if (map_IJs_proc.count(IJ))
        {
            const auto &source = map_IJs_proc[IJ];
            // exclude self
            if (!include_self && source == myid) continue;
            if (return_local)
            {
                // FIXME: now only works when myid == ad.myid
                const auto i_lo = as_size(ad.indx_g2l_r(as_int(id.first)));
                const auto j_lo = as_size(ad.indx_g2l_c(as_int(id.second)));
                map_proc_id2d_blacs[source].push_back({i_lo, j_lo});
            }
            else
            {
                map_proc_id2d_blacs[source].emplace_back(id);
            }
        }
    }
    Profiler::stop("map_proc_id2d_blacs");

    return {map_proc_id2d_ap, map_proc_id2d_blacs};
}

// indices computation backend for AP->BLACS conversion of symmetric matrix
static
std::pair<std::unordered_map<int, std::vector<std::pair<atpair_t, std::vector<std::pair<size_t, size_t>>>>>,
          std::unordered_map<int, std::vector<std::pair<size_t, size_t>>>>
_get_communicate_ids_list_ap_to_blacs_sy(std::unordered_map<int, std::vector<char>> &map_lconj,
                                          const int &myid,
                                          const char &uplo,
                                          const std::unordered_map<int, std::vector<atpair_t>> &map_proc_IJs_avail,
                                          const AtomicBasis &atbasis,
                                          const Array_Desc &ad, const bool row_fast,
                                          const bool return_local, const bool include_self = false)
{
    assert(ad.m() > 0 && ad.m() == ad.n() && atbasis.nb_total == as_size(ad.m()));
    assert(uplo == 'U' || uplo == 'u' || uplo == 'L' || uplo == 'l');

    // Objects to return
    // process ID -> list of basis indices (AP)
    std::unordered_map<int, std::vector<std::pair<atpair_t, std::vector<std::pair<size_t, size_t>>>>> map_proc_id2d_ap;
    // process ID -> list of basis indices (BLACS)
    std::unordered_map<int, std::vector<std::pair<size_t, size_t>>> map_proc_id2d_blacs;

    const bool ap_upper_half = uplo == 'U' || uplo == 'u';
    // manually check the keys to ensure they are either all in upper half or all in lower half
    for (const auto &proc_IJs: map_proc_IJs_avail)
    {
        for (const auto &IJ: proc_IJs.second)
        {
            const auto &I = IJ.first;
            const auto &J = IJ.second;
            assert((ap_upper_half && J >= I) || (!ap_upper_half && I >= J));
        }
    }

    int myprow, mypcol;
    int myprow_T, mypcol_T;
    int dummy = 0;

    // basis indices availble at myid in the atom-pair distribution
    const std::vector<atpair_t> &IJs =
        map_proc_IJs_avail.count(myid) ? map_proc_IJs_avail.at(myid) : std::vector<atpair_t>();
    const auto ids_ap = get_2d_mat_indices_atpair(atbasis, atbasis, IJs, row_fast);

    // global basis indices required by BLACS
    const auto ids_blacs = get_2d_mat_indices_blacs(ad, myid, row_fast);

    std::unordered_map<int, atom_mapping<std::vector<std::pair<size_t, size_t>>>::pair_t> map_proc_map_id2d_ap;
    std::unordered_map<int, std::vector<atpair_t>> map_proc_IJs_save_order;
    for (const auto &id: ids_ap)
    {
        const auto &i_gl = as_int(id.first);
        const auto &j_gl = as_int(id.second);
        // Only sending either half
        if (ap_upper_half && i_gl > j_gl) continue;
        if (!ap_upper_half && i_gl < j_gl) continue;
        // filter out the required BLACS indices that already owns by myid
        myprow = ScalapackConnector::indxg2p(i_gl, ad.mb(), dummy, ad.irsrc(), ad.nprows());
        mypcol = ScalapackConnector::indxg2p(j_gl, ad.nb(), dummy, ad.icsrc(), ad.npcols());
        myprow_T = ScalapackConnector::indxg2p(j_gl, ad.mb(), dummy, ad.irsrc(), ad.nprows());
        mypcol_T = ScalapackConnector::indxg2p(i_gl, ad.nb(), dummy, ad.icsrc(), ad.npcols());

        const auto target = Cblacs_pnum(ad.ictxt(), myprow, mypcol);
        const auto target_T = Cblacs_pnum(ad.ictxt(), myprow_T, mypcol_T);

        // if (myid == 0) std::cout << "id/target " << id << " " << target << std::endl;
        if (myid != target)
        {
            if (return_local)
            {
                int I, J;
                int i_lo, j_lo;
                atbasis.get_local_index(i_gl, I, i_lo);
                atbasis.get_local_index(j_gl, J, j_lo);
                auto &IJs_save_order = map_proc_IJs_save_order[target];
                const atpair_t IJ = {static_cast<atom_t>(I), static_cast<atom_t>(J)};
                if (std::find(IJs_save_order.cbegin(), IJs_save_order.cend(), IJ) == IJs_save_order.cend())
                {
                    IJs_save_order.emplace_back(IJ);
                }
                const atpair_t ij_lo = {static_cast<atom_t>(i_lo), static_cast<atom_t>(j_lo)};
                map_proc_map_id2d_ap[target][IJ].emplace_back(ij_lo);
            }
            else
            {
                // NOTE: global index requested, save all to the (0,0) dummy pair index
                const atpair_t IJ{0, 0};
                auto &IJs_save_order = map_proc_IJs_save_order[target];
                if (std::find(IJs_save_order.cbegin(), IJs_save_order.cend(), IJ) == IJs_save_order.cend())
                {
                    IJs_save_order.emplace_back(IJ);
                }
                map_proc_map_id2d_ap[target][IJ].emplace_back(id);
            }
        }

        // id_aps are only in either half, should complete the whole one except for the diagonals
        if (i_gl != j_gl && myid != target_T)
        {
            if (return_local)
            {
                int I, J;
                int i_lo, j_lo;
                // Sending the same element, but may be conjugate on the receiver side
                atbasis.get_local_index(i_gl, I, i_lo);
                atbasis.get_local_index(j_gl, J, j_lo);
                auto &IJs_save_order = map_proc_IJs_save_order[target_T];
                const atpair_t IJ = {static_cast<atom_t>(I), static_cast<atom_t>(J)};
                if (std::find(IJs_save_order.cbegin(), IJs_save_order.cend(), IJ) == IJs_save_order.cend())
                {
                    IJs_save_order.emplace_back(IJ);
                }
                const atpair_t ij_lo = {static_cast<atom_t>(i_lo), static_cast<atom_t>(j_lo)};
                map_proc_map_id2d_ap[target_T][IJ].emplace_back(ij_lo);
            }
            else
            {
                // NOTE: global index requested, save all to the (0,0) dummy pair index
                const atpair_t IJ{0, 0};
                auto &IJs_save_order = map_proc_IJs_save_order[target_T];
                if (std::find(IJs_save_order.cbegin(), IJs_save_order.cend(), IJ) == IJs_save_order.cend())
                {
                    IJs_save_order.emplace_back(IJ);
                }
                map_proc_map_id2d_ap[target_T][IJ].emplace_back(id);
            }

        }
    }

    for (const auto &target_map: map_proc_map_id2d_ap)
    {
        const auto &target = target_map.first;
        const auto &IJs_save_order = map_proc_IJs_save_order.at(target);
        for (const auto &IJ: IJs_save_order)
        {
            map_proc_id2d_ap[target].push_back({IJ, target_map.second.at(IJ)});
        }
    }

    // reverse map of map_proc_IJs_avail
    atom_mapping<int>::pair_t map_IJs_proc;
    for (const auto &proc_IJs: map_proc_IJs_avail)
    {
        const auto &proc = proc_IJs.first;
        for (const auto &IJ: proc_IJs.second)
        {
            map_IJs_proc[IJ] = proc;
        }
    }

    Cblacs_pcoord(ad.ictxt(), myid, &myprow, &mypcol);
    for (const auto &id: ids_blacs)
    {
        const bool id_in_upper = id.first <= id.second;
        const auto I = static_cast<atom_t>(atbasis.get_i_atom(id.first));
        const auto J = static_cast<atom_t>(atbasis.get_i_atom(id.second));
        const atpair_t IJ({I, J});

        // transpose case
        const atpair_t id_T{id.second, id.first};
        const auto I_T = static_cast<atom_t>(atbasis.get_i_atom(id_T.first));
        const auto J_T = static_cast<atom_t>(atbasis.get_i_atom(id_T.second));
        const atpair_t IJ_T({I_T, J_T});

        if ((ap_upper_half && id_in_upper) || (!ap_upper_half && !id_in_upper))
        {
            if (map_IJs_proc.count(IJ))
            {
                const auto &source = map_IJs_proc[IJ];
                // exclude self
                if (!include_self && source == myid) continue;
                if (return_local)
                {
                    // FIXME: now only works when myid == ad.myid
                    const auto i_lo = as_size(ad.indx_g2l_r(as_int(id.first)));
                    const auto j_lo = as_size(ad.indx_g2l_c(as_int(id.second)));
                    map_proc_id2d_blacs[source].push_back({i_lo, j_lo});
                }
                else
                {
                    map_proc_id2d_blacs[source].emplace_back(id);
                }
                map_lconj[source].push_back('n');
            }
        }
        else
        {
            // id not available in atom-pair blocks, must transpose before looking for source
            if (map_IJs_proc.count(IJ_T))
            {
                const auto &source = map_IJs_proc[IJ_T];
                // exclude self
                if (!include_self && source == myid) continue;
                if (return_local)
                {
                    // FIXME: now only works when myid == ad.myid
                    const auto i_lo = as_size(ad.indx_g2l_r(as_int(id.first)));
                    const auto j_lo = as_size(ad.indx_g2l_c(as_int(id.second)));
                    map_proc_id2d_blacs[source].push_back({i_lo, j_lo});
                }
                else
                {
                    map_proc_id2d_blacs[source].emplace_back(id);
                }
                map_lconj[source].push_back('c');
            }
        }
    }

    return {map_proc_id2d_ap, map_proc_id2d_blacs};
}

/* ================================================================================= */
/* For general matrix redistribution from atom-pair (AP) to BLACS block-cyclic forms */
/* ================================================================================= */

std::pair<std::unordered_map<int, std::vector<size_t>>, std::unordered_map<int, std::vector<size_t>>>
get_communicate_global_ids_list_ap_to_blacs(const int &myid,
                                            const std::unordered_map<int, std::vector<atpair_t>> &map_proc_IJs_avail,
                                            const AtomicBasis &atbasis_r,
                                            const AtomicBasis &atbasis_c,
                                            const Array_Desc &ad,
                                            bool row_fast, bool row_major, bool include_self)
{
    // Objects to return
    std::unordered_map<int, std::vector<size_t>> map_proc_id1d_send;
    std::unordered_map<int, std::vector<size_t>> map_proc_id1d_recv;

    const auto ids = _get_communicate_ids_list_ap_and_blacs(myid, map_proc_IJs_avail, atbasis_r,
                                                            atbasis_c, ad, row_fast, false, include_self);

    // process ID -> list of basis indices, sending from myid to others
    std::unordered_map<int, std::vector<std::pair<size_t, size_t>>> map_proc_id2d_send;
    for (const auto &proc_id2d: ids.first)
    {
        const auto &proc = proc_id2d.first;
        // only up to 1 member for the case of global indices, which is for the the (0,0) dummy atom pair
        assert(proc_id2d.second.size() == 1);
        const auto &IJ_id2d = proc_id2d.second[0];
        assert(IJ_id2d.first.first == 0);
        assert(IJ_id2d.first.second == 0);
        map_proc_id2d_send[proc] = IJ_id2d.second;
    }

    // process ID -> list of basis indices, receiving by myid from others
    const auto &map_proc_id2d_recv = ids.second;

    // flatten 2D indices to 1D
    for (const auto &proc_id2d: map_proc_id2d_send)
    {
        map_proc_id1d_send[proc_id2d.first] = flatten_2d_indices(proc_id2d.second, atbasis_r.nb_total, atbasis_c.nb_total, row_major);
    }
    for (const auto &proc_id2d: map_proc_id2d_recv)
    {
        map_proc_id1d_recv[proc_id2d.first] = flatten_2d_indices(proc_id2d.second, atbasis_r.nb_total, atbasis_c.nb_total, row_major);
    }

    return {map_proc_id1d_send, map_proc_id1d_recv};
}

std::pair<std::unordered_map<int, std::vector<std::pair<atpair_t, std::vector<size_t>>>>,
          std::unordered_map<int, std::vector<size_t>>>
get_communicate_local_ids_list_ap_to_blacs(const int &myid,
                                           const std::unordered_map<int, std::vector<atpair_t>> &map_proc_IJs_avail,
                                           const AtomicBasis &atbasis_r,
                                           const AtomicBasis &atbasis_c,
                                           const Array_Desc &ad,
                                           bool row_fast, bool row_major, bool include_self)
{
    const auto ids = _get_communicate_ids_list_ap_and_blacs(myid, map_proc_IJs_avail, atbasis_r,
                                                            atbasis_c, ad, row_fast, true, include_self);
    // process ID -> list of basis indices, sending from myid to others
    const auto &map_proc_id2d_send = ids.first;
    // process ID -> list of basis indices, receiving by myid from others
    const auto &map_proc_id2d_recv = ids.second;

    // flatten 2D indices to 1D
    std::unordered_map<int, std::vector<std::pair<atpair_t, std::vector<size_t>>>> map_proc_id1d_send;
    std::unordered_map<int, std::vector<size_t>> map_proc_id1d_recv;

    for (const auto &proc_IJ_id2d: map_proc_id2d_send)
    {
        const auto &proc = proc_IJ_id2d.first;
        for (const auto &IJ_id2d: proc_IJ_id2d.second)
        {
            auto &proc_vec = map_proc_id1d_send[proc];
            const auto &IJ = IJ_id2d.first;
            const auto &I = IJ.first;
            const auto &J = IJ.second;
            const auto id1d = flatten_2d_indices(IJ_id2d.second, atbasis_r.get_atom_nb(I), atbasis_c.get_atom_nb(J), row_major);
            const std::pair<atpair_t, std::vector<size_t>> v(IJ, std::move(id1d));
            proc_vec.emplace_back(v);
        }
    }
    for (const auto &proc_id2d: map_proc_id2d_recv)
    {
        map_proc_id1d_recv[proc_id2d.first] = flatten_2d_indices(proc_id2d.second, as_size(ad.m_loc()), as_size(ad.n_loc()), row_major);
    }
    return {map_proc_id1d_send, map_proc_id1d_recv};
}

/* ============================================================================================= */
/* For symmetric/Hermitian matrix redistribution from atom-pair (AP) to BLACS block-cyclic forms */
/* ============================================================================================= */

std::pair<std::unordered_map<int, std::vector<size_t>>,
          std::unordered_map<int, std::pair<std::vector<size_t>, std::vector<char>>>>
get_communicate_global_ids_list_ap_to_blacs_sy(const int &myid,
                                               const char &uplo,
                                               const std::unordered_map<int, std::vector<atpair_t>> &map_proc_IJs_avail,
                                               const AtomicBasis &atbasis,
                                               const Array_Desc &ad,
                                               bool row_fast, bool row_major, bool include_self)
{
    // Objects to return
    std::unordered_map<int, std::vector<size_t>> map_proc_id1d_send;
    std::unordered_map<int, std::pair<std::vector<size_t>, std::vector<char>>> map_proc_id1d_lconj_recv;

    std::unordered_map<int, std::vector<char>> map_lconj;
    const auto ids = _get_communicate_ids_list_ap_to_blacs_sy(map_lconj, myid, uplo, map_proc_IJs_avail, atbasis,
                                                              ad, row_fast, false, include_self);

    // process ID -> list of basis indices, sending from myid to others
    std::unordered_map<int, std::vector<std::pair<size_t, size_t>>> map_proc_id2d_send;
    for (const auto &proc_id2d: ids.first)
    {
        const auto &proc = proc_id2d.first;
        // only up to 1 member for the case of global indices, which is for the the (0,0) dummy atom pair
        assert(proc_id2d.second.size() == 1);
        const auto &IJ_id2d = proc_id2d.second[0];
        assert(IJ_id2d.first.first == 0);
        assert(IJ_id2d.first.second == 0);
        map_proc_id2d_send[proc] = IJ_id2d.second;
    }

    // process ID -> list of basis indices, receiving by myid from others
    const auto &map_proc_id2d_recv = ids.second;

    // flatten 2D indices to 1D
    for (const auto &proc_id2d: map_proc_id2d_send)
    {
        map_proc_id1d_send[proc_id2d.first] = flatten_2d_indices(proc_id2d.second, atbasis.nb_total, atbasis.nb_total, row_major);
    }
    for (const auto &proc_data: map_proc_id2d_recv)
    {
        const auto &proc = proc_data.first;
        const auto &id2d = proc_data.second;
        // indices
        map_proc_id1d_lconj_recv[proc].first = flatten_2d_indices(id2d, atbasis.nb_total, atbasis.nb_total, row_major);
        // flag for conjugate
        map_proc_id1d_lconj_recv[proc].second = std::move(map_lconj[proc]);
    }

    return {map_proc_id1d_send, map_proc_id1d_lconj_recv};
}

std::pair<std::unordered_map<int, std::vector<std::pair<atpair_t, std::vector<size_t>>>>,
          std::unordered_map<int, std::pair<std::vector<size_t>, std::vector<char>>>>
get_communicate_local_ids_list_ap_to_blacs_sy(const int &myid,
                                              const char &uplo,
                                              const std::unordered_map<int, std::vector<atpair_t>> &map_proc_IJs_avail,
                                              const AtomicBasis &atbasis,
                                              const Array_Desc &ad,
                                              bool row_fast, bool row_major, bool include_self)
{
    std::unordered_map<int, std::vector<char>> map_lconj;
    const auto ids = _get_communicate_ids_list_ap_to_blacs_sy(map_lconj, myid, uplo, map_proc_IJs_avail, atbasis,
                                                              ad, row_fast, true, include_self);
    // process ID -> list of basis indices, sending from myid to others
    const auto &map_proc_id2d_send = ids.first;
    // process ID -> list of basis indices, receiving by myid from others
    const auto &map_proc_id2d_recv = ids.second;

    // flatten 2D indices to 1D
    std::unordered_map<int, std::vector<std::pair<atpair_t, std::vector<size_t>>>> map_proc_id1d_send;
    std::unordered_map<int, std::pair<std::vector<size_t>, std::vector<char>>> map_proc_id1d_lconj_recv;

    for (const auto &proc_IJ_id2d: map_proc_id2d_send)
    {
        const auto &proc = proc_IJ_id2d.first;
        for (const auto &IJ_id2d: proc_IJ_id2d.second)
        {
            auto &proc_vec = map_proc_id1d_send[proc];
            const auto &IJ = IJ_id2d.first;
            const auto &I = IJ.first;
            const auto &J = IJ.second;
            const auto id1d = flatten_2d_indices(IJ_id2d.second, atbasis.get_atom_nb(I), atbasis.get_atom_nb(J), row_major);
            const std::pair<atpair_t, std::vector<size_t>> v(IJ, std::move(id1d));
            proc_vec.emplace_back(v);
        }
    }
    for (const auto &proc_data: map_proc_id2d_recv)
    {
        const auto &proc = proc_data.first;
        const auto &id2d = proc_data.second;
        // indices
        map_proc_id1d_lconj_recv[proc].first = flatten_2d_indices(id2d, as_size(ad.m_loc()), as_size(ad.n_loc()), row_major);
        // flag for conjugate
        map_proc_id1d_lconj_recv[proc].second = std::move(map_lconj[proc]);
    }
    return {map_proc_id1d_send, map_proc_id1d_lconj_recv};
}

/* ================================================================================= */
/* For general matrix redistribution from BLACS block-cyclic to atom-pair (AP) forms */
/* ================================================================================= */

std::pair<std::unordered_map<int, std::vector<size_t>>, std::unordered_map<int, std::vector<size_t>>>
get_communicate_global_ids_list_blacs_to_ap(const int &myid,
                                            const std::unordered_map<int, std::vector<atpair_t>> &map_proc_IJs_require,
                                            const AtomicBasis &atbasis_r,
                                            const AtomicBasis &atbasis_c,
                                            const Array_Desc &ad,
                                            bool row_fast, bool row_major, bool include_self)
{
    // Objects to return
    std::unordered_map<int, std::vector<size_t>> map_proc_id1d_send;
    std::unordered_map<int, std::vector<size_t>> map_proc_id1d_recv;

    const auto ids = _get_communicate_ids_list_ap_and_blacs(myid, map_proc_IJs_require, atbasis_r,
                                                           atbasis_c, ad, row_fast, false, include_self);

    std::unordered_map<int, std::vector<std::pair<size_t, size_t>>> map_proc_id2d_recv;
    for (const auto &proc_id2d: ids.first)
    {
        const auto &proc = proc_id2d.first;
        // only up to 1 member for the case of global indices, which is for the the (0,0) dummy atom pair
        assert(proc_id2d.second.size() == 1);
        const auto &IJ_id2d = proc_id2d.second[0];
        assert(IJ_id2d.first.first == 0);
        assert(IJ_id2d.first.second == 0);
        map_proc_id2d_recv[proc] = IJ_id2d.second;
    }

    // process ID -> list of basis indices, receiving by myid from others
    const auto &map_proc_id2d_send = ids.second;

    // flatten 2D indices to 1D
    for (const auto &proc_id2d: map_proc_id2d_send)
    {
        map_proc_id1d_send[proc_id2d.first] = flatten_2d_indices(proc_id2d.second, atbasis_r.nb_total, atbasis_c.nb_total, row_major);
    }
    for (const auto &proc_id2d: map_proc_id2d_recv)
    {
        map_proc_id1d_recv[proc_id2d.first] = flatten_2d_indices(proc_id2d.second, atbasis_r.nb_total, atbasis_c.nb_total, row_major);
    }

    return {map_proc_id1d_send, map_proc_id1d_recv};
}

std::pair<std::unordered_map<int, std::vector<size_t>>,
          std::unordered_map<int, std::vector<std::pair<atpair_t, std::vector<size_t>>>>>
get_communicate_local_ids_list_blacs_to_ap(const int &myid,
                                           const std::unordered_map<int, std::vector<atpair_t>> &map_proc_IJs_require,
                                           const AtomicBasis &atbasis_r,
                                           const AtomicBasis &atbasis_c,
                                           const Array_Desc &ad,
                                           bool row_fast, bool row_major, bool include_self)
{
    const auto ids = _get_communicate_ids_list_ap_and_blacs(myid, map_proc_IJs_require, atbasis_r,
                                                            atbasis_c, ad, row_fast, true, include_self);
    const auto &map_proc_id2d_send = ids.second;
    const auto &map_proc_id2d_recv = ids.first;

    std::unordered_map<int, std::vector<size_t>> map_proc_id1d_send;
    std::unordered_map<int, std::vector<std::pair<atpair_t, std::vector<size_t>>>> map_proc_id1d_recv;

    for (const auto &proc_IJ_id2d: map_proc_id2d_recv)
    {
        const auto &proc = proc_IJ_id2d.first;
        for (const auto &IJ_id2d: proc_IJ_id2d.second)
        {
            auto &proc_vec = map_proc_id1d_recv[proc];
            const auto &IJ = IJ_id2d.first;
            const auto &I = IJ.first;
            const auto &J = IJ.second;
            const auto id1d = flatten_2d_indices(IJ_id2d.second, atbasis_r.get_atom_nb(I), atbasis_c.get_atom_nb(J), row_major);
            const std::pair<atpair_t, std::vector<size_t>> v(IJ, std::move(id1d));
            proc_vec.emplace_back(v);
        }
    }
    for (const auto &proc_id2d: map_proc_id2d_send)
    {
        map_proc_id1d_send[proc_id2d.first] = flatten_2d_indices(proc_id2d.second, as_size(ad.m_loc()), as_size(ad.n_loc()), row_major);
    }
    return {map_proc_id1d_send, map_proc_id1d_recv};
}

} /* end of namespace utils */

} /* end of namespace LIBRPA */
