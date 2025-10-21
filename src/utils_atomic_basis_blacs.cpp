#include "utils_atomic_basis_blacs.h"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <map>
#include <stdexcept>

#include "base_utility.h"
#include "scalapack_connector.h"

#include "stl_io_helper.h"

namespace LIBRPA
{

namespace utils
{

std::set<std::pair<int, int>> get_necessary_IJ_from_block_2D(const AtomicBasis &atbasis_row, const AtomicBasis &atbasis_col, const Array_Desc& arrdesc)
{
    std::set<std::pair<int, int>> IJs;
    if (arrdesc.m() != atbasis_row.nb_total)
        throw std::invalid_argument("row basis and array descriptor inconsistent");
    if (arrdesc.n() != atbasis_col.nb_total)
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
    if (arrdesc.m() != atbasis.nb_total)
        throw std::invalid_argument("row basis and array descriptor inconsistent");
    if (arrdesc.n() != atbasis.nb_total)
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

static 
std::pair<std::map<int, std::vector<std::pair<atpair_t, std::vector<std::pair<size_t, size_t>>>>>,
          std::map<int, std::vector<std::pair<size_t, size_t>>>>
_get_communicate_ids_list_ap_to_blacs(const int &myid,
                                      const std::map<int, std::vector<atpair_t>> &map_proc_IJs_avail,
                                      const AtomicBasis &atbasis_r,
                                      const AtomicBasis &atbasis_c,
                                      const Array_Desc &ad, const bool row_fast, const bool return_local)
{
    assert(atbasis_r.nb_total == ad.m());
    assert(atbasis_c.nb_total == ad.n());

    // Objects to return
    // process ID -> list of basis indices, sending from myid to others
    std::map<int, std::vector<std::pair<atpair_t, std::vector<std::pair<size_t, size_t>>>>> map_proc_id2d_send;
    // process ID -> list of basis indices, receiving by myid from others
    std::map<int, std::vector<std::pair<size_t, size_t>>> map_proc_id2d_recv;

    int myprow, mypcol;
    int dummy = 0;

    // basis indices availble at myid in the atom-pair distribution
    const std::vector<atpair_t> &IJs_avail =
        map_proc_IJs_avail.count(myid) ? map_proc_IJs_avail.at(myid) : std::vector<atpair_t>();
    // if (myid == 0) std::cout << "IJs_avail " << IJs_avail << std::endl;
    const auto ids_ap = get_2d_mat_indices_atpair(atbasis_r, atbasis_c, IJs_avail, row_fast);

    // global basis indices required by BLACS
    const auto ids_blacs = get_2d_mat_indices_blacs(ad, myid, row_fast);

    std::map<int, std::map<atpair_t, std::vector<std::pair<size_t, size_t>>>> map_proc_map_id2d_send;
    std::map<int, std::vector<atpair_t>> map_proc_IJs_save_order;
    for (const auto &id: ids_ap)
    {
        const auto &i_gl = static_cast<int>(id.first);
        const auto &j_gl = static_cast<int>(id.second);
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
                map_proc_map_id2d_send[target][IJ].emplace_back(ij_lo);
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
                map_proc_map_id2d_send[target][IJ].emplace_back(id);
            }
        }
    }
    for (const auto &target_map: map_proc_map_id2d_send)
    {
        const auto &target = target_map.first;
        const auto &IJs_save_order = map_proc_IJs_save_order.at(target);
        for (const auto &IJ: IJs_save_order)
        {
            map_proc_id2d_send[target].push_back({IJ, target_map.second.at(IJ)});
        }
    }

    // reverse map of map_proc_IJs_avail
    std::map<atpair_t, int> map_IJs_proc;
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
        const auto I = static_cast<atom_t>(atbasis_r.get_i_atom(id.first));
        const auto J = static_cast<atom_t>(atbasis_c.get_i_atom(id.second));
        const atpair_t IJ({I, J});
        if (map_IJs_proc.count(IJ))
        {
            const auto &source = map_IJs_proc[IJ];
            // exclude self
            if (source == myid) continue;
            if (return_local)
            {
                // FIXME: wrong index now, only works when myid == ad.myid
                const auto i_lo = static_cast<size_t>(ad.indx_g2l_r(static_cast<int>(id.first)));
                const auto j_lo = static_cast<size_t>(ad.indx_g2l_c(static_cast<int>(id.second)));
                map_proc_id2d_recv[source].push_back({i_lo, j_lo});
            }
            else
            {
                map_proc_id2d_recv[source].emplace_back(id);
            }
        }
    }

    return {map_proc_id2d_send, map_proc_id2d_recv};
}

std::pair<std::map<int, std::vector<size_t>>, std::map<int, std::vector<size_t>>>
get_communicate_global_ids_list_ap_to_blacs(const int &myid,
                                            const std::map<int, std::vector<atpair_t>> &map_proc_IJs_avail,
                                            const AtomicBasis &atbasis_r,
                                            const AtomicBasis &atbasis_c,
                                            const Array_Desc &ad, bool row_fast, bool row_major)
{
    // Objects to return
    std::map<int, std::vector<size_t>> map_proc_id1d_send;
    std::map<int, std::vector<size_t>> map_proc_id1d_recv;

    const auto ids = _get_communicate_ids_list_ap_to_blacs(myid, map_proc_IJs_avail, atbasis_r,
                                                           atbasis_c, ad, row_fast, false);

    // process ID -> list of basis indices, sending from myid to others
    std::map<int, std::vector<std::pair<size_t, size_t>>> map_proc_id2d_send;
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

std::pair<std::map<int, std::vector<std::pair<atpair_t, std::vector<size_t>>>>,
          std::map<int, std::vector<size_t>>>
get_communicate_local_ids_list_ap_to_blacs(const int &myid,
                                           const std::map<int, std::vector<atpair_t>> &map_proc_IJs_avail,
                                           const AtomicBasis &atbasis_r,
                                           const AtomicBasis &atbasis_c,
                                           const Array_Desc &ad, bool row_fast, bool row_major)
{
    const auto ids = _get_communicate_ids_list_ap_to_blacs(myid, map_proc_IJs_avail, atbasis_r,
                                                           atbasis_c, ad, row_fast, true);
    // process ID -> list of basis indices, sending from myid to others
    const auto &map_proc_id2d_send = ids.first;
    // process ID -> list of basis indices, receiving by myid from others
    const auto &map_proc_id2d_recv = ids.second;

    // flatten 2D indices to 1D
    std::map<int, std::vector<std::pair<atpair_t, std::vector<size_t>>>> map_proc_id1d_send;
    std::map<int, std::vector<size_t>> map_proc_id1d_recv;

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
        map_proc_id1d_recv[proc_id2d.first] = flatten_2d_indices(proc_id2d.second, static_cast<size_t>(ad.m_loc()), static_cast<size_t>(ad.n_loc()), row_major);
    }
    return {map_proc_id1d_send, map_proc_id1d_recv};
}

} /* end of namespace utils */

} /* end of namespace LIBRPA */
