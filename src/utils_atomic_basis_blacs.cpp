#include "utils_atomic_basis_blacs.h"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <map>
#include <stdexcept>

#include "base_utility.h"
#include "scalapack_connector.h"

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

std::pair<std::map<int, std::vector<size_t>>, std::map<int, std::vector<size_t>>>
get_communicate_ids_list_ap_to_blacs(const int &myid,
                                     const std::map<int, std::vector<atpair_t>> &map_proc_IJs_avail,
                                     const AtomicBasis &atbasis_r,
                                     const AtomicBasis &atbasis_c,
                                     const Array_Desc &ad, bool row_fast, bool row_major)
{
    assert(atbasis_r.nb_total == ad.m());
    assert(atbasis_c.nb_total == ad.n());

    int myprow, mypcol;
    int dummy = 0;

    // basis indices availble at myid in the atom-pair distribution
    const std::vector<atpair_t> &IJs =
        map_proc_IJs_avail.count(myid) ? map_proc_IJs_avail.at(myid) : std::vector<atpair_t>();
    const auto ids_ap = get_2d_mat_indices_atpair(atbasis_r, atbasis_c, IJs, row_fast);

    // basis indices required by BLACS
    const auto ids_blacs = get_2d_mat_indices_blacs(ad, myid, row_fast);

    // process ID -> list of basis indices, sending from myid to others
    std::map<int, std::vector<std::pair<size_t, size_t>>> map_proc_id2d_send;
    for (const auto &id: ids_ap)
    {
        const auto &i = static_cast<int>(id.first);
        const auto &j = static_cast<int>(id.second);
        // filter out the required BLACS indices that already owns by myid
        myprow = ScalapackConnector::indxg2p(i, ad.mb(), dummy, ad.irsrc(), ad.nprows());
        mypcol = ScalapackConnector::indxg2p(j, ad.nb(), dummy, ad.icsrc(), ad.npcols());
        const auto target = Cblacs_pnum(ad.ictxt(), myprow, mypcol);

        if (myid != target)
        {
            map_proc_id2d_send[target].emplace_back(id);
        }
    }

    // process ID -> list of basis indices, receiving by myid from others
    std::map<int, std::vector<std::pair<size_t, size_t>>> map_proc_id2d_recv;
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
    for (const auto &id: ids_blacs)
    {
        const auto I = static_cast<atom_t>(atbasis_r.get_i_atom(id.first));
        const auto J = static_cast<atom_t>(atbasis_c.get_i_atom(id.second));
        const atpair_t IJ({I, J});
        if (map_IJs_proc.count(IJ))
        {
            const auto &source = map_IJs_proc[IJ];
            map_proc_id2d_recv[source].emplace_back(id);
        }
    }

    // flatten 2D indices to 1D
    std::map<int, std::vector<size_t>> map_proc_id1d_send;
    std::map<int, std::vector<size_t>> map_proc_id1d_recv;
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

} /* end of namespace utils */

} /* end of namespace LIBRPA */
