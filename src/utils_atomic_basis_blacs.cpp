#include "utils_atomic_basis_blacs.h"

#include <algorithm>
#include <stdexcept>
#include <map>
#include <iostream>

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
get_communicate_ids_list_ap_to_blacs(int myid,
                                     const std::vector<size_t> & ids_ap,
                                     const std::vector<size_t> & ids_blacs,
                                     const std::vector<int> &proc_ids)
{
    // process ID -> list of basis indices
    std::map<int, std::vector<size_t>> proc2idlist_send;
    std::map<int, std::vector<size_t>> proc2idlist_recv;

    // build receive list from desired indices
    // std::cout << myid << std::endl;
    for (const auto &i_blacs: ids_blacs)
    {
        const auto &proc = proc_ids[i_blacs];
        // if (myid == 0) std::cout << proc << " i_blacs " << i_blacs << std::endl;
        if (proc == myid) continue;
        // if (myid == 0) std::cout << proc << " i_blacs " << i_blacs << " recv accepted" << std::endl;
        proc2idlist_recv[proc].push_back(i_blacs);
    }

    // build send list from owned indices
    for (const auto &i_ap: ids_ap)
    {
        // if (myid == 0) std::cout << "i_ap " << i_ap << std::endl;
        const auto &proc = proc_ids[i_ap];
        if (proc == myid) continue;
        proc2idlist_send[proc].push_back(i_ap);
    }

    return {proc2idlist_send, proc2idlist_recv};
}

} /* end of namespace utils */

} /* end of namespace LIBRPA */
