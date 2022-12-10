#pragma once
#include <utility>
#include <set>
#include <array>
#include "atoms.h"
#include "vector3_order.h"

template <typename TA, typename Tcell>
struct libri_types
{
    typedef std::array<Tcell, 3> TC;
    typedef std::pair<TA, TC> TAC;
};

template <typename TA, typename Tcell, typename TAout = TA, typename Tcellout = Tcell>
std::pair<std::set<TAout>, std::set<typename libri_types<TAout, Tcellout>::TAC>>
get_s0_s1_for_comm_map2(const std::vector<std::pair<TA, TA>>& atpairs,
                        const std::vector<Vector3_Order<Tcell>>& cells)
{
    std::set<TAout> set_s0;
    std::set<typename libri_types<TAout, Tcellout>::TAC> set_s1;
    for (const auto& atpair: atpairs)
    {
        set_s0.insert(atpair.first);
        for (const auto& cell: cells)
        {
            typename libri_types<TA, Tcell>::TC c {cell.x, cell.y, cell.z};
            set_s1.insert({atpair.second, c});
        }
    }
    return {set_s0, set_s1};
}

template <typename TA, typename TAout = TA>
std::pair<std::set<TAout>, std::set<TAout>>
get_s0_s1_for_comm_map2_first(const std::vector<std::pair<TA, TA>>& atpairs)
{
    std::set<TAout> set_s0;
    std::set<TAout> set_s1;
    for (const auto& atpair: atpairs)
    {
        set_s0.insert(atpair.first);
        set_s1.insert(atpair.second);
    }
    return {set_s0, set_s1};
}

template <typename TA, typename TAout = TA>
std::pair<std::set<TAout>, std::set<TAout>>
get_s0_s1_for_comm_map2_first(const std::set<std::pair<TA, TA>>& atpairs)
{
    std::set<TAout> set_s0;
    std::set<TAout> set_s1;
    for (const auto& atpair: atpairs)
    {
        set_s0.insert(atpair.first);
        set_s1.insert(atpair.second);
    }
    return {set_s0, set_s1};
}
