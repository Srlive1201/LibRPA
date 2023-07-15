#pragma once
#include <utility>
#include <ostream>
#include <set>
#include <array>
#include "atoms.h"
#include "vector3_order.h"
#ifdef LIBRPA_USE_LIBRI
#include <RI/global/Tensor.h>
#endif

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

#ifdef LIBRPA_USE_LIBRI
template <typename T>
std::ostream& operator<<(std::ostream& os, const RI::Tensor<T>& t)
{
    switch (t.shape.size())
    {
        case 1:
        {
            for (size_t i0 = 0; i0 < t.shape[0]; ++i0) os << t(i0) << " ";
            os << std::endl;
            return os;
        }
        case 2:
        {
            for (size_t i0 = 0; i0 < t.shape[0]; ++i0)
            {
                for (size_t i1 = 0; i1 < t.shape[1]; ++i1)
                    os << (std::abs(t(i0, i1)) > 1E-10 ? t(i0, i1) : 0) << " ";
                os << std::endl;
            }
            return os;
        }
        case 3:
        {
            os << "[" << std::endl;
            for (size_t i0 = 0; i0 < t.shape[0]; ++i0)
            {
                for (size_t i1 = 0; i1 < t.shape[1]; ++i1)
                {
                    for (size_t i2 = 0; i2 < t.shape[2]; ++i2) os << t(i0, i1, i2) << " ";
                    os << std::endl;
                }
                os << std::endl;
            }
            os << "]" << std::endl;
            return os;
        }
        default:
            throw std::invalid_argument(std::string(__FILE__) + " line " +
                                        std::to_string(__LINE__));
    }
}
#endif
