#pragma once
#include <utility>
#include <ostream>
#include <set>
#include <array>
#include "../core/atom.h"
#include "../math/vector3_order.h"
#ifdef LIBRPA_USE_LIBRI
#include <RI/global/Tensor.h>
#else
#include "../utils/libri_stub.h"
#endif

namespace librpa_int {

template <typename TA, typename Tcell>
struct libri_types
{
    typedef std::array<Tcell, 3> TC;
    typedef std::pair<TA, TC> TAC;
};

template <typename TA, typename Tcell, typename Tdata>
size_t get_tensor_map_bytes(std::map<TA, std::map<std::pair<TA, Tcell>, RI::Tensor<Tdata>>> tensor_map)
{
    size_t msize = 0;
    for (const auto &I_JR_mat: tensor_map)
    {
        for (const auto &JR_mat: I_JR_mat.second)
        {
            const auto &mat = JR_mat.second;
            msize += mat.get_shape_all();
        }
    }
    return msize * sizeof(Tdata);
}

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

// Only works for real type
template <typename T>
T absmax(const RI::Tensor<T>& t)
{
    T maxval = T(-1.0);
    for (int i = 0; i < t.get_shape_all(); i++)
    {
        maxval = max(std::abs(t.ptr()[i]), maxval);
    }
    return maxval;
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

}
