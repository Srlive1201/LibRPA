/*
 * @file      libri_stub.h
 * @brief     Stubbing classes and functions for LibRI
 * @author    Min-Ye Zhang
 * @date      2024-07-08
 */
#pragma once
#ifndef LIBRPA_USE_LIBRI
#include <vector>
#include <memory>
#include <valarray>
#include <mpi.h>

namespace RI
{

template <typename Tdata>
class Tensor
{
public:
    Tensor()
    {};

    Tensor(const std::vector<int> &dimension, std::shared_ptr<std::valarray<Tdata>> data_in)
    {};

    Tensor(const std::initializer_list<std::size_t> &dimension, std::shared_ptr<std::valarray<Tdata>> data_in)
    {};

    inline Tdata operator() (const std::size_t i0, const std::size_t i1) const { return static_cast<Tdata>(0); };

    void clear() {};
};


namespace Communicate_Tensors_Map_Judge
{

template <typename TA, typename TAC, typename Tdata>
std::map<TA, std::map<TAC, Tensor<Tdata>>> comm_map2_first(
    const MPI_Comm &mpi_comm, const std::map<TA, std::map<TAC, Tensor<Tdata>>> &Ds_in,
    const std::set<TA> &s0, const std::set<TA> &s1)
{
    std::map<TA, std::map<TAC, Tensor<Tdata>>> dummy;
    return dummy;
}

} /* end of namespace Communicate_Tensors_Map_Judge */

} /* end of namespace RI */
#endif
