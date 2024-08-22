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

    void clear() {};
};

} /* end of namespace RI */
#endif
