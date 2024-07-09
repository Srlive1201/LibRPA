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
#include <stdexcept>

namespace RI
{

template <typename Tdata>
class Tensor
{
public:
    Tensor()
    { throw std::logic_error("stub Tensor constructor is called"); };

    Tensor(const std::vector<int> &dimension, std::shared_ptr<std::valarray<Tdata>> data_in)
    { throw std::logic_error("stub Tensor constructor is called"); };

    void clear() {};
};

} /* end of namespace RI */
#endif
