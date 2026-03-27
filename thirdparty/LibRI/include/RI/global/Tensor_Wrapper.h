// ===================
//  Author: Peize Lin
//  date: 2022.12.09
// ===================

#pragma once

#include "Global_Func-2.h"
#include <vector>

// Attention: very dangerous

namespace RI
{

template<typename T>
class Tensor_Wrapper
{
public:

	std::vector<std::size_t> shape;
	T *ptr_ = nullptr;

	Tensor_Wrapper()=default;
	explicit inline Tensor_Wrapper (const std::vector<std::size_t> &shape_in, T*const ptr_in) :shape(shape_in), ptr_(ptr_in){}

	T* ptr()const{ return this->ptr_; }

	// ||d||_p = (|d_1|^p+|d_2|^p+...)^{1/p}
	// if(p==std::numeric_limits<double>::max())    ||d||_max = max_i |d_i|
	Global_Func::To_Real_t<T> norm(const double p) const;
};

}

#include "Tensor_Wrapper.hpp"