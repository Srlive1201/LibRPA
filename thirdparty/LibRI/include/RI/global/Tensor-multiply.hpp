// ===================
//  Author: Peize Lin
//  date: 2021.10.31
// ===================

#pragma once

#include "Tensor.h"

#include "Blas_Interface-Tensor.h"
#include "Global_Func-2.h"

#include <cassert>
#include <stdexcept>
#include <string>

namespace RI
{

template<typename T,
	typename std::enable_if<!Global_Func::is_complex<T>::value,int>::type =0>
static T vector_dot_vector_tmp (const Tensor<T> &t1, const Tensor<T> &t2)
{
	return Blas_Interface::dot(t1, t2);
}
template<typename T,
	typename std::enable_if< Global_Func::is_complex<T>::value,int>::type =0>
static T vector_dot_vector_tmp (const Tensor<T> &t1, const Tensor<T> &t2)
{
	throw std::invalid_argument("complex vector * complex vector is forbidden, for ambigous dotu or dotc.\n"
		+std::string(__FILE__)+" line "+std::to_string(__LINE__));
}

template<typename T>
Tensor<T> operator* (const Tensor<T> &t1, const Tensor<T> &t2)
{
	switch(t1.shape.size())
	{
		case 1:
		{
			switch(t2.shape.size())
			{
				case 1:
				{
					assert(t1.shape[0] == t2.shape[0]);
					Tensor<T> t({1});
					//t(0) = Blas_Interface::dot(t1, t2);
					t(0) = vector_dot_vector_tmp(t1, t2);
					return t;
				}
				case 2:
				{
					assert(t1.shape[0] == t2.shape[0]);
					return Blas_Interface::gemv('T', T(1), t2, t1);
				}
				default:
					throw std::invalid_argument("Tensor multiply:\t"+std::to_string(t1.shape.size())+" * "+std::to_string(t2.shape.size())+".\n"+std::string(__FILE__)+" line "+std::to_string(__LINE__));
			}
		}
		case 2:
		{
			switch(t2.shape.size())
			{
				case 1:
				{
					assert(t1.shape[1] == t2.shape[0]);
					return Blas_Interface::gemv('N', T(1), t1, t2);
				}
				case 2:
				{
					assert(t1.shape[1] == t2.shape[0]);
					return Blas_Interface::gemm('N', 'N', T(1), t1, t2);
				}
				default:
					throw std::invalid_argument("Tensor multiply:\t"+std::to_string(t1.shape.size())+" * "+std::to_string(t2.shape.size())+".\n"+std::string(__FILE__)+" line "+std::to_string(__LINE__));
			}
		}
		default:
			throw std::invalid_argument("Tensor multiply:\t"+std::to_string(t1.shape.size())+" * "+std::to_string(t2.shape.size())+".\n"+std::string(__FILE__)+" line "+std::to_string(__LINE__));
	}
}

}