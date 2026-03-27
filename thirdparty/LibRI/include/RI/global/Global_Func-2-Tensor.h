// ===================
//  Author: Peize Lin
//  date: 2022.06.02
// ===================

#pragma once

#include "Global_Func-2.h"
#include "Tensor.h"

#include <complex>
#include <type_traits>
#include <omp.h>

namespace RI
{

namespace Global_Func
{
	// t = convert(t)
	template<
		typename Tout, typename Tin,
		typename std::enable_if<std::is_same<Tin,Tout>::value,int>::type =0>
	Tensor<Tout> convert_Tensor(const Tensor<Tin> &t)
	{ return t; }

	// t = convert(t)
	template<
		typename Tout, typename Tin,
		typename std::enable_if<!std::is_same<Tin,Tout>::value,int>::type =0>
	Tensor<Tout> convert_Tensor(const Tensor<Tin> &t_in)
	{
		Tensor<Tout> t_out(t_in.shape);
		const Tin*const ptr_in = t_in.ptr();
		Tout*const ptr_out = t_out.ptr();
		const std::size_t size = t_in.shape.get_shape_all();
		#pragma omp parallel for
		for(std::size_t i=0; i<size; ++i)
			ptr_out[i] = Global_Func::convert<Tout>(ptr_in[i]);
		return t_out;
	}
}

}