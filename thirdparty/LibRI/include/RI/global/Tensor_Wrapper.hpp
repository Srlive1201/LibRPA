// ===================
//  Author: Peize Lin
//  date: 2022.12.09
// ===================

#pragma once

#include "Tensor_Wrapper.h"

namespace RI
{
	template<typename T>
	Global_Func::To_Real_t<T> Tensor_Wrapper<T>::norm(const double p) const
	{
		using T_res = Global_Func::To_Real_t<T>;
		const std::size_t shape_all = this->shape.get_shape_all();
		if(p==2)
			return Blas_Interface::nrm2(*this);
		else if(p==1)
		{
			T_res s = 0;
			for(std::size_t i=0; i<shape_all; ++i)
				s += std::abs(this->ptr_[i]);
			return s;
		}
		else if(p==std::numeric_limits<double>::max())
		{
			T_res s = 0;
			for(std::size_t i=0; i<shape_all; ++i)
				s = std::max(std::real(s), std::abs(this->ptr_[i]));
			return s;
		}
		else
		{
			T_res s = 0;
			for(std::size_t i=0; i<shape_all; ++i)
				s += std::pow(std::abs(this->ptr_[i]), p);
			return std::pow(s,1.0/p);
		}
	}
}