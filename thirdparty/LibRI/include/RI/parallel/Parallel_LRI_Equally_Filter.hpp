// ===================
//  Author: Peize Lin
//  date: 2023.03.07
// ===================

#pragma once

#include "Parallel_LRI_Equally_Filter.h"

namespace RI
{

template<typename TA, typename Tcell, std::size_t Ndim, typename Tdata>
void Parallel_LRI_Equally_Filter<TA,Tcell,Ndim,Tdata>::filter_Ab2 (
	const std::map<TA, std::map<TAC, Tensor<Tdata>>> &Ds_b)
{
	using namespace Array_Operator;
	this->list_Ab2_filter.clear();
	for(const TAC &Ab01 : this->list_Ab01)
	{
		for(const TAC &Ab2 : this->list_Ab2)
		{
			const Tensor<Tdata> D_b = Global_Func::find(
				Ds_b,
				Ab01.first, TAC{Ab2.first, (Ab2.second-Ab01.second)%this->period});
			if(!D_b.empty())
				this->list_Ab2_filter[Ab01].push_back(Ab2);
		}
	}
}

}