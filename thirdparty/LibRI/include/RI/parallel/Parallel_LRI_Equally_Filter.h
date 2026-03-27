// ===================
//  Author: Peize Lin
//  date: 2023.03.07
// ===================

#pragma once

#include "Parallel_LRI_Equally.h"

namespace RI
{

template<typename TA, typename Tcell, std::size_t Ndim, typename Tdata>
class Parallel_LRI_Equally_Filter: public Parallel_LRI_Equally<TA,Tcell,Ndim,Tdata>
{
public:
	using TC = std::array<Tcell,Ndim>;
	using TAC = std::pair<TA,TC>;
	using Tatom_pos = std::array<double,Ndim>;		// tmp

	const std::vector<TAC>& get_list_Ab2 (const TA &Aa01, const TAC &Aa2, const TAC &Ab01) const override
	{
		return Global_Func::find(this->list_Ab2_filter,Ab01);
	}

	void filter_Ab2 (const std::map<TA, std::map<TAC, Tensor<Tdata>>> &Ds_b);

// private:
public:
	std::map<TAC,std::vector<TAC>> list_Ab2_filter;
};

}

#include "Parallel_LRI_Equally_Filter.hpp"