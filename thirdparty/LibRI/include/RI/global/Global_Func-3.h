// ===================
//  Author: Peize Lin
//  date: 2022.07.13
// ===================

#pragma once

#include "Array_Operator.h"
#include <vector>
#include <array>

namespace RI
{

namespace Global_Func
{
	template<typename Tcell, std::size_t Ndim>
	std::vector<std::array<Tcell, Ndim>> mod_period(
		const std::vector<std::array<Tcell, Ndim>> &cells_origin,
		const std::array<Tcell,Ndim> &period)
	{
		using namespace Array_Operator;
		std::vector<std::array<Tcell, Ndim>> cells_mod(cells_origin.size());
		for(std::size_t i=0; i<cells_origin.size(); ++i)
			cells_mod[i] = cells_origin[i] % period;
		return cells_mod;
	}
}

}