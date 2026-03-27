// ====================
//  Author: Peize Lin
//  Date:   2022.07.27
// ====================

#pragma once

#include <map>
#include <functional>

namespace RI
{

namespace Map_Operator
{
	template<typename Tkey, typename Tvalue, typename Tdata>
	extern std::map<Tkey,Tvalue> zip_union(
		const std::map<Tkey,Tvalue> &m1,
		const std::map<Tkey,Tvalue> &m2,
		const std::function<Tdata(const Tdata&,const Tdata&)> &func);

	template<typename Tkey,
		typename Tvalue1, typename Tvalue2, typename Tvalue_return,
		typename Tdata1, typename Tdata2, typename Tdata_return>
	extern std::map<Tkey,Tvalue_return> zip_intersection(
		const std::map<Tkey,Tvalue1> &m1,
		const std::map<Tkey,Tvalue2> &m2,
		const std::function<Tdata_return(const Tdata1&,const Tdata2&)> &func);
}

}

#include "Map_Operator-2.hpp"