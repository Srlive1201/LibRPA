// ====================
//  Author: Peize Lin
//  Date:   2022.08.11
// ====================

#pragma once

#include <map>
#include <functional>

namespace RI
{

namespace Map_Operator
{
	template<typename Tkey, typename Tvalue, typename Tdata>
	extern std::map<Tkey,Tvalue> transform(
		const std::map<Tkey,Tvalue> &m_in,
		const std::function<Tdata(const Tdata&)> &func);

	template<typename Tkey, typename Tvalue, typename Tdata>
	extern void foreach(
		const std::map<Tkey,Tvalue> &m_in,
		const std::function<void(const Tdata&)> &func);

	template<typename Tkey, typename Tvalue, typename Tdata>
	extern Tdata reduce(
		const std::map<Tkey,Tvalue> &m,
		const Tdata &data_init,
		const std::function<Tdata(const Tdata&, const Tdata&)> &func);
}

}

#include "Map_Operator-3.hpp"