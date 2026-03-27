// ===================
//  Author: Peize Lin
//  Date:  2022.07.25
// ===================

#pragma once

#include <map>

namespace RI
{

namespace Map_Operator
{
	// m1+m2
	template<typename Tkey, typename Tvalue>
	extern std::map<Tkey,Tvalue> operator+(const std::map<Tkey,Tvalue> &m1, const std::map<Tkey,Tvalue> &m2);

	// m2 cover m1
	template<typename Tkey, typename Tvalue>
	extern std::map<Tkey,Tvalue> cover(const std::map<Tkey,Tvalue> &m1, const std::map<Tkey,Tvalue> &m2);

	// -m_in
	template<typename Tkey, typename Tvalue>
	extern std::map<Tkey,Tvalue> operator-(const std::map<Tkey,Tvalue> &m_in);

	// m1-m2
	template<typename Tkey, typename Tvalue>
	extern std::map<Tkey,Tvalue> operator-(const std::map<Tkey,Tvalue> &m1, const std::map<Tkey,Tvalue> &m2);
}

}

#include "Map_Operator.hpp"