// ====================
//  Author: Peize Lin
//  Date:   2022.08.11
// ====================

#pragma once

#include "Map_Operator-3.h"

namespace RI
{

namespace Map_Operator
{
	template<typename Tkey, typename Tvalue>
	std::map<Tkey,Tvalue> transform_prototype(
		const std::map<Tkey,Tvalue> &m_in,
		const std::function<Tvalue(const Tvalue&)> &func_prototype)
	{
		std::map<Tkey,Tvalue> m_out;
		for(const auto &item : m_in)
			m_out.emplace_hint(m_out.end(), item.first, func_prototype(item.second));
		return m_out;
	}

	template<typename Tkey, typename Tdata>
	std::map<Tkey,Tdata> transform(
		const std::map<Tkey,Tdata> &m_in,
		const std::function<Tdata(const Tdata&)> &func)
	{
		const std::function<Tdata(const Tdata&)>
			func_prototype = [&func](const Tdata& v) -> Tdata
			{ return func(v); };
		return transform_prototype(m_in, func_prototype);
	}

	template<typename Tkey1, typename Tkey2, typename Tvalue, typename Tdata>
	std::map<Tkey1,std::map<Tkey2,Tvalue>> transform(
		const std::map<Tkey1,std::map<Tkey2,Tvalue>> &m_in,
		const std::function<Tdata(const Tdata&)> &func)
	{
		using Tvalue1 = std::map<Tkey2,Tvalue>;
		const std::function<Tvalue1(const Tvalue1&)>
			func_prototype = [&func](const Tvalue1& v) -> Tvalue1
			{ return transform(v, func); };
		return transform_prototype(m_in, func_prototype);
	}

	template<typename Tkey, typename Tvalue>
	void for_each_prototype(
		std::map<Tkey,Tvalue> &m,
		const std::function<void(Tvalue&)> &func_prototype)
	{
		for(auto &item : m)
			func_prototype(item.second);
	}

	template<typename Tkey, typename Tdata>
	void for_each(
		std::map<Tkey,Tdata> &m,
		const std::function<void(Tdata&)> &func)
	{
		const std::function<void(Tdata&)>
			func_prototype = [&func](Tdata& v)
			{ return func(v); };
		for_each_prototype(m, func_prototype);
	}

	template<typename Tkey1, typename Tkey2, typename Tvalue, typename Tdata>
	void for_each(
		std::map<Tkey1,std::map<Tkey2,Tvalue>> &m,
		const std::function<void(Tdata&)> &func)
	{
		using Tvalue1 = std::map<Tkey2,Tvalue>;
		const std::function<void(Tvalue1&)>
			func_prototype = [&func](Tvalue1& v)
			{ return for_each(v, func); };
		for_each_prototype(m, func_prototype);
	}

	template<typename Tkey, typename Tdata>
	Tdata reduce(
		const std::map<Tkey,Tdata> &m,
		const Tdata &data_init,
		const std::function<Tdata(const Tdata&, const Tdata&)> &func)
	{
		Tdata data = data_init;
		for(const auto &item : m)
			data = func(item.second, data);
		return data;
	}

	template<typename TkeyA, typename TkeyB, typename Tvalue, typename Tdata>
	Tdata reduce(
		const std::map<TkeyA,std::map<TkeyB,Tvalue>> &m,
		const Tdata &data_init,
		const std::function<Tdata(const Tdata&, const Tdata&)> &func)
	{
		Tdata data = data_init;
		for(const auto &item : m)
			data = reduce(item.second, data, func);
		return data;
	}
}

}

#include "Map_Operator-2.hpp"