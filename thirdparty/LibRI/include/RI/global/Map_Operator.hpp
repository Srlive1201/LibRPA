// ===================
//  Author: Peize Lin
//  Date:  2022.07.25
// ===================

#pragma once

#include "Map_Operator.h"

namespace RI
{

namespace Map_Operator
{
	// m1+m2
	template<typename Tkey, typename Tvalue>
	std::map<Tkey,Tvalue> operator+(const std::map<Tkey,Tvalue> &m1, const std::map<Tkey,Tvalue> &m2)
	{
		std::map<Tkey,Tvalue> m;
		auto ptr1 = m1.begin();
		auto ptr2 = m2.begin();
		while(ptr1!=m1.end() && ptr2!=m2.end())
		{
			if(ptr1->first == ptr2->first)
			{
				m.emplace_hint(m.end(), ptr1->first, ptr1->second + ptr2->second);
				++ptr1;
				++ptr2;
			}
			else if(ptr1->first < ptr2->first)
			{
				m.emplace_hint(m.end(), ptr1->first, ptr1->second);
				++ptr1;
			}
			else
			{
				m.emplace_hint(m.end(), ptr2->first, ptr2->second);
				++ptr2;
			}
		}
		m.insert(ptr1, m1.end());
		m.insert(ptr2, m2.end());
		return m;
	}

	// m2 cover m1
	template<typename Tkey, typename Tvalue>
	std::map<Tkey,Tvalue> cover(const std::map<Tkey,Tvalue> &m1, const std::map<Tkey,Tvalue> &m2)
	{
		std::map<Tkey,Tvalue> m = m1;
		for(const auto &item : m2)
			m[item.first] = item.second;
		return m;
	}

	// -m_in
	template<typename Tkey, typename Tvalue>
	std::map<Tkey,Tvalue> operator-(const std::map<Tkey,Tvalue> &m_in)
	{
		std::map<Tkey,Tvalue> m;
		for(const auto &item : m_in)
			m.emplace_hint(m.end(), item.first, -item.second);
		return m;
	}

	// m1-m2
	template<typename Tkey, typename Tvalue>
	std::map<Tkey,Tvalue> operator-(const std::map<Tkey,Tvalue> &m1, const std::map<Tkey,Tvalue> &m2)
	{
		std::map<Tkey,Tvalue> m;
		auto ptr1 = m1.begin();
		auto ptr2 = m2.begin();
		while(ptr1!=m1.end() && ptr2!=m2.end())
		{
			if(ptr1->first == ptr2->first)
			{
				m.emplace_hint(m.end(), ptr1->first, ptr1->second - ptr2->second);
				++ptr1;
				++ptr2;
			}
			else if(ptr1->first < ptr2->first)
			{
				m.emplace_hint(m.end(), ptr1->first, ptr1->second);
				++ptr1;
			}
			else
			{
				m.emplace_hint(m.end(), ptr2->first, -ptr2->second);
				++ptr2;
			}
		}
		m.insert(ptr1, m1.end());
		while(ptr2!=m2.end())
		{
			m.emplace_hint(m.end(), ptr2->first, -ptr2->second);
			++ptr2;
		}
		return m;
	}
}

}