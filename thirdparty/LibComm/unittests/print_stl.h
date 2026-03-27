// ===================
//  Author: Peize Lin
//  date: 2021.08.21
// ===================

#pragma once

#include <array>
#include <vector>
#include <map>
#include <set>
#include <tuple>
#include <iostream>

template<typename T>
std::ostream &operator<<(std::ostream &os, const std::vector<T> &v)
{
	for(size_t i=0; i<v.size(); ++i)
		os<<v[i]<<"|\t";
	return os;
}

template<typename T>
std::ostream &operator<<(std::ostream &os, const std::set<T> &s)
{
	for(const T &i : s)
		os<<i<<"|\t";
	return os;
}

template<typename T, size_t N>
std::ostream &operator<<(std::ostream &os, const std::array<T,N> &v)
{
	for(size_t i=0; i<v.size(); ++i)
		os<<v[i]<<"\t";
	return os;
}

template<typename Tkey, typename Tvalue>
std::ostream &operator<<(std::ostream &os, const std::pair<Tkey,Tvalue> &p)
{
	os<<"{ "<<p.first<<", "<<p.second<<" }";
	return os;
}

template<typename Tkey, typename Tvalue>
std::ostream &operator<<(std::ostream &os, const std::map<Tkey,Tvalue> &m)
{
	for(const auto &i : m)
		os<<i.first<<"\t"<<i.second<<std::endl;
	return os;
}

template<typename T0, typename T1>
std::ostream &operator<<(std::ostream &os, const std::tuple<T0,T1> &t)
{
	os<<"[ "<<std::get<0>(t)<<", "<<std::get<1>(t)<<" ]";
	return os;
}