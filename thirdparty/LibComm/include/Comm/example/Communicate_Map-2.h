//=======================
// AUTHOR : Peize Lin
// DATE :   2022-07-06
//=======================

#pragma once

#include <map>
#include <tuple>
#include <set>
#include <cereal/archives/binary.hpp>

namespace Comm
{

namespace Communicate_Map
{
	template<typename Tkey, typename Tvalue>
	void traverse_keys(
		const std::map<Tkey,Tvalue> &datas,
		std::function<void(const Tkey&)> &func)
	{
		for(const auto &data : datas)
			func(data.first);
	}
	template<typename Tkey0, typename Tkey1, typename Tvalue>
	void traverse_keys(
		const std::map<Tkey0, std::map<Tkey1,Tvalue>> &datas,
		std::function<void(const std::tuple<Tkey0,Tkey1> &)> &func)
	{
		for(const auto &data0 : datas)
			for(const auto &data1 : data0.second)
				func(std::make_tuple(std::cref(data0.first), std::cref(data1.first)));
	}
	template<typename Tkey0, typename Tkey1, typename Tkey2, typename Tvalue>
	void traverse_keys(
		const std::map<Tkey0, std::map<Tkey1, std::map<Tkey2,Tvalue>>> &datas,
		std::function<void(const std::tuple<Tkey0,Tkey1,Tkey2> &)> &func)
	{
		for(const auto &data0 : datas)
			for(const auto &data1 : data0.second)
				for(const auto &data2 : data1.second)
					func(std::make_tuple(std::cref(data0.first), std::cref(data1.first), std::cref(data2.first)));
	}

	template<typename Tkey, typename Tvalue>
	const Tvalue &get_value(
		const Tkey &key,
		const std::map<Tkey, Tvalue> &m)
	{
		return m.at(key);
	}
	template<typename Tkey0, typename Tkey1, typename Tvalue>
	const Tvalue &get_value(
		const std::tuple<Tkey0, Tkey1> &key,
		const std::map<Tkey0, std::map<Tkey1, Tvalue>> &m)
	{
		return m.at(std::get<0>(key)).at(std::get<1>(key));
	}
	template<typename Tkey0, typename Tkey1, typename Tkey2, typename Tvalue>
	const Tvalue &get_value(
		const std::tuple<Tkey0, Tkey1, Tkey2> &key,
		const std::map<Tkey0, std::map<Tkey1, std::map<Tkey2, Tvalue>>> &m)
	{
		return m.at(std::get<0>(key)).at(std::get<1>(key)).at(std::get<2>(key));
	}

	template<typename Tkey>
	class Judge_Map
	{
	public:
		bool judge(const Tkey &key) const
		{
			return s.find(key)!=s.end();
		}
		std::set<Tkey> s;
		template <class Archive> void serialize( Archive & ar ){ ar(s); }
	};

	template<typename Tkey0, typename Tkey1>
	class Judge_Map2
	{
	public:
		bool judge(const std::tuple<Tkey0,Tkey1> &key) const
		{
			return (s0.find(std::get<0>(key))!=s0.end())
				&& (s1.find(std::get<1>(key))!=s1.end());
		}
		std::set<Tkey0> s0;
		std::set<Tkey1> s1;
		template <class Archive> void serialize( Archive & ar ){ ar(s0); ar(s1); }
	};

	template<typename Tkey0, typename Tkey1, typename Tkey2>
	class Judge_Map3
	{
	public:
		bool judge(const std::tuple<Tkey0,Tkey1,Tkey2> &key) const
		{
			return (s0.find(std::get<0>(key))!=s0.end())
				&& (s1.find(std::get<1>(key))!=s1.end())
				&& (s2.find(std::get<2>(key))!=s2.end());
		}
		std::set<Tkey0> s0;
		std::set<Tkey1> s1;
		std::set<Tkey2> s2;
		template <class Archive> void serialize( Archive & ar ){ ar(s0); ar(s1); ar(s2); }
	};
}

}