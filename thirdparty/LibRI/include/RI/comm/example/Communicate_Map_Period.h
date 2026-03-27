//=======================
// AUTHOR : Peize Lin
// DATE :   2022-07-06
//=======================

#pragma once

#include "../../global/Array_Operator.h"

#include <map>
#include <tuple>
#include <set>
#include <cereal/archives/binary.hpp>

namespace RI
{

namespace Communicate_Map_Period
{
	template<typename TA>
	class Judge_Map2_First
	{
	public:
		template<typename TC>
		bool judge(const std::tuple<TA,std::pair<TA,TC>> &key) const
		{
			return (s0.find(std::get<0>(key))!=s0.end())
				&& (s1.find(std::get<1>(key).first)!=s1.end());
		}
		std::set<TA> s0;
		std::set<TA> s1;
		template <class Archive> void serialize( Archive & ar ){ ar(s0); ar(s1); }
	};

	template<typename TA>
	class Judge_Map3_First
	{
	public:
		template<typename TP1, typename TP2>
		bool judge(const std::tuple<TA,std::pair<TA,TP1>,std::pair<TA,TP2>> &key) const
		{
			return (s0.find(std::get<0>(key))!=s0.end())
				&& (s1.find(std::get<1>(key).first)!=s1.end())
				&& (s2.find(std::get<2>(key).first)!=s2.end());
		}
		std::set<TA> s0;
		std::set<TA> s1;
		std::set<TA> s2;
		template <class Archive> void serialize( Archive & ar ){ ar(s0); ar(s1); ar(s2); }
	};

	template<typename TA, typename TC>
	class Judge_Map2_Period
	{
		using TAC = std::pair<TA,TC>;
	public:
		bool judge(const std::tuple<TA,TAC> &key) const
		{
			using namespace Array_Operator;
			for(const TAC &i0 : s0)
			{
				if(i0.first != std::get<0>(key))	continue;
				for(const TAC &i1 : s1)
				{
					if(i1.first != std::get<1>(key).first)	continue;
					if((i1.second-i0.second)%period != std::get<1>(key).second)	continue;
					return true;
				}
			}
			return false;
		}
		TC period;
		std::set<TAC> s0;
		std::set<TAC> s1;
		template <class Archive> void serialize( Archive & ar ){ ar(s0); ar(s1); ar(period); }
	};

	template<typename TA, typename TC>
	class Judge_Map3_Period
	{
		using TAC = std::pair<TA,TC>;
	public:
		bool judge(const std::tuple<TA,TAC,TAC> &key) const
		{
			using namespace Array_Operator;
			for(const TAC &i0 : s0)
			{
				if(i0.first != std::get<0>(key))	continue;
				for(const TAC &i1 : s1)
				{
					if(i1.first != std::get<1>(key).first)	continue;
					if((i1.second-i0.second)%period != std::get<1>(key).second)	continue;
					for(const TAC &i2 : s2)
					{
						if(i2.first != std::get<2>(key).first)	continue;
						if((i2.second-i0.second)%period != std::get<2>(key).second)	continue;
						return true;
					}
				}
			}
			return false;
		}
		TC period;
		std::set<TAC> s0;
		std::set<TAC> s1;
		std::set<TAC> s2;
		template <class Archive> void serialize( Archive & ar ){ ar(s0); ar(s1); ar(s2); ar(period); }
	};
}

}