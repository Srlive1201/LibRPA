//=======================
// AUTHOR : Peize Lin
// DATE :   2023-08-31
//=======================

#pragma once

#include "../../global/Cereal_Types.h"

#include <tuple>

namespace RI
{

namespace Communicate_Map_Combine
{
	template<typename TA, typename TC,
		template<typename T0, typename T1> class T_Judge_Map2_x,
		template<typename T0, typename T1> class T_Judge_Map2_y>
	class Judge_Map2_Combine2
	{
		using TAC = std::pair<TA,TC>;
	public:
		bool judge(const std::tuple<TA,TAC> &key) const
		{
			for(const auto &judge_i : std::get<0>(this->judge_list))
				if(judge_i.judge(key))
					return true;
			for(const auto &judge_i : std::get<1>(this->judge_list))
				if(judge_i.judge(key))
					return true;
			return false;
		}
		std::tuple<
			std::vector<T_Judge_Map2_x<TA,TAC>>,
			std::vector<T_Judge_Map2_y<TA,TC>> > judge_list;
		template <class Archive> void serialize( Archive & ar ){ ar(judge_list); }
	};
}

}