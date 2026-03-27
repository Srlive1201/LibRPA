// ===================
//  Author: Peize Lin
//  date: 2022.07.01
// ===================

#pragma once

#include "RI_Tools.h"
#include "../global/Array_Operator.h"

namespace RI
{

namespace RI_Tools
{
	template<typename TA, typename TAC, typename Tdata>
	std::map<TA, std::map<TAC, Tensor<Tdata>>> filter(
		const std::map<TA, std::map<TAC, Tensor<Tdata>>> &Ds,
		const T_filter_func<Tdata> &filter_func,
		const Global_Func::To_Real_t<Tdata> &threshold)
	{
		std::map<TA, std::map<TAC, Tensor<Tdata>>> Ds_filter;
		for(const auto &Ds_tmp0 : Ds)
			for(const auto &Ds_tmp1 : Ds_tmp0.second)
				if(filter_func( Ds_tmp1.second, threshold ))
					Ds_filter[Ds_tmp0.first][Ds_tmp1.first] = Ds_tmp1.second;
		return Ds_filter;
	}
	template<typename TA, typename TAC, typename Tdata>
	std::map<TA, std::map<TAC, Tensor<Tdata>>> filter(
		const std::map<TA, std::map<TAC, Tensor<Tdata>>> &Ds,
		const std::vector<T_filter_func<Tdata>> &filter_func_list,
		const Global_Func::To_Real_t<Tdata> &threshold)
	{
		std::map<TA, std::map<TAC, Tensor<Tdata>>> Ds_filter;
		for(const auto &Ds_tmp0 : Ds)
			for(const auto &Ds_tmp1 : Ds_tmp0.second)
				for(const auto &filter_func : filter_func_list)
					if(filter_func( Ds_tmp1.second, threshold ))
					{
						Ds_filter[Ds_tmp0.first][Ds_tmp1.first] = Ds_tmp1.second;
						break;
					}
		return Ds_filter;
	}


	template<typename TA, typename TC, typename Tdata>
	std::map<TA, std::map<std::pair<TA,TC>, Tensor<Tdata>>> cal_period(
		const std::map<TA, std::map<std::pair<TA,TC>, Tensor<Tdata>>> &Ds,
		const TC &period)
	{
		using namespace Array_Operator;
		std::map<TA, std::map<std::pair<TA,TC>, Tensor<Tdata>>> Ds_period;
		for(const auto &Ds_tmp0 : Ds)
		{
			for(const auto &Ds_tmp1 : Ds_tmp0.second)
			{
				Tensor<Tdata> &D_period = Ds_period[Ds_tmp0.first][{Ds_tmp1.first.first, Ds_tmp1.first.second % period}];
				if(D_period.empty())
					D_period = Ds_tmp1.second;					// share memory
				else
					D_period = D_period + Ds_tmp1.second;		// new tensor
			}
		}
		return Ds_period;
	}

	template<typename TA, typename TAC, typename Tvalue>
	std::vector<std::set<TA>> get_index(const std::map<TA, std::map<TAC, Tvalue>> &Ds)
	{
		std::vector<std::set<TA>> index(1);
		std::set<TA> index_i;
		for(const auto &Ds_A : Ds)
			for(const auto &Ds_B : Ds_A.second)
				index_i.insert(Ds_B.first.first);
		index[0] = std::move(index_i);
		return index;
	}
}

}