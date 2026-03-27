// ===================
//  Author: Peize Lin
//  date: 2022.07.01
// ===================

#pragma once

#include "../global/Tensor.h"
#include "../global/Global_Func-2.h"
#include <map>

namespace RI
{

namespace RI_Tools
{
	template<typename Tdata>
	using T_filter_func =
		std::function<bool(
			const Tensor<Tdata> &D,
			const Global_Func::To_Real_t<Tdata> &threshold)>;

	template<typename TA, typename TAC, typename Tdata>
	extern std::map<TA, std::map<TAC, Tensor<Tdata>>> filter(
		const std::map<TA, std::map<TAC, Tensor<Tdata>>> &Ds,
		const T_filter_func<Tdata> &filter_func,
		const Global_Func::To_Real_t<Tdata> &threshold);
	template<typename TA, typename TAC, typename Tdata>
	extern std::map<TA, std::map<TAC, Tensor<Tdata>>> filter(
		const std::map<TA, std::map<TAC, Tensor<Tdata>>> &Ds,
		const std::vector<T_filter_func<Tdata>> &filter_func_list,
		const Global_Func::To_Real_t<Tdata> &threshold);

	template<typename TA, typename TC, typename Tdata>
	extern std::map<TA, std::map<std::pair<TA,TC>, Tensor<Tdata>>> cal_period(
		const std::map<TA, std::map<std::pair<TA,TC>, Tensor<Tdata>>> &Ds,
		const TC &period);
}

}

#include "RI_Tools.hpp"