// ===================
//  Author: Peize Lin
//  date: 2022.06.06
// ===================

#pragma once

#include "Label.h"
#include "../global/Global_Func-2.h"
#include "../global/Tensor.h"

#include <map>

namespace RI
{

namespace CS_Matrix_Tools
{
	enum class Uplimit_Type
	{
		square_two,			norm_two,
		square_three_0,		norm_three_0,
		square_three_1,		norm_three_1,
		square_three_2,		norm_three_2
	};

	extern inline std::pair<Label::ab, Label::ab> split_label(const Label::ab_ab &label);

	template<typename Tkey, typename Tvalue>
	extern auto cal_uplimit(
		const Uplimit_Type &uplimit_type,
		const std::map<Tkey,Tvalue> &Ds)
	-> std::map<Tkey, decltype(cal_uplimit(uplimit_type,Ds.begin()->second))>;

	template<typename Tdata>
	extern Global_Func::To_Real_t<Tdata> cal_uplimit(
		const Uplimit_Type &uplimit_type,
		const Tensor<Tdata> &D);
}

}

#include "CS_Matrix_Tools.hpp"