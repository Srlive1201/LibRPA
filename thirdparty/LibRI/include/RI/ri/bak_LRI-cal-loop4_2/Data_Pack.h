// ===================
//  Author: Peize Lin
//  date: 2023.08.25
// ===================

#pragma once

#include "../global/Tensor.h"
#include "CS_Matrix.h"
#include "../global/Global_Func-1.h"

#include <vector>
#include <map>
#include <set>

namespace RI
{

template<typename TA, typename TC, typename Tdata>
struct Data_Pack
{
	using TAC = std::pair<TA,TC>;
	using Tdata_real = Global_Func::To_Real_t<Tdata>;

	std::map<TA, std::map<TAC, Tensor<Tdata>>> Ds_ab;					// Ds_ab[A0][{A1,C1}]
	std::vector<std::set<TA>> index_Ds_ab;								// index_Ds_ab[0]=A1
	typename CS_Matrix<TA,TC,Tdata_real>::Uplimits csm_uplimits;
};


template<typename TA, typename TC, typename Tdata>
class Data_Pack_Wrapper
{
public:
	Data_Pack_Wrapper(
		std::map<std::string, Data_Pack<TA,TC,Tdata>> &data_pool_in,
		std::unordered_map<Label::ab, std::string> &data_ab_name_in)
			:data_pool(data_pool_in), data_ab_name(data_ab_name_in){}
	inline Data_Pack<TA,TC,Tdata> &operator()(const Label::ab &label) const
	{
		return this->data_pool.at( this->data_ab_name.at(label) );
	}

	std::map<std::string, Data_Pack<TA,TC,Tdata>> &data_pool;
	std::unordered_map<Label::ab, std::string> &data_ab_name;
};

}