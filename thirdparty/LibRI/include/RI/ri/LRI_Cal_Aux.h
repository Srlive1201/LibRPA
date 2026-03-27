// ===================
//  Author: Peize Lin
//  date: 2022.08.12
// ===================

#pragma once

#include "../global/Tensor.h"
#include "Label.h"
#include "../global/Map_Operator.h"
#include "../global/Array_Operator.h"
#include "../parallel/Parallel_LRI.h"

#include <map>
#include <memory.h>
#include <cassert>
#include <stdexcept>
#include <omp.h>

namespace RI
{

namespace LRI_Cal_Aux
{
	template<typename Tdata>
	inline Tensor<Tdata> tensor3_merge(const Tensor<Tdata> &D, const bool flag_01_2)
	{
		assert(D.shape.size()==3);
		if(flag_01_2)
			return Tensor<Tdata>( {D.shape[0]*D.shape[1], D.shape[2]}, D.data );
		else
			return Tensor<Tdata>( {D.shape[0], D.shape[1]*D.shape[2]}, D.data );
	}

	template<typename Tdata>
	Tensor<Tdata> tensor3_transpose(const Tensor<Tdata> &D)
	{
		assert(D.shape.size()==3);
		Tensor<Tdata> D_new({D.shape[1], D.shape[0], D.shape[2]});
		for(std::size_t i0=0; i0<D.shape[0]; ++i0)
			for(std::size_t i1=0; i1<D.shape[1]; ++i1)
			{
				memcpy(
					D_new.ptr()+(i1*D.shape[0]+i0)*D.shape[2],
					D.ptr()+(i0*D.shape[1]+i1)*D.shape[2],
					D.shape[2]*sizeof(Tdata));
			}
		return D_new;
	}

	template<typename Tdata>
	inline void add_D(const Tensor<Tdata> &D_add, Tensor<Tdata> &D_result)
	{
		if(D_result.empty())
			D_result = D_add;
		else
			D_result = D_result + D_add;
	}

	template<typename Tdata>
	inline void add_Ds(
		Tensor<Tdata> &&D_add,
		Tensor<Tdata> &D_result,
		const double fac = 1.0)
	{
		if(D_result.empty())
		{
			if(1.0==fac)
				D_result = std::move(D_add);
			else if(-1.0==fac)
				D_result = -D_add;
			else
				D_result = Tdata(fac) * D_add;
		}
		else
		{
			if(1.0==fac)
				D_result += D_add;
			else if(-1.0==fac)
				D_result -= D_add;
			else
				D_result += Tdata(fac) * D_add;
		}
	}

	template<typename Tkey, typename Tvalue>
	void add_Ds(
		std::map<Tkey, Tvalue> &&Ds_add,
		std::map<Tkey, Tvalue> &Ds_result,
		const double fac = 1.0)
	{
		if(Ds_add.empty())
			return;
		if(Ds_result.empty() && 1.0==fac)
			Ds_result = std::move(Ds_add);
		else
		{
			for(auto &&Ds_add_A : Ds_add)
				add_Ds(std::move(Ds_add_A.second), Ds_result[Ds_add_A.first], fac);
			Ds_add.clear();
		}
	}

	/*
	template<typename Tvalue>
	void add_Ds(
		std::vector<Tvalue> &&Ds_add,
		std::vector<Tvalue> &Ds_result,
		const double fac = 1.0)
	{
		if(Ds_add.empty())
			return;
		if(Ds_result.empty() && 1.0==fac)
		{
			Ds_result = std::move(Ds_add);
			Ds_add.resize(Ds_result.size());
		}
		else
		{
			if(Ds_result.empty())
				Ds_result.resize(Ds_add.size());
			else
				assert(Ds_add.size()==Ds_result.size());
			for(std::size_t i=0; i<Ds_result.size(); ++i)
				add_Ds(std::move(Ds_add[i]), Ds_result[i], fac);
			Ds_add.clear();
			Ds_add.resize(Ds_result.size());
		}
	}
	*/

	template<typename TA, typename TAC, typename Tdata>
	void add_Ds_omp_try(
		std::map<TA, std::map<TAC, Tensor<Tdata>>> &Ds_result_thread,
		std::map<TA, std::map<TAC, Tensor<Tdata>>> &Ds_result,
		omp_lock_t &lock_Ds_result_add,
		const double &fac)
	{
		if( !Ds_result_thread.empty() && omp_test_lock(&lock_Ds_result_add) )
		{
			LRI_Cal_Aux::add_Ds(std::move(Ds_result_thread), Ds_result, fac);
			omp_unset_lock(&lock_Ds_result_add);
			Ds_result_thread.clear();
		}
	}

	template<typename TA, typename TAC, typename Tdata>
	void add_Ds_omp_wait(
		std::map<TA, std::map<TAC, Tensor<Tdata>>> &&Ds_result_thread,
		std::map<TA, std::map<TAC, Tensor<Tdata>>> &Ds_result,
		omp_lock_t &lock_Ds_result_add,
		const double &fac)
	{
		if(!Ds_result_thread.empty())
		{
			omp_set_lock(&lock_Ds_result_add);
			LRI_Cal_Aux::add_Ds(std::move(Ds_result_thread), Ds_result, fac);
			omp_unset_lock(&lock_Ds_result_add);
			Ds_result_thread.clear();
		}
	}

	template<typename TA, typename TAC, typename Tdata>
	void add_Ds_omp_try_map(
		std::map<TA, std::map<TAC, Tensor<Tdata>>> &Ds_result_thread,
		std::map<TA, std::map<TAC, Tensor<Tdata>>> &Ds_result,
		std::map<TA, omp_lock_t> &lock_Ds_result_add_map,
		const double &fac)
	{
		for(auto ptr=Ds_result_thread.begin(); ptr!=Ds_result_thread.end(); )
		{
			const TA key = ptr->first;
			if(omp_test_lock(&lock_Ds_result_add_map.at(key)))
			{
				LRI_Cal_Aux::add_Ds(std::move(ptr->second), Ds_result.at(key), fac);
				omp_unset_lock(&lock_Ds_result_add_map.at(key));
				ptr = Ds_result_thread.erase(ptr);
			}
			else
			{
				++ptr;
			}
		}
	}

	template<typename TA, typename TAC, typename Tdata>
	void add_Ds_omp_wait_map(
		std::map<TA, std::map<TAC, Tensor<Tdata>>> &Ds_result_thread,
		std::map<TA, std::map<TAC, Tensor<Tdata>>> &Ds_result,
		std::map<TA, omp_lock_t> &lock_Ds_result_add_map,
		const double &fac)
	{
		if(Ds_result_thread.empty())
			return;
		while(true)
		{
			add_Ds_omp_try_map(Ds_result_thread, Ds_result, lock_Ds_result_add_map, fac);
			if(Ds_result_thread.empty())
				return;
			#pragma omp taskyield
		}
	}

	/*
	template<typename T>
	inline bool judge_Ds_empty(const std::vector<T> &Ds)
	{
		for(const T &D : Ds)
			if(!D.empty())
				return false;
		return true;
	}
	*/

	inline int judge_x(const Label::ab_ab &label)
	{
		switch(label)
		{
			case Label::ab_ab::a1b0_a2b2:	case Label::ab_ab::a1b2_a2b0:
			case Label::ab_ab::a0b0_a2b2:	case Label::ab_ab::a0b2_a2b0:
			case Label::ab_ab::a0b0_a1b2:	case Label::ab_ab::a0b2_a1b0:
				return 0;
			case Label::ab_ab::a1b1_a2b2:	case Label::ab_ab::a1b2_a2b1:
			case Label::ab_ab::a0b1_a2b2:	case Label::ab_ab::a0b2_a2b1:
			case Label::ab_ab::a0b1_a1b2:	case Label::ab_ab::a0b2_a1b1:
				return 1;
			default:
				return -1;
		}
	}

	inline Label::ab get_abx(const Label::ab_ab &label)
	{
		switch(label)
		{
			case Label::ab_ab::a0b0_a1b2:	case Label::ab_ab::a0b0_a2b2:	return Label::ab::a0b0;
			case Label::ab_ab::a0b1_a1b2:	case Label::ab_ab::a0b1_a2b2:	return Label::ab::a0b1;
			case Label::ab_ab::a0b2_a1b0:	case Label::ab_ab::a1b0_a2b2:	return Label::ab::a1b0;
			case Label::ab_ab::a0b2_a1b1:	case Label::ab_ab::a1b1_a2b2:	return Label::ab::a1b1;
			case Label::ab_ab::a0b2_a2b0:	case Label::ab_ab::a1b2_a2b0:	return Label::ab::a2b0;
			case Label::ab_ab::a0b2_a2b1:	case Label::ab_ab::a1b2_a2b1:	return Label::ab::a2b1;
			default:	throw std::invalid_argument("get_abx");
		}
	}

	template<typename Tkey0, typename Tkey1, typename Tdata>
	std::map<Tkey0, std::map<Tkey1, Tensor<Tdata>>>
	cal_Ds_transpose(const std::map<Tkey0, std::map<Tkey1, Tensor<Tdata>>> &Ds)
	{
		std::map<Tkey0, std::map<Tkey1, Tensor<Tdata>>> Ds_transpose;
		for(const auto &Ds0 : Ds)
			for(const auto &Ds1 : Ds0.second)
				Ds_transpose[Ds0.first][Ds1.first] = tensor3_transpose(Ds1.second);
		return Ds_transpose;
	}

	// D[{A1,C1}] + C0 -> D[{A1,C1-C0}]
	template<typename TA, typename TC, typename Tdata>
	auto Ds_translate(
		std::map<std::pair<TA,TC>,Tensor<Tdata>> &&Ds_in,
		const TC &translation, const TC &period)
	-> std::map<std::pair<TA,TC>,Tensor<Tdata>>
	{
		using TAC = std::pair<TA,TC>;
		using namespace Array_Operator;
		std::map<TAC,Tensor<Tdata>> Ds_out;
		for(auto &&D_in : Ds_in)
		{
			const TAC key = {D_in.first.first, (D_in.first.second-translation)%period};
			Ds_out[key] = std::move(D_in.second);
		}
		return Ds_out;
	}
	// D[{A0,C0}] + {A1,C1} -> D[A0][{A1,C1-C0}]
	template<typename TA, typename TC, typename Tdata>
	auto Ds_exchange(
		std::map<std::pair<TA,TC>,Tensor<Tdata>> &&Ds_in,
		const std::pair<TA,TC> &key1_origin, const TC &period)
	-> std::map<TA, std::map<std::pair<TA,TC>,Tensor<Tdata>>>
	{
		using TAC = std::pair<TA,TC>;
		using namespace Array_Operator;
		std::map<TA, std::map<std::pair<TA,TC>,Tensor<Tdata>>> Ds_out;
		for(auto &&D_in : Ds_in)
		{
			const TA key0 = D_in.first.first;
			const TAC key1 = {key1_origin.first, (key1_origin.second-D_in.first.second)%period};
			Ds_out[key0][key1] = std::move(D_in.second);
		}
		return Ds_out;
	}

	template<typename TA, typename Tvalue>
	std::vector<TA> filter_list_map(
		const std::vector<TA> &list_in,
		const std::map<TA, Tvalue> &Ds)
	{
		std::vector<TA> list_filter;
		for(const TA &item : list_in)
			if(Ds.find(item) != Ds.end())
				list_filter.push_back(item);
		return list_filter;
	}

	template<typename TA, typename TC, typename Tvalue>
	std::vector<std::pair<TA,TC>> filter_list_map(
		const std::vector<std::pair<TA,TC>> &list_in,
		const std::map<TA, Tvalue> &Ds)
	{
		std::vector<std::pair<TA,TC>> list_filter;
		for(const std::pair<TA,TC> &item : list_in)
			if(Ds.find(item.first) != Ds.end())
				list_filter.push_back(item);
		return list_filter;
	}

	template<typename TA>
	std::vector<TA> filter_list_set(
		const std::vector<TA> &list_in,
		const std::set<TA> &index)
	{
		std::vector<TA> list_filter;
		for(const TA &item : list_in)
			if(index.find(item) != index.end())
				list_filter.push_back(item);
		return list_filter;
	}

	template<typename TA, typename TC>
	std::vector<std::pair<TA,TC>> filter_list_set(
		const std::vector<std::pair<TA,TC>> &list_in,
		const std::set<TA> &index)
	{
		std::vector<std::pair<TA,TC>> list_filter;
		for(const std::pair<TA,TC> &item : list_in)
			if(index.find(item.first) != index.end())
				list_filter.push_back(item);
		return list_filter;
	}

	template<typename TA, typename TAC, typename Tdata>
	std::map<TA, omp_lock_t> init_lock_result(
		const std::vector<Label::ab_ab> &labels,
		const std::unordered_map<Label::Aab_Aab, List_A<TA,TAC>> &list_A,
		std::map<TA, std::map<TAC, Tensor<Tdata>>> &Ds_result)
	{
		std::map<TA, omp_lock_t> lock_Ds_result_add_map;
		for(const Label::ab_ab &label : labels)
		{
			switch(Label_Tools::to_Aab_Aab(label))
			{
				case Label::Aab_Aab::a01b01_a2b01:
				case Label::Aab_Aab::a01b01_a2b2:
				case Label::Aab_Aab::a01b2_a2b01:
					for(const TA &Aa01 : list_A.at(Label_Tools::to_Aab_Aab(label)).a01)
					{
						Ds_result[Aa01];
						lock_Ds_result_add_map[Aa01];
					}
					break;
				case Label::Aab_Aab::a01b01_a01b01:
				case Label::Aab_Aab::a01b01_a01b2:
					for(const TAC &Aa2 : list_A.at(Label_Tools::to_Aab_Aab(label)).a2)
					{
						Ds_result[Aa2.first];
						lock_Ds_result_add_map[Aa2.first];
					}
					break;
				default:
					throw std::invalid_argument(std::string(__FILE__)+" line "+std::to_string(__LINE__));
			}
		}
		for(auto &lock : lock_Ds_result_add_map)
			omp_init_lock(&lock.second);
		return lock_Ds_result_add_map;
	}

	template<typename TA, typename Tvalue>
	void destroy_lock_result(
		std::map<TA, omp_lock_t> &locks,
		std::map<TA, Tvalue> &Ds_result)
	{
		for(auto &lock : locks)
			omp_destroy_lock(&lock.second);
		for(auto ptr=Ds_result.begin(); ptr!=Ds_result.end(); )
		{
			if(ptr->second.empty())
				ptr = Ds_result.erase(ptr);
			else
				++ptr;
		}
	}
}

}