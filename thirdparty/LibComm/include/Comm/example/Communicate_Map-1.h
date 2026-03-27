//=======================
// AUTHOR : Peize Lin
// DATE :   2022-01-05
//=======================

#pragma once

#include <map>
#include <tuple>
#include <functional>

namespace Comm
{

namespace Communicate_Map
{
	/*
	template<typename Tkey0, typename Tkey1, typename Tvalue>
	inline Tvalue &get_value(
		const std::tuple<Tkey0,Tkey1> &key)
	{
		return data[std::get<0>(key)][std::get<1>(key)];
	}

	Tvalue &get_value(
		const std::tuple<Tkey0,Tkey1,Tkey2> &key)
	{
		return data[std::get<0>(key)][std::get<1>(key)][std::get<2>(key)];
	}
	*/

	// 等于
	/*
	template<typename Tkey0, typename Tkey1, typename Tvalue>
	void set_value_assignment(
		const std::tuple<Tkey0,Tkey1> &key,
		const Tvalue &value,
		std::map<Tkey0,std::map<Tkey1,Tvalue>> &data)
	{
		data[std::get<0>(key)][std::get<1>(key)] = value;
	}
	*/

	template<typename Tkey, typename Tvalue>
	void set_value_assignment(
		Tkey &&key,
		Tvalue &&value,
		std::map<Tkey,Tvalue> &data)
	{
		data[key] = std::move(value);
	}

	template<typename Tkey0, typename Tkey1, typename Tvalue>
	void set_value_assignment(
		std::tuple<Tkey0,Tkey1> &&key,
		Tvalue &&value,
		std::map<Tkey0,std::map<Tkey1,Tvalue>> &data)
	{
		data[std::get<0>(key)][std::get<1>(key)] = std::move(value);
	}

	template<typename Tkey0, typename Tkey1, typename Tkey2, typename Tvalue>
	void set_value_assignment(
		std::tuple<Tkey0,Tkey1,Tkey2> &&key,
		Tvalue &&value,
		std::map<Tkey0,std::map<Tkey1,std::map<Tkey2,Tvalue>>> &data)
	{
		data[std::get<0>(key)][std::get<1>(key)][std::get<2>(key)] = std::move(value);
	}

	// 加
	/*
	template<typename Tkey0, typename Tkey1, typename Tvalue>
	void set_value_add(
		const std::tuple<Tkey0,Tkey1> &key,
		const Tvalue &value,
		std::map<Tkey0,std::map<Tkey1,Tvalue>> &data)
	{
		Tvalue &value_tmp = data[std::get<0>(key)][std::get<1>(key)];
		if (value_tmp.empty())
			value_tmp = value;
		else
			value_tmp = value_tmp + value;
	}
	*/

	template<typename Tkey, typename Tvalue>
	void set_value_add(
		Tkey &&key,
		Tvalue &&value,
		std::map<Tkey,Tvalue> &data)
	{
		auto ptr = data.find(key);
		if(ptr==data.end())
			data[key] = std::move(value);
		else
			ptr->second = ptr->second + std::move(value);
//		Tvalue &value_tmp = data[key];
//		if (!value_tmp)
//			value_tmp = std::move(value);
//		else
//			value_tmp = value_tmp + std::move(value);
	}

	template<typename Tkey0, typename Tkey1, typename Tvalue>
	void set_value_add(
		std::tuple<Tkey0,Tkey1> &&key,
		Tvalue &&value,
		std::map<Tkey0,std::map<Tkey1,Tvalue>> &data)
	{
		set_value_add(std::move(std::get<1>(key)), std::move(value), data[std::get<0>(key)]);
//		Tvalue &value_tmp = data[std::get<0>(key)][std::get<1>(key)];
//		if (!value_tmp)
//			value_tmp = std::move(value);
//		else
//			value_tmp = value_tmp + std::move(value);
	}

	template<typename Tkey0, typename Tkey1, typename Tkey2, typename Tvalue>
	void set_value_add(
		std::tuple<Tkey0,Tkey1,Tkey2> &&key,
		Tvalue &&value,
		std::map<Tkey0,std::map<Tkey1,std::map<Tkey2,Tvalue>>> &data)
	{
		set_value_add(std::move(std::get<2>(key)), std::move(value), data[std::get<0>(key)][std::get<1>(key)]);
//		Tvalue &value_tmp = data[std::get<0>(key)][std::get<1>(key)][std::get<2>(key)];
//		if (!value_tmp)
//			value_tmp = std::move(value);
//		else
//			value_tmp = value_tmp + std::move(value);
	}

	/*
	// 等于
	template<typename Tkey0, typename Tkey1, typename Tkey2, typename Tvalue>
	set_value(
		const std::tuple<Tkey0,Tkey1,Tkey2> &key,
		const Tvalue &value,
		std::map<Tkey0,std::map<Tkey1,std::map<Tkey2,Tvalue>>> &data)
	{
		data[std::get<0>(key)][std::get<1>(key)][std::get<2>(key)] = value;
	}

	// 加
	set_value(
		const std::tuple<Tkey0,Tkey1,Tkey2> &key,
		const Tvalue &value,
		std::map<Tkey0,std::map<Tkey1,std::map<Tkey2,Tvalue>>> &data)
	{
		Tensor<> &data_tmp = data[std::get<0>(key)][std::get<1>(key)][std::get<2>(key)];
		if (data_tmp.c)
			data_tmp = data_tmp + value;
		else
			data_tmp = value;
	}
	*/

	// 无筛选，全部遍历
	template<typename Tkey, typename Tvalue>
	void traverse_datas_all(
		const std::map<Tkey,Tvalue> &data,
		const int rank_isend,
		std::function<void(const Tkey&, const Tvalue&)> &func)
	{
		for (const auto &item : data)
			func (item.first, item.second);
	}

	template<typename Tkey0, typename Tkey1, typename Tvalue>
	void traverse_datas_all(
		const std::map<Tkey0,std::map<Tkey1,Tvalue>> &data,
		const int rank_isend,
		std::function<void(const std::tuple<Tkey0,Tkey1>&, const Tvalue&)> &func)
	{
		for (const auto &dataA : data)
			for (const auto &dataB : dataA.second)
				func (std::make_tuple(dataA.first,dataB.first), dataB.second);
	}



	// 加横纵接收筛选。进程给出的label不做筛选。即，当set_value为=时需所有进程给出的label重复很少，或set_value为+。
	template<typename Tkey0, typename Tkey1, typename Tvalue>
	void traverse_datas_mask(
		const std::map<Tkey0,std::map<Tkey1,Tvalue>> &data,
		const int rank_isend,
		std::function<void(const std::tuple<Tkey0,Tkey1>&, const Tvalue&)> &func,
		const std::function<bool(const int, const Tkey0&)> &mask0,
		const std::function<bool(const int, const Tkey1&)> &mask1)
	{
		for (const auto &dataA : data)
		{
			const Tkey0 &key0 = dataA.first;
			if (!mask0(rank_isend,key0))	continue;
			for (const auto &dataB : dataA.second)
			{
				const Tkey1 &key1 = dataB.first;
				if (!mask1(rank_isend,key1))	continue;
				func (std::make_tuple(key0,key1), dataB.second);
			}
		}
	}


	template<typename Tkey, typename Tvalue>
	std::map<Tkey,Tvalue> init_datas_local(const int rank_recv)
	{
		return std::map<Tkey,Tvalue>();
	}

	template<typename Tkey0, typename Tkey1, typename Tvalue>
	std::map<Tkey0,std::map<Tkey1,Tvalue>> init_datas_local(const int rank_recv)
	{
		return std::map<Tkey0,std::map<Tkey1,Tvalue>>();
	}

	template<typename Tkey0, typename Tkey1, typename Tkey2, typename Tvalue>
	std::map<Tkey0,std::map<Tkey1,std::map<Tkey2,Tvalue>>> init_datas_local(const int rank_recv)
	{
		return std::map<Tkey0,std::map<Tkey1,std::map<Tkey2,Tvalue>>>();
	}


	template<typename Tkey, typename Tvalue>
	void add_datas(
		std::map<Tkey,Tvalue> &&data_local,
		std::map<Tkey,Tvalue> &data_recv)
	{
		auto ptr_local=data_local.begin();
		auto ptr_recv=data_recv.begin();
		for(; ptr_local!=data_local.end() && ptr_recv!=data_recv.end(); )
		{
			const Tkey &key_local = ptr_local->first;
			const Tkey &key_recv = ptr_recv->first;
			if(key_local == key_recv)
			{
				ptr_recv->second = ptr_recv->second + std::move(ptr_local->second);
				++ptr_local;
				++ptr_recv;
			}
			else if(key_local < key_recv)
			{
				ptr_recv = data_recv.emplace_hint(ptr_recv, key_local, std::move(ptr_local->second));
				++ptr_local;
			}
			else
			{
				++ptr_recv;
			}
		}
		for(; ptr_local!=data_local.end(); ++ptr_local)
		{
			ptr_recv = data_recv.emplace_hint(ptr_recv, ptr_local->first, std::move(ptr_local->second));
		}

//		for (auto &&data_local_A : data_local)
//		{
//			auto ptr = data_recv.find(data_local_A.first);
//			if(ptr==data_recv.end())
//				data_recv[data_local_A.first] = std::move(data_local_A.second);
//			else
//				ptr->second = ptr->second + std::move(data_local_A.second);
//		}

//		for (auto &&data_local_A : data_local)
//		{
//			Tvalue &data = data_recv[std::move(data_local_A.first)];
//			if(!data)
//				data = std::move(data_local_A.second);
//			else
//				data = data + std::move(data_local_A.second);
//		}
	}

	template<typename Tkey0, typename Tkey1, typename Tvalue>
	void add_datas(
		std::map<Tkey0,std::map<Tkey1,Tvalue>> &&data_local,
		std::map<Tkey0,std::map<Tkey1,Tvalue>> &data_recv)
	{
		auto ptr_local=data_local.begin();
		auto ptr_recv=data_recv.begin();
		for(; ptr_local!=data_local.end() && ptr_recv!=data_recv.end(); )
		{
			const Tkey0 &key_local = ptr_local->first;
			const Tkey0 &key_recv = ptr_recv->first;
			if(key_local == key_recv)
			{
				add_datas(std::move(ptr_local->second), ptr_recv->second);
				++ptr_local;
				++ptr_recv;
			}
			else if(key_local < key_recv)
			{
				ptr_recv = data_recv.emplace_hint(ptr_recv, key_local, std::move(ptr_local->second));
				++ptr_local;
			}
			else
			{
				++ptr_recv;
			}
		}
		for(; ptr_local!=data_local.end(); ++ptr_local)
		{
			ptr_recv = data_recv.emplace_hint(ptr_recv, ptr_local->first, std::move(ptr_local->second));
		}

//		for (auto &&data_local_A : data_local)
//		{
//			for (auto &&data_local_B : data_local_A.second)
//			{
//				Tvalue &data = data_recv[std::move(data_local_A.first)][std::move(data_local_B.first)];
//				if(!data)
//					data = std::move(data_local_B.second);
//				else
//					data = data + std::move(data_local_B.second);
//			}
//		}
	}

//	template<typename Tkey0, typename Tkey1, typename Tkey2, typename Tvalue>
//	void add_datas(
//		std::map<Tkey0,std::map<Tkey1,std::map<Tkey2,Tvalue>>> &&data_local,
//		std::map<Tkey0,std::map<Tkey1,std::map<Tkey2,Tvalue>>> &data_recv)
//	{
//		for (auto &&data_local_A : data_local)
//		{
//			for (auto &&data_local_B : data_local_A.second)
//			{
//				for (auto &&data_local_C : data_local_B.second)
//				{
//					Tvalue &data = data_recv[std::move(data_local_A.first)][std::move(data_local_B.first)][std::move(data_local_C.first)];
//					if(!data)
//						data = std::move(data_local_C.second);
//					else
//						data = data + std::move(data_local_C.second);
//				}
//			}
//		}
//	}
}

}