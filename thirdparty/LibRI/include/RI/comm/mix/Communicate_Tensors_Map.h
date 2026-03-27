//=======================
// AUTHOR : Peize Lin
// DATE :   2022-07-20
//=======================

#pragma once

#include "../../global/Cereal_Types.h"
#include "../../global/Tensor.h"
#include <Comm/Comm_Assemble/Comm_Assemble.h>
#include <Comm/example/Communicate_Map-1.h>
#include <Comm/example/Communicate_Map-2.h>

#include <map>
#include <mpi.h>

namespace RI
{

namespace Communicate_Tensors_Map
{
	template<typename TA, typename Tdata, typename Tjudge>
	std::map<TA,Tensor<Tdata>>
	comm_map(
		const MPI_Comm &mpi_comm,
		const std::map<TA,Tensor<Tdata>> &Ds_in,
		const Tjudge &judge)
	{
		Comm::Comm_Assemble<
			TA,
			Tensor<Tdata>,
			std::map<TA,Tensor<Tdata>>,
			Tjudge,
			std::map<TA,Tensor<Tdata>>
		> com(mpi_comm);

		com.traverse_keys_provide = Comm::Communicate_Map::traverse_keys<TA,Tensor<Tdata>>;
		com.get_value_provide = Comm::Communicate_Map::get_value<TA,Tensor<Tdata>>;
		com.set_value_require = Comm::Communicate_Map::set_value_add<TA,Tensor<Tdata>>;
		com.flag_lock_set_value = Comm::Comm_Tools::Lock_Type::Copy_merge;
		com.init_datas_local = Comm::Communicate_Map::init_datas_local<TA,Tensor<Tdata>>;
		com.add_datas = Comm::Communicate_Map::add_datas<TA,Tensor<Tdata>>;

		std::map<TA,Tensor<Tdata>> Ds_out;
		com.communicate( Ds_in, judge, Ds_out );
		return Ds_out;
	}

	template<typename TA, typename TAC, typename Tdata, typename Tjudge>
	std::map<TA,std::map<TAC,Tensor<Tdata>>>
	comm_map2(
		const MPI_Comm &mpi_comm,
		const std::map<TA,std::map<TAC,Tensor<Tdata>>> &Ds_in,
		const Tjudge &judge)
	{
		Comm::Comm_Assemble<
			std::tuple<TA,TAC>,
			Tensor<Tdata>,
			std::map<TA,std::map<TAC,Tensor<Tdata>>>,
			Tjudge,
			std::map<TA,std::map<TAC,Tensor<Tdata>>>
		> com(mpi_comm);

		com.traverse_keys_provide = Comm::Communicate_Map::traverse_keys<TA,TAC,Tensor<Tdata>>;
		com.get_value_provide = Comm::Communicate_Map::get_value<TA,TAC,Tensor<Tdata>>;
		com.set_value_require = Comm::Communicate_Map::set_value_add<TA,TAC,Tensor<Tdata>>;
		com.flag_lock_set_value = Comm::Comm_Tools::Lock_Type::Copy_merge;
		com.init_datas_local = Comm::Communicate_Map::init_datas_local<TA,TAC,Tensor<Tdata>>;
		com.add_datas = Comm::Communicate_Map::add_datas<TA,TAC,Tensor<Tdata>>;

		std::map<TA,std::map<TAC,Tensor<Tdata>>> Ds_out;
		com.communicate( Ds_in, judge, Ds_out );
		return Ds_out;
	}

	template<typename TA, typename TAC, typename Tdata, typename Tjudge>
	std::map<TA,std::map<TAC,std::map<TAC,Tensor<Tdata>>>>
	comm_map3(
		const MPI_Comm &mpi_comm,
		const std::map<TA,std::map<TAC,std::map<TAC,Tensor<Tdata>>>> &Ds_in,
		const Tjudge &judge)
	{
		Comm::Comm_Assemble<
			std::tuple<TA,TAC,TAC>,
			Tensor<Tdata>,
			std::map<TA,std::map<TAC,std::map<TAC,Tensor<Tdata>>>>,
			Tjudge,
			std::map<TA,std::map<TAC,std::map<TAC,Tensor<Tdata>>>>
		> com(mpi_comm);

		com.traverse_keys_provide = Comm::Communicate_Map::traverse_keys<TA,TAC,TAC,Tensor<Tdata>>;
		com.get_value_provide = Comm::Communicate_Map::get_value<TA,TAC,TAC,Tensor<Tdata>>;
		com.set_value_require = Comm::Communicate_Map::set_value_add<TA,TAC,TAC,Tensor<Tdata>>;
		com.flag_lock_set_value = Comm::Comm_Tools::Lock_Type::Copy_merge;
		com.init_datas_local = Comm::Communicate_Map::init_datas_local<TA,TAC,TAC,Tensor<Tdata>>;
		com.add_datas = Comm::Communicate_Map::add_datas<TA,TAC,std::map<TAC,Tensor<Tdata>>>;

		std::map<TA,std::map<TAC,std::map<TAC,Tensor<Tdata>>>> Ds_out;
		com.communicate( Ds_in, judge, Ds_out );
		return Ds_out;
	}
}

}