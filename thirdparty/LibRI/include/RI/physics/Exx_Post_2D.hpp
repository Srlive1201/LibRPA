#pragma once

#include "Exx_Post_2D.h"

#include "../global/Global_Func-1.h"
#include "../global/Map_Operator-2.h"
#include "../global/Map_Operator-3.h"
#include "../global/MPI_Wrapper.h"
#include "../comm/mix/Communicate_Tensors_Map_Judge.h"
#include "../distribute/Distribute_Equally.h"

#include "Comm/Comm_Trans/Comm_Trans.h"
#include "Comm/example/Communicate_Map-1.h"

#include <vector>
#include <mpi.h>
#include <string>
#include <stdexcept>

#define MPI_CHECK(x) if((x)!=MPI_SUCCESS)	throw std::runtime_error(std::string(__FILE__)+" line "+std::to_string(__LINE__));

namespace RI
{

template<typename TA, typename TC, typename Tdata> template<typename Tatom_pos>
void Exx_Post_2D<TA,TC,Tdata>::set_parallel(
	const MPI_Comm &mpi_comm_in,
	const std::map<TA,Tatom_pos> &atoms_pos,
	const TC &period)
{
	this->mpi_comm = mpi_comm_in;

	constexpr std::size_t num_index = 2;
	const std::vector<TA> atoms_vec = Global_Func::map_key_to_vec(atoms_pos);

	const std::pair<std::vector<TA>, std::vector<std::vector<TAC>>>
		atoms_split_list = Distribute_Equally::distribute_atoms_periods(
			this->mpi_comm, atoms_vec, period, num_index, false);

	this->list_Aa = Global_Func::to_set(atoms_split_list.first);
	this->list_Ab = Global_Func::to_set(atoms_split_list.second[0]);
}


template<typename TA, typename TC, typename Tdata>
auto Exx_Post_2D<TA,TC,Tdata>::set_tensors_map2( const std::map<TA,std::map<TAC,Tensor<Tdata>>> &Ds_in ) const
-> std::map<TA,std::map<TAC,Tensor<Tdata>>>
{
	return Communicate_Tensors_Map_Judge::comm_map2(this->mpi_comm, Ds_in, this->list_Aa, this->list_Ab);
}


template<typename TA, typename TC, typename Tdata>
auto Exx_Post_2D<TA,TC,Tdata>::cal_zip_dotc(
	const std::map<TA,std::map<TAC,Tensor<Tdata>>> &Ds,
	const std::map<TA,std::map<TAC,Tensor<Tdata>>> &Hs) const
-> std::map<TA, std::map<TAC, Tdata>>
{
	typedef Tdata(*Tfunc_dotc)(const Tensor<Tdata>&, const Tensor<Tdata>&);
	const std::function<Tdata(const Tensor<Tdata>&, const Tensor<Tdata>&)>
		dotc = static_cast<Tfunc_dotc>(Blas_Interface::dotc);
	const std::map<TA, std::map<TAC, Tdata>>
		E_map = Map_Operator::zip_intersection(Ds, Hs, dotc);
	return E_map;
}


template<typename TA, typename TC, typename Tdata>
Tdata Exx_Post_2D<TA,TC,Tdata>::cal_energy(
	const std::map<TA,std::map<TAC,Tensor<Tdata>>> &Ds,
	const std::map<TA,std::map<TAC,Tensor<Tdata>>> &Hs) const
{
	const std::map<TA, std::map<TAC, Tdata>>
		E_map = cal_zip_dotc(Ds, Hs);

	const std::function<Tdata(const Tdata&, const Tdata&)>
		plus = std::plus<Tdata>();
	Tdata E = Map_Operator::reduce(E_map, Tdata(0), plus);

	MPI_Wrapper::mpi_allreduce( E, MPI_SUM, this->mpi_comm );
	return E;
}

template<typename TA, typename TC, typename Tdata>
void Exx_Post_2D<TA,TC,Tdata>::cal_force(
	const std::map<TA,std::map<TAC,Tensor<Tdata>>> &Ds,
	const std::map<TA,std::map<TAC,Tensor<Tdata>>> &Hs,
	const bool flag_add,
	std::map<TA,Tdata> &force) const
{
	const std::map<TA, std::map<TAC, Tdata>>
		E_map = cal_zip_dotc(Ds, Hs);

	for(const auto &E_map_A : E_map)
		for(const auto &E_map_B : E_map_A.second)
			if(flag_add)
				force[E_map_A.first] += E_map_B.second;
			else
				force[std::get<0>(E_map_B.first)] += E_map_B.second;
}

template<typename TA, typename TC, typename Tdata>
std::map<TA,Tdata> Exx_Post_2D<TA,TC,Tdata>::reduce_force(
	const std::map<TA,Tdata> &F_local) const
{
	MPI_CHECK( MPI_Barrier(this->mpi_comm) );
	// add MPI_Barrier() to cover up an upresolved error.
	// Error message:
	// 		terminate called after throwing an instance of 'cereal::Exception'
	// 		what():  Failed to read 12 bytes from input stream! Read 4

	Comm::Comm_Trans<TA, Tdata, std::map<TA,Tdata>, std::map<TA,Tdata>> comm(this->mpi_comm);
	comm.traverse_isend = Comm::Communicate_Map::traverse_datas_all<TA,Tdata>;
	comm.set_value_recv = Comm::Communicate_Map::set_value_add<TA,Tdata>;
	comm.flag_lock_set_value = Comm::Comm_Tools::Lock_Type::Copy_merge;
	comm.init_datas_local = Comm::Communicate_Map::init_datas_local<TA,Tdata>;
	comm.add_datas = Comm::Communicate_Map::add_datas<TA,Tdata>;

	std::map<TA,Tdata> F_global;
	comm.communicate(F_local, F_global);
	return F_global;
}

}

#undef MPI_CHECK