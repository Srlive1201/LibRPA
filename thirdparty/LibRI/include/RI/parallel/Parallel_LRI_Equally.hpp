// ===================
//  Author: Peize Lin
//  date: 2022.07.23
// ===================

#pragma once

#include "Parallel_LRI_Equally.h"
#include "../global/Global_Func-1.h"
#include "../distribute/Distribute_Equally.h"
#include "../comm/mix/Communicate_Tensors_Map_Judge.h"

namespace RI
{

template<typename TA, typename Tcell, std::size_t Ndim, typename Tdata>
void Parallel_LRI_Equally<TA,Tcell,Ndim,Tdata>::set_parallel(
	const MPI_Comm &mpi_comm_in,
	const std::map<TA,Tatom_pos> &atoms_pos,
	const std::array<Tatom_pos,Ndim> &latvec,
	const std::array<Tcell,Ndim> &period_in,
	const std::set<Label::Aab_Aab> &labels)
{
	this->mpi_comm = mpi_comm_in;
	this->period = period_in;
	const std::vector<TA> atoms_vec = Global_Func::map_key_to_vec(atoms_pos);

	this->set_parallel_loop4(atoms_vec);
	this->set_parallel_loop3(atoms_vec, labels);
}

template<typename TA, typename Tcell, std::size_t Ndim, typename Tdata>
void Parallel_LRI_Equally<TA,Tcell,Ndim,Tdata>::set_parallel_loop4(
	const std::vector<TA> &atoms_vec)
{
	constexpr std::size_t num_index = 4;

	const std::pair<std::vector<TA>, std::vector<std::vector<std::pair<TA,TC>>>>
		atoms_split_list = Distribute_Equally::distribute_atoms_periods(
			this->mpi_comm, atoms_vec, this->period, num_index, false);

	this->list_Aa01 = atoms_split_list.first;
	this->list_Aa2  = atoms_split_list.second[0];
	this->list_Ab01 = atoms_split_list.second[1];
	this->list_Ab2  = atoms_split_list.second[2];
}

template<typename TA, typename Tcell, std::size_t Ndim, typename Tdata>
void Parallel_LRI_Equally<TA,Tcell,Ndim,Tdata>::set_parallel_loop3(
	const std::vector<TA> &atoms_vec,
	const std::set<Label::Aab_Aab> &labels)
{
	constexpr std::size_t num_index = 2;
	const std::vector<TAC> atoms_period_vec = Divide_Atoms::traversal_atom_period(atoms_vec, this->period);

	const std::pair<std::vector<TA>, std::vector<std::vector<std::pair<TA,TC>>>>
		atoms_split_list1 = Distribute_Equally::distribute_atoms_periods(
			this->mpi_comm, atoms_vec, this->period, num_index, false);
	const std::vector<std::vector<std::pair<TA,TC>>>
		atoms_split_list2 = Distribute_Equally::distribute_periods(
			this->mpi_comm, atoms_vec, this->period, num_index, false);
	for(const Label::Aab_Aab &label : labels)
	{
		List_A<TA,TAC> &atoms = this->list_A[label];
		atoms.a01 = atoms_vec;
		atoms.a2 = atoms_period_vec;
		atoms.b01 = atoms_period_vec;
		atoms.b2 = atoms_period_vec;
		switch(label)
		{
			case Label::Aab_Aab::a01b01_a01b01:
				atoms.a2  = atoms_split_list2[0];
				atoms.b01 = atoms_split_list2[1];
				break;
			case Label::Aab_Aab::a01b01_a2b01:
				atoms.a01 = atoms_split_list1.first;
				atoms.b01 = atoms_split_list1.second[0];
				break;
			case Label::Aab_Aab::a01b01_a01b2:
				atoms.b01 = atoms_split_list1.second[0];
				atoms.a01 = atoms_split_list1.first;
				break;
			case Label::Aab_Aab::a01b01_a2b2:
				atoms.a01 = atoms_split_list1.first;
				atoms.b2  = atoms_split_list1.second[0];
				break;
			case Label::Aab_Aab::a01b2_a2b01:
				atoms.a01 = atoms_split_list1.first;
				atoms.b01 = atoms_split_list1.second[0];
				break;
			default:
				throw std::invalid_argument(std::string(__FILE__)+" line "+std::to_string(__LINE__));
		}
	}
}


/*
template<typename TA, typename Tcell, std::size_t Ndim, typename Tdata>
auto Parallel_LRI_Equally<TA,Tcell,Ndim,Tdata>::comm_tensors_map2(
	const Label::ab &label,
	const std::map<TA,std::map<TAC,Tensor<Tdata>>> &Ds) const
-> std::map<TA,std::map<TAC,Tensor<Tdata>>>
{
	switch(label)
	{
		case Label::ab::a:
			return Communicate_Tensors_Map_Judge::comm_map2(this->mpi_comm, Ds, Global_Func::to_set(this->list_Aa01), Global_Func::to_set(this->list_Aa2));
		case Label::ab::b:
			return Communicate_Tensors_Map_Judge::comm_map2_period(this->mpi_comm, Ds, Global_Func::to_set(this->list_Ab01), Global_Func::to_set(this->list_Ab2), this->period);
		case Label::ab::a0b0:	case Label::ab::a0b1:
		case Label::ab::a1b0:	case Label::ab::a1b1:
			return Communicate_Tensors_Map_Judge::comm_map2(this->mpi_comm, Ds, Global_Func::to_set(this->list_Aa01), Global_Func::to_set(this->list_Ab01));
		case Label::ab::a0b2:	case Label::ab::a1b2:
			return Communicate_Tensors_Map_Judge::comm_map2(this->mpi_comm, Ds, Global_Func::to_set(this->list_Aa01), Global_Func::to_set(this->list_Ab2));
		case Label::ab::a2b0:	case Label::ab::a2b1:
			return Communicate_Tensors_Map_Judge::comm_map2_period(this->mpi_comm, Ds, Global_Func::to_set(this->list_Aa2), Global_Func::to_set(this->list_Ab01), this->period);
		case Label::ab::a2b2:
			return Communicate_Tensors_Map_Judge::comm_map2_period(this->mpi_comm, Ds, Global_Func::to_set(this->list_Aa2), Global_Func::to_set(this->list_Ab2), this->period);
		default:
			throw std::invalid_argument(std::string(__FILE__)+" line "+std::to_string(__LINE__));
	}
}
*/

template<typename TA, typename Tcell, std::size_t Ndim, typename Tdata>
auto Parallel_LRI_Equally<TA,Tcell,Ndim,Tdata>::comm_tensors_map2(
	const std::vector<Label::ab> &label_list,
	const std::map<TA,std::map<TAC,Tensor<Tdata>>> &Ds) const
-> std::map<TA,std::map<TAC,Tensor<Tdata>>>
{
	std::tuple<
		std::vector<std::tuple< std::set<TA>, std::set<std::pair<TA,TC>> >>,
		std::vector<std::tuple< std::set<std::pair<TA,TC>>, std::set<std::pair<TA,TC>> >>
		> s_list;

	std::vector<bool> flags(6, false);
	for(const Label::ab &label : label_list)
	{
		switch(label)
		{
			case Label::ab::a:
				flags[0]=true;	break;
			case Label::ab::b:
				flags[1]=true;	break;
			case Label::ab::a0b0:	case Label::ab::a0b1:
			case Label::ab::a1b0:	case Label::ab::a1b1:
				flags[2]=true;	break;
			case Label::ab::a0b2:	case Label::ab::a1b2:
				flags[3]=true;	break;
			case Label::ab::a2b0:	case Label::ab::a2b1:
				flags[4]=true;	break;
			case Label::ab::a2b2:
				flags[5]=true;	break;
			default:
				throw std::invalid_argument(std::string(__FILE__)+" line "+std::to_string(__LINE__));
		}
	}

	// used in loop3
	for(const auto &list_atom : this->list_A)
	{
		if(flags[0])	std::get<0>(s_list).push_back(std::make_tuple( Global_Func::to_set(list_atom.second.a01), Global_Func::to_set(list_atom.second.a2)  ));
		if(flags[1])	std::get<1>(s_list).push_back(std::make_tuple( Global_Func::to_set(list_atom.second.b01), Global_Func::to_set(list_atom.second.b2)  ));
		if(flags[2])	std::get<0>(s_list).push_back(std::make_tuple( Global_Func::to_set(list_atom.second.a01), Global_Func::to_set(list_atom.second.b01) ));
		if(flags[3])	std::get<0>(s_list).push_back(std::make_tuple( Global_Func::to_set(list_atom.second.a01), Global_Func::to_set(list_atom.second.b2)  ));
		if(flags[4])	std::get<1>(s_list).push_back(std::make_tuple( Global_Func::to_set(list_atom.second.a2),  Global_Func::to_set(list_atom.second.b01) ));
		if(flags[5])	std::get<1>(s_list).push_back(std::make_tuple( Global_Func::to_set(list_atom.second.a2),  Global_Func::to_set(list_atom.second.b2)  ));
	}

	// used in loop4
	if(flags[0])	std::get<0>(s_list).push_back(std::make_tuple( Global_Func::to_set(this->list_Aa01), Global_Func::to_set(this->list_Aa2) ));
	if(flags[1])	std::get<1>(s_list).push_back(std::make_tuple( Global_Func::to_set(this->list_Ab01), Global_Func::to_set(this->list_Ab2) ));
	if(flags[2])	std::get<0>(s_list).push_back(std::make_tuple( Global_Func::to_set(this->list_Aa01), Global_Func::to_set(this->list_Ab01) ));
	if(flags[3])	std::get<0>(s_list).push_back(std::make_tuple( Global_Func::to_set(this->list_Aa01), Global_Func::to_set(this->list_Ab2) ));
	if(flags[4])	std::get<1>(s_list).push_back(std::make_tuple( Global_Func::to_set(this->list_Aa2), Global_Func::to_set(this->list_Ab01) ));
	if(flags[5])	std::get<1>(s_list).push_back(std::make_tuple( Global_Func::to_set(this->list_Aa2), Global_Func::to_set(this->list_Ab2) ));

	return Communicate_Tensors_Map_Judge::comm_map2_combine_origin_period(this->mpi_comm, Ds, s_list, this->period);
}

}