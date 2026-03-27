// ===================
//  Author: Peize Lin
//  date: 2022.06.02
// ===================

#pragma once

#include "LRI.h"
#include "RI_Tools.h"
#include "Label_Tools.h"
#include "CS_Matrix.h"
#include "../global/Map_Operator.h"
#include <algorithm>

namespace RI
{

template<typename TA, typename Tcell, std::size_t Ndim, typename Tdata>
void LRI<TA,Tcell,Ndim,Tdata>::set_parallel(
	const MPI_Comm &mpi_comm_in,
	const std::map<TA,Tatom_pos> &atoms_pos,
	const std::array<Tatom_pos,Ndim> &latvec,
	const std::array<Tcell,Ndim> &period_in,
	const std::vector<Label::ab_ab> &labels_all)
{
	this->mpi_comm = mpi_comm_in;
	this->period = period_in;
	this->parallel->set_parallel(
		this->mpi_comm, atoms_pos, latvec, this->period,
		Label_Tools::to_Aab_Aab_set(labels_all));
}

/*
template<typename TA, typename Tcell, std::size_t Ndim, typename Tdata>
void LRI<TA,Tcell,Ndim,Tdata>::set_tensors_map2(
	const std::map<TA, std::map<TAC, Tensor<Tdata>>> &Ds_local,
	const Label::ab &label,
	const std::map<std::string, double> &para_in,
	const std::string &save_name_in)
{
	const std::map<std::string, double> para_default = {
		{"flag_period",      true},
		{"flag_comm",        true},
		{"flag_filter",      true},
		{"threshold_filter", 0.0}};
	const std::map<std::string, double> para = Map_Operator::cover(para_default, para_in);

	const std::string save_name =
		save_name_in=="default"
		? Label_Tools::get_name(label)
		: save_name_in;
	this->data_ab_name[label] = save_name;

	std::map<TA, std::map<TAC, Tensor<Tdata>>> Ds_new =
		para.at("flag_period")
		? RI_Tools::cal_period(Ds_local, this->period)
		: Ds_local;

	if(para.at("flag_comm"))
		Ds_new = this->parallel->comm_tensors_map2(label, std::move(Ds_new));

	if(para.at("flag_filter"))
		Ds_new = RI_Tools::filter(std::move(Ds_new), filter_funcs[label], para.at("threshold_filter"));

	this->data_pool[save_name].Ds_ab = std::move(Ds_new);

	this->data_pool[save_name].index_Ds_ab = RI_Tools::get_index(this->data_pool[save_name].Ds_ab);

	this->data_pool[save_name].csm_uplimits = this->csm.cal_uplimits(label, this->data_pool[save_name].Ds_ab);
}
*/

template<typename TA, typename Tcell, std::size_t Ndim, typename Tdata>
void LRI<TA,Tcell,Ndim,Tdata>::set_tensors_map2(
	const std::map<TA, std::map<TAC, Tensor<Tdata>>> &Ds_local,
	const std::vector<Label::ab> &label_list,
	const std::map<std::string, double> &para_in,
	const std::string &save_name_in)
{
	const std::map<std::string, double> para_default = {
		{"flag_period",      true},
		{"flag_comm",        true},
		{"flag_filter",      true},
		{"threshold_filter", 0.0}};
	const std::map<std::string, double> para = Map_Operator::cover(para_default, para_in);

	std::map<TA, std::map<TAC, Tensor<Tdata>>> Ds_new =
		para.at("flag_period")
		? RI_Tools::cal_period(Ds_local, this->period)
		: Ds_local;

	if(para.at("flag_comm"))
		Ds_new = this->parallel->comm_tensors_map2(label_list, std::move(Ds_new));

	if(para.at("flag_filter"))
	{
		std::vector<RI_Tools::T_filter_func<Tdata>> filter_func_list;
		for(const Label::ab &label : label_list)
			filter_func_list.push_back(this->filter_funcs[label]);
		Ds_new = RI_Tools::filter(std::move(Ds_new), filter_func_list, para.at("threshold_filter"));
	}

	const std::string save_name =
		save_name_in!="default"
		? save_name_in
		: Label_Tools::get_name(label_list);
	for(const Label::ab &label : label_list)
		this->data_ab_name[label] = save_name;

	this->data_pool[save_name].Ds_ab = std::move(Ds_new);

	this->data_pool[save_name].index_Ds_ab = RI_Tools::get_index(this->data_pool[save_name].Ds_ab);

	this->data_pool[save_name].csm_uplimits = this->csm.cal_uplimits(label_list, this->data_pool[save_name].Ds_ab);
}

/*
template<typename TA, typename Tcell, std::size_t Ndim, typename Tdata>
void RI<TA,Tcell,Ndim,Tdata>::set_tensors_map3(
	const std::map<TA, std::map<TAC, std::map<TAC, Tensor<Tdata>>>> &Ds_local,
	const Label::ab &label,
	const Tdata_real &threshold)
{
	if(label==Label::ab::a)
	{
		//this->Ds_ab[label] = Communicate::communicate(Ds_ab, threshold);
		if(threshold)
			this->Ds_a = LRI::filter(Ds_local, filter_func, threshold);
		else
			this->Ds_a = Ds_local;
		//this->cs_matrix.set_tensor2(this->Ds_ab[label], label);
	}
	else if(label==Label::ab::b)
	{
		//this->Ds_ab[label] = Communicate::communicate(Ds_ab, threshold);
		if(threshold)
			this->Ds_b = LRI::filter(Ds_local, filter_func, threshold);
		else
			this->Ds_b = Ds_local;
		//this->cs_matrix.set_tensor2(this->Ds_ab[label], label);
	}
}
*/


/*
template<typename TA, typename Tcell, std::size_t Ndim, typename Tdata>
std::map<TA,std::map<TAC,std::map<TAC,Tensor<Tdata>>>>
RI<TA,Tcell,Ndim,Tdata>::comm_tensors_map3(
	const Label::ab &label,
	const std::map<TA,std::map<TAC,Tensor<Tdata>>> &Ds) const
{
	switch(label)
	{
		case Label::ab:a:
			return Map3(this->mpi_comm, Ds, this->list_Aa0(), this->list_Aa1(), this->list_Aa2());
		case Label::ab:b:
			return Map3_period(this->mpi_comm, Ds, this->list_Ab0(), this->list_Ab1(), this->list_Ab2());
		default:
			throw std::invalid_argument(std::string(__FILE__)+" line "+std::to_string(__LINE__));
	}
}
*/

}