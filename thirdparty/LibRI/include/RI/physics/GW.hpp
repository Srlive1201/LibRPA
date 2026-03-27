// ===================
//  Author: Minye Zhang, almost completely copied from Exx.hpp
//  date: 2022.12.18
// ===================
#pragma once

#include "GW.h"
#include "../ri/Label.h"
#include "./symmetry/Filter_Atom_Symmetry.h"

namespace RI
{

template<typename TA, typename Tcell, std::size_t Ndim, typename Tdata>
void GW<TA,Tcell,Ndim,Tdata>::set_parallel(
	const MPI_Comm &mpi_comm,
	const std::map<TA,Tatom_pos> &atoms_pos,
	const std::array<Tatom_pos,Ndim> &latvec,
	const std::array<Tcell,Ndim> &period)
{
	this->lri.set_parallel(
		mpi_comm, atoms_pos, latvec, period,
		{Label::ab_ab::a0b0_a1b1, Label::ab_ab::a0b0_a1b2, Label::ab_ab::a0b0_a2b1, Label::ab_ab::a0b0_a2b2});
	this->flag_finish.stru = true;
	//if()
		// this->post_2D.set_parallel(this->mpi_comm, this->atoms_pos, this->period);
}

template<typename TA, typename Tcell, std::size_t Ndim, typename Tdata>
void GW<TA,Tcell,Ndim,Tdata>::set_symmetry(
	const bool flag_symmetry,
	const std::map<std::pair<TA,TA>, std::set<TC>> &irreducible_sector)
{
	if(flag_symmetry)
		this->lri.filter_atom = std::make_shared<Filter_Atom_Symmetry<TA,TC,Tdata>>(
			this->lri.period, irreducible_sector);
	else
		this->lri.filter_atom = std::make_shared<Filter_Atom<TA,TAC>>();
}

template<typename TA, typename Tcell, std::size_t Ndim, typename Tdata>
void GW<TA,Tcell,Ndim,Tdata>::set_Cs(
	const std::map<TA, std::map<TAC, Tensor<Tdata>>> &Cs,
	const Tdata_real &threshold,
	const std::string &save_name_suffix)
{
	this->lri.set_tensors_map2(
		Cs,
		{Label::ab::a, Label::ab::b},
		{{"threshold_filter", threshold}},
		"Cs_"+save_name_suffix );
	this->flag_finish.Cs = true;
}
template<typename TA, typename Tcell, std::size_t Ndim, typename Tdata>
void GW<TA,Tcell,Ndim,Tdata>::free_Cs(const std::string &save_name_suffix)
{
	this->lri.free_tensors_map2("Cs_"+save_name_suffix);
	this->flag_finish.Cs = false;
}

template<typename TA, typename Tcell, std::size_t Ndim, typename Tdata>
void GW<TA,Tcell,Ndim,Tdata>::set_Ws(
	const std::map<TA, std::map<TAC, Tensor<Tdata>>> &Ws,
	const Tdata_real &threshold,
	const std::string &save_name_suffix)
{
	this->lri.set_tensors_map2(
		Ws,
		{Label::ab::a0b0},
		{{"threshold_filter", threshold}},
		"Ws_"+save_name_suffix );
	this->flag_finish.Ws = true;
}
template<typename TA, typename Tcell, std::size_t Ndim, typename Tdata>
void GW<TA,Tcell,Ndim,Tdata>::free_Ws(const std::string &save_name_suffix)
{
	this->lri.free_tensors_map2("Ws_"+save_name_suffix);
	this->flag_finish.Ws = false;
}

template<typename TA, typename Tcell, std::size_t Ndim, typename Tdata>
void GW<TA,Tcell,Ndim,Tdata>::set_Gs(
	const std::map<TA, std::map<TAC, Tensor<Tdata>>> &Gs,
	const Tdata_real &threshold,
	const std::string &save_name_suffix)
{
	this->lri.set_tensors_map2(
		Gs,
		{Label::ab::a1b1, Label::ab::a1b2, Label::ab::a2b1, Label::ab::a2b2},
		{{"threshold_filter", threshold}},
		"Gs_"+save_name_suffix );
	this->flag_finish.Gs = true;
}
template<typename TA, typename Tcell, std::size_t Ndim, typename Tdata>
void GW<TA,Tcell,Ndim,Tdata>::free_Gs(const std::string &save_name_suffix)
{
	this->lri.free_tensors_map2("Gs_"+save_name_suffix);
	this->flag_finish.Gs = false;
}

template <typename TA, typename Tcell, std::size_t Ndim, typename Tdata>
void GW<TA, Tcell, Ndim, Tdata>::cal_Sigmas(
	const std::array<std::string,3> &save_names_suffix)						// "Cs","Ws","Gs"
{
	assert(this->flag_finish.stru);

	assert(this->flag_finish.Cs);
	for(const Label::ab label : {Label::ab::a, Label::ab::b})
		this->lri.data_ab_name[label] = "Cs_"+save_names_suffix[0];

	assert(this->flag_finish.Ws);
	this->lri.data_ab_name[Label::ab::a0b0] = "Ws_"+save_names_suffix[1];

	assert(this->flag_finish.Gs);
	for(const Label::ab label : {Label::ab::a1b1, Label::ab::a1b2, Label::ab::a2b1, Label::ab::a2b2})
		this->lri.data_ab_name[label] = "Gs_"+save_names_suffix[2];

	this->Sigmas.clear();
	this->lri.cal_loop3(
		{Label::ab_ab::a0b0_a1b1,
		 Label::ab_ab::a0b0_a1b2,
		 Label::ab_ab::a0b0_a2b1,
		 Label::ab_ab::a0b0_a2b2},
		this->Sigmas);
}

} // namespace RI