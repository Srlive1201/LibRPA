// ==========================================
//  Author: Peize Lin, Rong Shi, Minye Zhang
//  Date:   2022.07.25
// ==========================================

#pragma once

#include "RPA.h"
#include "../ri/Label.h"
#include "./symmetry/Filter_Atom_Symmetry.h"

#include <cassert>

namespace RI
{

template<typename TA, typename Tcell, std::size_t Ndim, typename Tdata>
void RPA<TA,Tcell,Ndim,Tdata>::set_parallel(
	const MPI_Comm &mpi_comm,
	const std::map<TA,Tatom_pos> &atoms_pos,
	const std::array<Tatom_pos,Ndim> &latvec,
	const std::array<Tcell,Ndim> &period)
{
	this->lri.set_parallel(
		mpi_comm, atoms_pos, latvec, period,
		{Label::ab_ab::a1b1_a2b2, Label::ab_ab::a1b2_a2b1});
	this->flag_finish.stru = true;
}

template<typename TA, typename Tcell, std::size_t Ndim, typename Tdata>
void RPA<TA,Tcell,Ndim,Tdata>::set_symmetry(
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
void RPA<TA,Tcell,Ndim,Tdata>::set_Cs(
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
void RPA<TA,Tcell,Ndim,Tdata>::free_Cs(const std::string &save_name_suffix)
{
	this->lri.free_tensors_map2("Cs_"+save_name_suffix);
	this->flag_finish.Cs = false;
}

template<typename TA, typename Tcell, std::size_t Ndim, typename Tdata>
void RPA<TA,Tcell,Ndim,Tdata>::set_Gs_pos(
	const std::map<TA, std::map<TAC, Tensor<Tdata>>> &Gs_pos,
	const Tdata_real &threshold,
	const std::string &save_name_suffix)
{
	this->lri.set_tensors_map2(
		Gs_pos,
		{Label::ab::a1b1, Label::ab::a1b2, Label::ab::a2b1, Label::ab::a2b2},
		{{"threshold_filter", threshold}},
		"Gs_pos_"+save_name_suffix );
	this->flag_finish.Gs_pos = true;
}
template<typename TA, typename Tcell, std::size_t Ndim, typename Tdata>
void RPA<TA,Tcell,Ndim,Tdata>::free_Gs_pos(const std::string &save_name_suffix)
{
	this->lri.free_tensors_map2("Gs_pos_"+save_name_suffix);
	this->flag_finish.Gs_pos = false;
}

template<typename TA, typename Tcell, std::size_t Ndim, typename Tdata>
void RPA<TA,Tcell,Ndim,Tdata>::set_Gs_neg(
	const std::map<TA, std::map<TAC, Tensor<Tdata>>> &Gs_neg,
	const Tdata_real &threshold,
	const std::string &save_name_suffix)
{
	this->lri.set_tensors_map2(
		Gs_neg,
		{Label::ab::a1b1, Label::ab::a1b2, Label::ab::a2b1, Label::ab::a2b2},
		{{"threshold_filter", threshold}},
		"Gs_neg_"+save_name_suffix );
	this->flag_finish.Gs_neg = true;
}
template<typename TA, typename Tcell, std::size_t Ndim, typename Tdata>
void RPA<TA,Tcell,Ndim,Tdata>::free_Gs_neg(const std::string &save_name_suffix)
{
	this->lri.free_tensors_map2("Gs_neg_"+save_name_suffix);
	this->flag_finish.Gs_neg = false;
}


template<typename TA, typename Tcell, std::size_t Ndim, typename Tdata>
void RPA<TA,Tcell,Ndim,Tdata>::cal_chi0s(
	const std::array<std::string,3> &save_names_suffix)						// "Cs","Gs_pos","Gs_neg"
{
	assert(this->flag_finish.stru);
	assert(this->flag_finish.Cs);
	assert(this->flag_finish.Gs_pos);
	assert(this->flag_finish.Gs_neg);

	for(const Label::ab label : {Label::ab::a, Label::ab::b})
		this->lri.data_ab_name[label] = "Cs_"+save_names_suffix[0];

	this->chi0s.clear();

	for(const Label::ab label : {Label::ab::a1b1, Label::ab::a1b2})
		this->lri.data_ab_name[label] = "Gs_pos_"+save_names_suffix[1];
	for(const Label::ab label : {Label::ab::a2b1, Label::ab::a2b2})
		this->lri.data_ab_name[label] = "Gs_neg_"+save_names_suffix[2];

	this->lri.cal_loop3({
		Label::ab_ab::a1b1_a2b2,
		Label::ab_ab::a1b2_a2b1},
		this->chi0s);

	// conj?
	for(const Label::ab label : {Label::ab::a1b1, Label::ab::a1b2})
		this->lri.data_ab_name[label] = "Gs_neg_"+save_names_suffix[2];
	for(const Label::ab label : {Label::ab::a2b1, Label::ab::a2b2})
		this->lri.data_ab_name[label] = "Gs_pos_"+save_names_suffix[1];

	this->lri.cal_loop3({
		Label::ab_ab::a1b1_a2b2,
		Label::ab_ab::a1b2_a2b1},
		this->chi0s);
}

}