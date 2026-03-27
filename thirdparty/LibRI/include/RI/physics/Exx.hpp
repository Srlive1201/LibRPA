// ===================
//  Author: Peize Lin
//  date: 2022.06.02
// ===================

#pragma once

#include "Exx.h"
#include "../ri/Cell_Nearest.h"
#include "../ri/Label.h"
#include "../global/Map_Operator.h"
#include "./symmetry/Filter_Atom_Symmetry.h"

#include <cassert>

namespace RI
{

template<typename TA, typename Tcell, std::size_t Ndim, typename Tdata>
void Exx<TA,Tcell,Ndim,Tdata>::set_parallel(
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
		this->post_2D.set_parallel(mpi_comm, atoms_pos, period);
}

template<typename TA, typename Tcell, std::size_t Ndim, typename Tdata>
void Exx<TA,Tcell,Ndim,Tdata>::set_symmetry(
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
void Exx<TA,Tcell,Ndim,Tdata>::set_Cs(
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
void Exx<TA,Tcell,Ndim,Tdata>::free_Cs(const std::string &save_name_suffix)
{
	this->lri.free_tensors_map2("Cs_"+save_name_suffix);
	this->flag_finish.Cs = false;
}

template<typename TA, typename Tcell, std::size_t Ndim, typename Tdata>
void Exx<TA,Tcell,Ndim,Tdata>::set_Vs(
	const std::map<TA, std::map<TAC, Tensor<Tdata>>> &Vs,
	const Tdata_real &threshold,
	const std::string &save_name_suffix)
{
	this->lri.set_tensors_map2(
		Vs,
		{Label::ab::a0b0},
		{{"threshold_filter", threshold}},
		"Vs_"+save_name_suffix );
	this->flag_finish.Vs = true;
}
template<typename TA, typename Tcell, std::size_t Ndim, typename Tdata>
void Exx<TA,Tcell,Ndim,Tdata>::free_Vs(const std::string &save_name_suffix)
{
	this->lri.free_tensors_map2("Vs_"+save_name_suffix);
	this->flag_finish.Vs = false;
}

template<typename TA, typename Tcell, std::size_t Ndim, typename Tdata>
void Exx<TA,Tcell,Ndim,Tdata>::set_Ds(
	const std::map<TA, std::map<TAC, Tensor<Tdata>>> &Ds,
	const Tdata_real &threshold,
	const std::string &save_name_suffix)
{
	this->lri.set_tensors_map2(
		Ds,
		{Label::ab::a1b1, Label::ab::a1b2, Label::ab::a2b1, Label::ab::a2b2},
		{{"threshold_filter", threshold}},
		"Ds_"+save_name_suffix );
	this->flag_finish.Ds = true;
	this->flag_finish.Ds_delta = false;

	//if()
		this->post_2D.saves["Ds_"+save_name_suffix] = this->post_2D.set_tensors_map2(Ds);
}
template<typename TA, typename Tcell, std::size_t Ndim, typename Tdata>
void Exx<TA,Tcell,Ndim,Tdata>::free_Ds(const std::string &save_name_suffix)
{
	this->lri.free_tensors_map2("Ds_"+save_name_suffix);
	this->flag_finish.Ds = false;
}

template<typename TA, typename Tcell, std::size_t Ndim, typename Tdata>
void Exx<TA,Tcell,Ndim,Tdata>::set_Ds_delta(
	const std::map<TA, std::map<TAC, Tensor<Tdata>>> &Ds,
	const Tdata_real &threshold,
	const std::string &save_name_suffix)
{
	using namespace Map_Operator;

	assert(flag_finish.Ds);
	this->lri.set_tensors_map2(
		Ds,
		{Label::ab::a1b1, Label::ab::a1b2, Label::ab::a2b1, Label::ab::a2b2},
		{{"flag_filter", false}},
		"Ds_tmp" );
	this->lri.set_tensors_map2(
		this->lri.data_pool["Ds_tmp"].Ds_ab - this->lri.data_pool["Ds_"+save_name_suffix].Ds_ab,
		{Label::ab::a1b1, Label::ab::a1b2, Label::ab::a2b1, Label::ab::a2b2},
		{{"flag_period", false}, {"flag_comm", false}, {"flag_filter", true}, {"threshold_filter", threshold}},
		"Ds_delta_"+save_name_suffix);
	this->lri.data_pool.erase("Ds_tmp");
	this->lri.set_tensors_map2(
		this->lri.data_pool["Ds_delta_"+save_name_suffix].Ds_ab + this->lri.data_pool["Ds_"+save_name_suffix].Ds_ab,
		{Label::ab::a1b1, Label::ab::a1b2, Label::ab::a2b1, Label::ab::a2b2},
		{{"flag_period", false}, {"flag_comm", false}, {"flag_filter", false}},
		"Ds_"+save_name_suffix);
	this->flag_finish.Ds_delta = true;
	this->flag_finish.Ds = true;

	//if()
		this->post_2D.saves["Ds_"+save_name_suffix] = this->post_2D.set_tensors_map2(Ds);
}
template<typename TA, typename Tcell, std::size_t Ndim, typename Tdata>
void Exx<TA,Tcell,Ndim,Tdata>::free_Ds_delta(const std::string &save_name_suffix)
{
	this->lri.free_tensors_map2("Ds_delta_"+save_name_suffix);
	this->flag_finish.Ds_delta = false;
}

template<typename TA, typename Tcell, std::size_t Ndim, typename Tdata>
void Exx<TA,Tcell,Ndim,Tdata>::set_dCs(
	const std::array<std::map<TA, std::map<TAC, Tensor<Tdata>>>,Npos> &dCs,
	const Tdata_real &threshold,
	const std::string &save_name_suffix)
{
	for(std::size_t ipos=0; ipos<Npos; ++ipos)
		this->lri.set_tensors_map2(
			dCs[ipos],
			{Label::ab::a, Label::ab::b},
			{{"threshold_filter", threshold}},
			"dCs_"+std::to_string(ipos)+"_"+save_name_suffix );
	this->flag_finish.dCs = true;
}
template<typename TA, typename Tcell, std::size_t Ndim, typename Tdata>
void Exx<TA,Tcell,Ndim,Tdata>::free_dCs(const std::string &save_name_suffix)
{
	for(std::size_t ipos=0; ipos<Npos; ++ipos)
		this->lri.free_tensors_map2("dCs_"+std::to_string(ipos)+"_"+save_name_suffix);
	this->flag_finish.dCs = false;
}

template<typename TA, typename Tcell, std::size_t Ndim, typename Tdata>
void Exx<TA,Tcell,Ndim,Tdata>::set_dVs(
	const std::array<std::map<TA, std::map<TAC, Tensor<Tdata>>>,Npos> &dVs,
	const Tdata_real &threshold,
	const std::string &save_name_suffix)
{
	for(std::size_t ipos=0; ipos<Npos; ++ipos)
		this->lri.set_tensors_map2(
			dVs[ipos],
			{Label::ab::a0b0},
			{{"threshold_filter", threshold}},
			"dVs_"+std::to_string(ipos)+"_"+save_name_suffix );
	this->flag_finish.dVs = true;
}
template<typename TA, typename Tcell, std::size_t Ndim, typename Tdata>
void Exx<TA,Tcell,Ndim,Tdata>::free_dVs(const std::string &save_name_suffix)
{
	for(std::size_t ipos=0; ipos<Npos; ++ipos)
		this->lri.free_tensors_map2("dVs_"+std::to_string(ipos)+"_"+save_name_suffix);
	this->flag_finish.dVs = false;
}

template<typename TA, typename Tcell, std::size_t Ndim, typename Tdata>
void Exx<TA,Tcell,Ndim,Tdata>::set_dCRs(
	const std::array<std::array<std::map<TA, std::map<TAC, Tensor<Tdata>>>,Npos>,Npos> &dCRs,
	const Tdata_real &threshold,
	const std::string &save_name_suffix)
{
	for(std::size_t ipos0=0; ipos0<Npos; ++ipos0)
		for(std::size_t ipos1=0; ipos1<Npos; ++ipos1)
			this->lri.set_tensors_map2(
				dCRs[ipos0][ipos1],
				{Label::ab::a, Label::ab::b},
				{{"threshold_filter", threshold}},
				"dCRs_"+std::to_string(ipos0)+"_"+std::to_string(ipos1)+"_"+save_name_suffix );
	this->flag_finish.dCRs = true;
}
template<typename TA, typename Tcell, std::size_t Ndim, typename Tdata>
void Exx<TA,Tcell,Ndim,Tdata>::free_dCRs(const std::string &save_name_suffix)
{
	for(std::size_t ipos0=0; ipos0<Npos; ++ipos0)
		for(std::size_t ipos1=0; ipos1<Npos; ++ipos1)
			this->lri.free_tensors_map2("dCRs_"+std::to_string(ipos0)+"_"+std::to_string(ipos1)+"_"+save_name_suffix);
	this->flag_finish.dCRs = false;
}

template<typename TA, typename Tcell, std::size_t Ndim, typename Tdata>
void Exx<TA,Tcell,Ndim,Tdata>::set_dVRs(
	const std::array<std::array<std::map<TA, std::map<TAC, Tensor<Tdata>>>,Npos>,Npos> &dVRs,
	const Tdata_real &threshold,
	const std::string &save_name_suffix)
{
	for(std::size_t ipos0=0; ipos0<Npos; ++ipos0)
		for(std::size_t ipos1=0; ipos1<Npos; ++ipos1)
			this->lri.set_tensors_map2(
				dVRs[ipos0][ipos1],
				{Label::ab::a0b0},
				{{"threshold_filter", threshold}},
				"dVRs_"+std::to_string(ipos0)+"_"+std::to_string(ipos1)+"_"+save_name_suffix );
	this->flag_finish.dVRs = true;
}
template<typename TA, typename Tcell, std::size_t Ndim, typename Tdata>
void Exx<TA,Tcell,Ndim,Tdata>::free_dVRs(const std::string &save_name_suffix)
{
	for(std::size_t ipos0=0; ipos0<Npos; ++ipos0)
		for(std::size_t ipos1=0; ipos1<Npos; ++ipos1)
			this->lri.free_tensors_map2("dVRs_"+std::to_string(ipos0)+"_"+std::to_string(ipos1)+"_"+save_name_suffix);
	this->flag_finish.dVRs = false;
}


template<typename TA, typename Tcell, std::size_t Ndim, typename Tdata>
void Exx<TA,Tcell,Ndim,Tdata>::cal_Hs(
	const std::array<std::string,3> &save_names_suffix)						// "Cs","Vs","Ds"
{
	using namespace Map_Operator;

	assert(this->flag_finish.stru);

	assert(this->flag_finish.Cs);
	this->lri.data_ab_name[Label::ab::a] = this->lri.data_ab_name[Label::ab::b] = "Cs_"+save_names_suffix[0];

	assert(this->flag_finish.Vs);
	this->lri.data_ab_name[Label::ab::a0b0] = "Vs_"+save_names_suffix[1];

	if(!this->flag_finish.Ds_delta)
	{
		assert(this->flag_finish.Ds);
		for(const Label::ab label : {Label::ab::a1b1, Label::ab::a1b2, Label::ab::a2b1, Label::ab::a2b2})
			this->lri.data_ab_name[label] = "Ds_"+save_names_suffix[2];
	}
	else
	{
		for(const Label::ab label : {Label::ab::a1b1, Label::ab::a1b2, Label::ab::a2b1, Label::ab::a2b2})
			this->lri.data_ab_name[label] = "Ds_delta_"+save_names_suffix[2];
	}

	if(!this->flag_finish.Ds_delta)
		this->Hs.clear();
	this->lri.cal_loop3(
		{Label::ab_ab::a0b0_a1b1,
		 Label::ab_ab::a0b0_a1b2,
		 Label::ab_ab::a0b0_a2b1,
		 Label::ab_ab::a0b0_a2b2},
		this->Hs);

	//if()
		this->energy = this->post_2D.cal_energy(
			this->post_2D.saves["Ds_"+save_names_suffix[2]],
			this->post_2D.set_tensors_map2(this->Hs) );

	if(!this->flag_save_result.Hs)
		this->Hs.clear();
}

template<typename TA, typename Tcell, std::size_t Ndim, typename Tdata>
void Exx<TA,Tcell,Ndim,Tdata>::cal_force(
	const std::array<std::string,5> &save_names_suffix)						// "Cs","Vs","Ds","dCs","dVs"
{
	assert(this->flag_finish.stru);
	assert(this->flag_finish.Cs);
	assert(this->flag_finish.Vs);
	assert(this->flag_finish.Ds);
	assert(this->flag_finish.dCs);
	assert(this->flag_finish.dVs);

	for(const Label::ab label : {Label::ab::a1b1, Label::ab::a1b2, Label::ab::a2b1, Label::ab::a2b2})
		this->lri.data_ab_name[label] = "Ds_"+save_names_suffix[2];
	for(std::size_t ipos=0; ipos<Npos; ++ipos)
	{
		std::map<TA,Tdata> force_ipos;

		{
			this->dHs[ipos][0].clear();

			this->lri.data_ab_name[Label::ab::a   ] = "dCs_"+std::to_string(ipos)+"_"+save_names_suffix[3];
			this->lri.data_ab_name[Label::ab::a0b0] = "Vs_"+save_names_suffix[1];
			this->lri.data_ab_name[Label::ab::b   ] = "Cs_"+save_names_suffix[0];

			this->lri.cal_loop3(
				{Label::ab_ab::a0b0_a1b1,
				 Label::ab_ab::a0b0_a1b2,},
				this->dHs[ipos][0],
				-1.0);

			this->lri.cal_loop3(
				{Label::ab_ab::a0b0_a2b1,
				 Label::ab_ab::a0b0_a2b2},
				this->dHs[ipos][0],
				1.0);

			this->lri.data_ab_name[Label::ab::a   ] = "Cs_"+save_names_suffix[0];
			this->lri.data_ab_name[Label::ab::a0b0] = "dVs_" +std::to_string(ipos)+"_"+save_names_suffix[4];

			this->lri.cal_loop3(
				{Label::ab_ab::a0b0_a2b2,
				 Label::ab_ab::a0b0_a2b1},
				this->dHs[ipos][0],
				1.0);

			this->post_2D.cal_force(
				this->post_2D.saves["Ds_"+save_names_suffix[2]],
				this->post_2D.set_tensors_map2(this->dHs[ipos][0]),
				true,
				force_ipos );

			//mul(D)
			//this->Fs[ipos] = this->post_2D.cal_F(dHs);
			//this->stress[ipos] = this->post_2D.cal_sttress(dHs);
			//this->Fs[ipos][I] = \sum_J \sum_{i,j} dHs(i,j) * D(i,j)

			if(!this->flag_save_result.dHs)
				this->dHs[ipos][0].clear();
		}

		{
			this->dHs[ipos][1].clear();

			this->lri.cal_loop3(
				{Label::ab_ab::a0b0_a2b2,
				 Label::ab_ab::a0b0_a1b2},
				this->dHs[ipos][1],
				-1.0);

			this->lri.data_ab_name[Label::ab::a0b0] = "Vs_"+save_names_suffix[1];
			this->lri.data_ab_name[Label::ab::b   ] = "dCs_"+std::to_string(ipos)+"_"+save_names_suffix[3];

			this->lri.cal_loop3(
				{Label::ab_ab::a0b0_a1b1,
				 Label::ab_ab::a0b0_a2b1},
				this->dHs[ipos][1],
				-1.0);

			this->lri.cal_loop3(
				{Label::ab_ab::a0b0_a1b2,
				 Label::ab_ab::a0b0_a2b2},
				this->dHs[ipos][1],
				1.0);

			this->post_2D.cal_force(
				this->post_2D.saves["Ds_"+save_names_suffix[2]],
				this->post_2D.set_tensors_map2(this->dHs[ipos][1]),
				false,
				force_ipos );

			//mul(D)
			//this->Fs[ipos] -= this->post_2D.cal_F(dHs);
			//this->stress[ipos] -= this->post_2D.cal_sttress(dHs);
			//this->Fs[ipos][J] = \sum_I \sum_{i,j} dHs(i,j) * D(i,j)

			if(!this->flag_save_result.dHs)
				this->dHs[ipos][1].clear();
		}
		this->force[ipos] = this->post_2D.reduce_force(force_ipos);
	} // end for(ipos)
}



template<typename TA, typename Tcell, std::size_t Ndim, typename Tdata>
void Exx<TA,Tcell,Ndim,Tdata>::cal_stress(
	const std::array<std::string,5> &save_names_suffix)						// "Cs","Vs","Ds","dCRs","dVRs"
{
	assert(this->flag_finish.stru);
	assert(this->flag_finish.Cs);
	assert(this->flag_finish.Vs);
	assert(this->flag_finish.Ds);
	assert(this->flag_finish.dCRs);
	assert(this->flag_finish.dVRs);

	this->stress = Tensor<Tdata>({Npos, Npos});

	this->lri.data_ab_name[Label::ab::a1b1] = this->lri.data_ab_name[Label::ab::a2b1] = "Ds_"+save_names_suffix[2];
	for(std::size_t ipos0=0; ipos0<Npos; ++ipos0)
		for(std::size_t ipos1=0; ipos1<Npos; ++ipos1)
		{
			this->dHRs[ipos0][ipos1].clear();

			this->lri.data_ab_name[Label::ab::a   ] = "dCRs_"+std::to_string(ipos0)+"_"+std::to_string(ipos1)+"_"+save_names_suffix[3];
			this->lri.data_ab_name[Label::ab::a0b0] = "Vs_"+save_names_suffix[1];
			this->lri.data_ab_name[Label::ab::b   ] = "Cs_"+save_names_suffix[0];

			this->lri.cal_loop3(
				{Label::ab_ab::a0b0_a1b1,
				Label::ab_ab::a0b0_a2b1},
				this->dHRs[ipos0][ipos1],
				1.0);

			this->lri.data_ab_name[Label::ab::a   ] = "Cs_"+save_names_suffix[0];
			this->lri.data_ab_name[Label::ab::a0b0] = "dVRs_"+std::to_string(ipos0)+"_"+std::to_string(ipos1)+"_"+save_names_suffix[4];

			this->lri.cal_loop3(
				{Label::ab_ab::a0b0_a1b1,
				Label::ab_ab::a0b0_a2b1},
				this->dHRs[ipos0][ipos1],
				1.0);

			this->lri.data_ab_name[Label::ab::a0b0] = "Vs_"+save_names_suffix[1];
			this->lri.data_ab_name[Label::ab::b   ] = "dCRs_"+std::to_string(ipos0)+"_"+std::to_string(ipos1)+"_"+save_names_suffix[3];

			this->lri.cal_loop3(
				{Label::ab_ab::a0b0_a1b1,
				Label::ab_ab::a0b0_a2b1},
				this->dHRs[ipos0][ipos1],
				1.0);

			this->stress(ipos0,ipos1) = post_2D.cal_energy(
				this->post_2D.saves["Ds_"+save_names_suffix[2]],
				this->post_2D.set_tensors_map2(this->dHRs[ipos0][ipos1]));

			if(!this->flag_save_result.dHRs)
				this->dHRs[ipos0][ipos1].clear();
		}
}

}