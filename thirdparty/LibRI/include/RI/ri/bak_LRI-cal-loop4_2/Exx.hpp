// ===================
//  Author: Peize Lin
//  date: 2022.06.02
// ===================

#pragma once

#include "Exx.h"
#include "../ri/Cell_Nearest.h"
#include "../ri/Label.h"
#include "../global/Map_Operator.h"

#include <cassert>

namespace RI
{

template<typename TA, typename Tcell, std::size_t Ndim, typename Tdata>
void Exx<TA,Tcell,Ndim,Tdata>::set_parallel(
	const MPI_Comm &mpi_comm_in,
	const std::map<TA,Tatom_pos> &atoms_pos_in,
	const std::array<Tatom_pos,Ndim> &latvec_in,
	const std::array<Tcell,Ndim> &period_in)
{
	this->mpi_comm = mpi_comm_in;
	this->atoms_pos = atoms_pos_in;
	this->latvec = latvec_in;
	this->period = period_in;

	this->lri.set_parallel(
		this->mpi_comm, this->atoms_pos, this->latvec, this->period,
		{Label::ab_ab::a0b0_a1b1, Label::ab_ab::a0b0_a1b2, Label::ab_ab::a0b0_a2b1, Label::ab_ab::a0b0_a2b2});
	this->flag_finish.stru = true;
	//if()
		this->post_2D.set_parallel(this->mpi_comm, this->atoms_pos, this->period);
}

template<typename TA, typename Tcell, std::size_t Ndim, typename Tdata>
void Exx<TA,Tcell,Ndim,Tdata>::set_Cs(
	const std::map<TA, std::map<TAC, Tensor<Tdata>>> &Cs,
	const Tdata_real &threshold_C,
	const std::string &save_name_suffix)
{
	this->lri.set_tensors_map2(
		Cs,
		{Label::ab::a, Label::ab::b},
		{{"threshold_filter", threshold_C}},
		"Cs_"+save_name_suffix );
	this->flag_finish.C = true;
}

template<typename TA, typename Tcell, std::size_t Ndim, typename Tdata>
void Exx<TA,Tcell,Ndim,Tdata>::set_Vs(
	const std::map<TA, std::map<TAC, Tensor<Tdata>>> &Vs,
	const Tdata_real &threshold_V,
	const std::string &save_name_suffix)
{
	this->lri.set_tensors_map2(
		Vs,
		{Label::ab::a0b0},
		{{"threshold_filter", threshold_V}},
		"Vs_"+save_name_suffix );
	this->flag_finish.V = true;
}

template<typename TA, typename Tcell, std::size_t Ndim, typename Tdata>
void Exx<TA,Tcell,Ndim,Tdata>::set_Ds(
	const std::map<TA, std::map<TAC, Tensor<Tdata>>> &Ds,
	const Tdata_real &threshold_D,
	const std::string &save_name_suffix)
{
	this->lri.set_tensors_map2(
		Ds,
		{Label::ab::a1b1, Label::ab::a1b2, Label::ab::a2b1, Label::ab::a2b2},
		{{"threshold_filter", threshold_D}},
		"Ds_"+save_name_suffix );
	this->flag_finish.D = true;
	this->flag_finish.D_delta = false;

	//if()
		this->post_2D.saves["Ds_"+save_name_suffix] = this->post_2D.set_tensors_map2(Ds);
}

template<typename TA, typename Tcell, std::size_t Ndim, typename Tdata>
void Exx<TA,Tcell,Ndim,Tdata>::set_Ds_delta(
	const std::map<TA, std::map<TAC, Tensor<Tdata>>> &Ds,
	const Tdata_real &threshold_D,
	const std::string &save_name_suffix)
{
	using namespace Map_Operator;

	assert(flag_finish.D);
	this->lri.set_tensors_map2(
		Ds,
		{Label::ab::a1b1, Label::ab::a1b2, Label::ab::a2b1, Label::ab::a2b2},
		{{"flag_filter", false}},
		"Ds_tmp" );
	this->lri.set_tensors_map2(
		this->lri.data_pool["Ds_tmp"].Ds_ab - this->lri.data_pool["Ds_"+save_name_suffix].Ds_ab,
		{Label::ab::a1b1, Label::ab::a1b2, Label::ab::a2b1, Label::ab::a2b2},
		{{"flag_period", false}, {"flag_comm", false}, {"flag_filter", true}, {"threshold_filter", threshold_D}},
		"Ds_delta_"+save_name_suffix);
	this->lri.data_pool.erase("Ds_tmp");
	this->lri.set_tensors_map2(
		this->lri.data_pool["Ds_delta_"+save_name_suffix].Ds_ab + this->lri.data_pool["Ds_"+save_name_suffix].Ds_ab,
		{Label::ab::a1b1, Label::ab::a1b2, Label::ab::a2b1, Label::ab::a2b2},
		{{"flag_period", false}, {"flag_comm", false}, {"flag_filter", false}},
		"Ds_"+save_name_suffix);
	this->flag_finish.D_delta = true;
	this->flag_finish.D = true;

	//if()
		this->post_2D.saves["Ds_"+save_name_suffix] = this->post_2D.set_tensors_map2(Ds);
}

template<typename TA, typename Tcell, std::size_t Ndim, typename Tdata>
void Exx<TA,Tcell,Ndim,Tdata>::set_dCs(
	const std::array<std::map<TA, std::map<TAC, Tensor<Tdata>>>,Npos> &dCs,
	const Tdata_real &threshold_dC,
	const std::string &save_name_suffix)
{
	for(std::size_t ipos=0; ipos<Npos; ++ipos)
		this->lri.set_tensors_map2(
			dCs[ipos],
			{Label::ab::a, Label::ab::b},
			{{"threshold_filter", threshold_dC}},
			"dCs_"+std::to_string(ipos)+"_"+save_name_suffix );
	this->flag_finish.dC = true;
}

template<typename TA, typename Tcell, std::size_t Ndim, typename Tdata>
void Exx<TA,Tcell,Ndim,Tdata>::set_dVs(
	const std::array<std::map<TA, std::map<TAC, Tensor<Tdata>>>,Npos> &dVs,
	const Tdata_real &threshold_dV,
	const std::string &save_name_suffix)
{
	for(std::size_t ipos=0; ipos<Npos; ++ipos)
		this->lri.set_tensors_map2(
			dVs[ipos],
			{Label::ab::a0b0},
			{{"threshold_filter", threshold_dV}},
			"dVs_"+std::to_string(ipos)+"_"+save_name_suffix );
	this->flag_finish.dV = true;
}
template<typename TA, typename Tcell, std::size_t Ndim, typename Tdata>
void Exx<TA,Tcell,Ndim,Tdata>::set_dCRs(
	const std::array<std::array<std::map<TA, std::map<TAC, Tensor<Tdata>>>,Npos>,Npos> &dCRs,
	const Tdata_real &threshold_dCR,
	const std::string &save_name_suffix)
{
	for(std::size_t ipos0=0; ipos0<Npos; ++ipos0)
		for(std::size_t ipos1=0; ipos1<Npos; ++ipos1)
			this->lri.set_tensors_map2(
				dCRs[ipos0][ipos1],
				{Label::ab::a, Label::ab::b},
				{{"threshold_filter", threshold_dCR}},
				"dCRs_"+std::to_string(ipos0)+"_"+std::to_string(ipos1)+"_"+save_name_suffix );
	this->flag_finish.dCR = true;
}

template<typename TA, typename Tcell, std::size_t Ndim, typename Tdata>
void Exx<TA,Tcell,Ndim,Tdata>::set_dVRs(
	const std::array<std::array<std::map<TA, std::map<TAC, Tensor<Tdata>>>,Npos>,Npos> &dVRs,
	const Tdata_real &threshold_dVR,
	const std::string &save_name_suffix)
{
	for(std::size_t ipos0=0; ipos0<Npos; ++ipos0)
		for(std::size_t ipos1=0; ipos1<Npos; ++ipos1)
			this->lri.set_tensors_map2(
				dVRs[ipos0][ipos1],
				{Label::ab::a0b0},
				{{"threshold_filter", threshold_dVR}},
				"dVRs_"+std::to_string(ipos0)+"_"+std::to_string(ipos1)+"_"+save_name_suffix );
	this->flag_finish.dVR = true;
}


template<typename TA, typename Tcell, std::size_t Ndim, typename Tdata>
void Exx<TA,Tcell,Ndim,Tdata>::cal_Hs(
	const std::array<std::string,3> &save_names_suffix)						// "Cs","Vs","Ds"
{
	using namespace Map_Operator;

	assert(this->flag_finish.stru);

	assert(this->flag_finish.C);
	this->lri.data_ab_name[Label::ab::a] = this->lri.data_ab_name[Label::ab::b] = "Cs_"+save_names_suffix[0];

	assert(this->flag_finish.V);
	this->lri.data_ab_name[Label::ab::a0b0] = "Vs_"+save_names_suffix[1];

	if(!this->flag_finish.D_delta)
	{
		assert(this->flag_finish.D);
		for(const Label::ab label : {Label::ab::a1b1, Label::ab::a1b2, Label::ab::a2b1, Label::ab::a2b2})
			this->lri.data_ab_name[label] = "Ds_"+save_names_suffix[2];
	}
	else
	{
		for(const Label::ab label : {Label::ab::a1b1, Label::ab::a1b2, Label::ab::a2b1, Label::ab::a2b2})
			this->lri.data_ab_name[label] = "Ds_delta_"+save_names_suffix[2];
	}

	std::vector<std::map<TA, std::map<TAC, Tensor<Tdata>>>> Hs_vec(1);
	this->lri.coefficients = {nullptr};
	this->lri.cal_loop3(
		{Label::ab_ab::a0b0_a1b1,
		 Label::ab_ab::a0b0_a1b2,
		 Label::ab_ab::a0b0_a2b1,
		 Label::ab_ab::a0b0_a2b2},
		Hs_vec);
	if(this->flag_finish.D_delta)
		this->Hs = this->Hs + Hs_vec[0];
	else
		this->Hs = std::move(Hs_vec[0]);

	//if()
		this->energy = this->post_2D.cal_energy(
			this->post_2D.saves["Ds_"+save_names_suffix[2]],
			this->post_2D.set_tensors_map2(this->Hs) );
}

template<typename TA, typename Tcell, std::size_t Ndim, typename Tdata>
void Exx<TA,Tcell,Ndim,Tdata>::cal_force(
	const std::array<std::string,5> &save_names_suffix)						// "Cs","Vs","Ds","dCs","dVs"
{
	assert(this->flag_finish.stru);
	assert(this->flag_finish.C);
	assert(this->flag_finish.V);
	assert(this->flag_finish.D);
	assert(this->flag_finish.dC);
	assert(this->flag_finish.dV);

	for(const Label::ab label : {Label::ab::a1b1, Label::ab::a1b2, Label::ab::a2b1, Label::ab::a2b2})
		this->lri.data_ab_name[label] = "Ds_"+save_names_suffix[2];
	for(std::size_t ipos=0; ipos<Npos; ++ipos)
	{
		std::map<TA,Tdata> force_ipos;

		{
			std::vector<std::map<TA,std::map<TAC,Tensor<Tdata>>>> dHs_vec(1);

			this->lri.data_ab_name[Label::ab::a   ] = "dCs_"+std::to_string(ipos)+"_"+save_names_suffix[3];
			this->lri.data_ab_name[Label::ab::a0b0] = "Vs_"+save_names_suffix[1];
			this->lri.data_ab_name[Label::ab::b   ] = "Cs_"+save_names_suffix[0];

			this->lri.cal_loop3(
				{Label::ab_ab::a0b0_a1b1,
				 Label::ab_ab::a0b0_a1b2,},
				dHs_vec,
				-1.0);

			this->lri.cal_loop3(
				{Label::ab_ab::a0b0_a2b1,
				 Label::ab_ab::a0b0_a2b2},
				dHs_vec,
				1.0);

			this->lri.data_ab_name[Label::ab::a   ] = "Cs_"+save_names_suffix[0];
			this->lri.data_ab_name[Label::ab::a0b0] = "dVs_" +std::to_string(ipos)+"_"+save_names_suffix[4];

			this->lri.cal_loop3(
				{Label::ab_ab::a0b0_a2b2,
				 Label::ab_ab::a0b0_a2b1},
				dHs_vec,
				1.0);

			this->post_2D.cal_force(
				this->post_2D.saves["Ds_"+save_names_suffix[2]],
				this->post_2D.set_tensors_map2(std::move(dHs_vec[0])),
				true,
				force_ipos );

//			mul(D)
//			this->Fs[ipos] = this->post_2D.cal_F(dHs);
//			this->stress[ipos] = this->post_2D.cal_sttress(dHs);
	//		this->Fs[ipos][I] = \sum_J \sum_{i,j} dHs(i,j) * D(i,j)
		}

		{
			std::vector<std::map<TA,std::map<TAC,Tensor<Tdata>>>> dHs_vec(1);

			this->lri.cal_loop3(
				{Label::ab_ab::a0b0_a2b2,
				 Label::ab_ab::a0b0_a1b2},
				dHs_vec,
				1.0);

			this->lri.data_ab_name[Label::ab::a0b0] = "Vs_"+save_names_suffix[1];
			this->lri.data_ab_name[Label::ab::b   ] = "dCs_"+std::to_string(ipos)+"_"+save_names_suffix[3];

			this->lri.cal_loop3(
				{Label::ab_ab::a0b0_a1b1,
				 Label::ab_ab::a0b0_a2b1},
				dHs_vec,
				1.0);

			this->lri.cal_loop3(
				{Label::ab_ab::a0b0_a1b2,
				 Label::ab_ab::a0b0_a2b2},
				dHs_vec,
				-1.0);

			this->post_2D.cal_force(
				this->post_2D.saves["Ds_"+save_names_suffix[2]],
				this->post_2D.set_tensors_map2(std::move(dHs_vec[0])),
				false,
				force_ipos );

//			mul(D)
//			this->Fs[ipos] -= this->post_2D.cal_F(dHs);
//			this->stress[ipos] -= this->post_2D.cal_sttress(dHs);
	//		this->Fs[ipos][J] = \sum_I \sum_{i,j} dHs(i,j) * D(i,j)
		}
		this->force[ipos] = this->post_2D.reduce_force(force_ipos);
	} // end for(ipos)
}


/*
template<typename TA, typename Tcell, std::size_t Ndim, typename Tdata>
void Exx<TA,Tcell,Ndim,Tdata>::cal_stress(
	const std::array<std::string,5> &save_names_suffix)						// "Cs","Vs","Ds","dCs","dVs"
{
	assert(this->flag_finish.stru);
	assert(this->flag_finish.C);
	assert(this->flag_finish.V);
	assert(this->flag_finish.D);
	assert(this->flag_finish.dC);
	assert(this->flag_finish.dV);

	Cell_Nearest<TA,Tcell,Ndim,Tpos,Npos> cell_nearest;
	cell_nearest.init(this->atoms_pos, this->latvec, this->period);

	using namespace Array_Operator;
	auto get_delta_pos = [this, &cell_nearest](
		const TA &Ax, const TA &Ay, const std::array<Tcell,Ndim> &celly, const std::size_t &ipos1)
		-> Tpos
	{
		Tpos delta_pos = this->atoms_pos[Ay][ipos1] - this->atoms_pos[Ax][ipos1];
		const std::array<Tcell,Ndim> celly_nearest = cell_nearest.get_cell_nearest_discrete(Ax, Ay, celly);
		for(std::size_t idim=0; idim<Ndim; ++idim)
			delta_pos += celly_nearest[idim] * this->latvec[idim][ipos1];
		return delta_pos;
	};

	this->stress = Tensor<Tdata>({Npos, Npos});

	this->lri.data_ab_name[Label::ab::a1b1] = this->lri.data_ab_name[Label::ab::a2b1] = "Ds_"+save_names_suffix[2];
	for(std::size_t ipos0=0; ipos0<Npos; ++ipos0)
	{
		std::vector<std::map<TA,std::map<TAC,Tensor<Tdata>>>> dHs_vec(Npos);
		this->lri.coefficients.resize(Npos);

		this->lri.data_ab_name[Label::ab::a   ] = "dCs_"+std::to_string(ipos0)+"_"+save_names_suffix[3];
		this->lri.data_ab_name[Label::ab::a0b0] = "Vs_"+save_names_suffix[1];
		this->lri.data_ab_name[Label::ab::b   ] = "Cs_"+save_names_suffix[0];

		for(std::size_t ipos1=0; ipos1<Npos; ++ipos1)
			this->lri.coefficients[ipos1] =
				[this,ipos1,&get_delta_pos](
					const Label::ab_ab &label, const TA &Aa01, const TAC &Aa2, const TAC &Ab01, const TAC &Ab2) -> Tdata
				{	return get_delta_pos(Aa01, Aa2.first, Aa2.second, ipos1);	};

		this->lri.cal(
			{Label::ab_ab::a0b0_a1b1,
			 Label::ab_ab::a0b0_a2b1},
			dHs_vec);

		this->lri.data_ab_name[Label::ab::a   ] = "Cs_"+save_names_suffix[0];
		this->lri.data_ab_name[Label::ab::a0b0] = "dVs_"+std::to_string(ipos0)+"_"+save_names_suffix[4];

		for(std::size_t ipos1=0; ipos1<Npos; ++ipos1)
			this->lri.coefficients[ipos1] =
				[this,ipos1,&get_delta_pos](
					const Label::ab_ab &label, const TA &Aa01, const TAC &Aa2, const TAC &Ab01, const TAC &Ab2) -> Tdata
				{	return get_delta_pos(Aa01, Ab01.first, Ab01.second, ipos1);	};

		this->lri.cal(
			{Label::ab_ab::a0b0_a1b1,
			 Label::ab_ab::a0b0_a2b1},
			dHs_vec);

		this->lri.data_ab_name[Label::ab::a0b0] = "Vs_"+save_names_suffix[1];
		this->lri.data_ab_name[Label::ab::b   ] = "dCs_"+std::to_string(ipos0)+"_"+save_names_suffix[3];

		for(std::size_t ipos1=0; ipos1<Npos; ++ipos1)
			this->lri.coefficients[ipos1] =
				[this,ipos1,&get_delta_pos](
					const Label::ab_ab &label, const TA &Aa01, const TAC &Aa2, const TAC &Ab01, const TAC &Ab2) -> Tdata
				{	return get_delta_pos(Ab01.first, Ab2.first, (Ab2.second-Ab01.second)%this->period, ipos1); };

		this->lri.cal(
			{Label::ab_ab::a0b0_a1b1,
			 Label::ab_ab::a0b0_a2b1},
			dHs_vec);

		for(std::size_t ipos1=0; ipos1<Npos; ++ipos1)
			this->stress(ipos0,ipos1) = post_2D.cal_energy(
				this->post_2D.saves["Ds_"+save_names_suffix[2]],
				this->post_2D.set_tensors_map2(dHs_vec[ipos1]));
	}
}
*/


template<typename TA, typename Tcell, std::size_t Ndim, typename Tdata>
void Exx<TA,Tcell,Ndim,Tdata>::cal_stress(
	const std::array<std::string,5> &save_names_suffix)						// "Cs","Vs","Ds","dCRs","dVRs"
{
	assert(this->flag_finish.stru);
	assert(this->flag_finish.C);
	assert(this->flag_finish.V);
	assert(this->flag_finish.D);
	assert(this->flag_finish.dCR);
	assert(this->flag_finish.dVR);

	this->stress = Tensor<Tdata>({Npos, Npos});

	this->lri.data_ab_name[Label::ab::a1b1] = this->lri.data_ab_name[Label::ab::a2b1] = "Ds_"+save_names_suffix[2];
	for(std::size_t ipos0=0; ipos0<Npos; ++ipos0)
		for(std::size_t ipos1=0; ipos1<Npos; ++ipos1)
		{
			std::vector<std::map<TA,std::map<TAC,Tensor<Tdata>>>> dHs_vec(1);

			this->lri.data_ab_name[Label::ab::a   ] = "dCRs_"+std::to_string(ipos0)+"_"+std::to_string(ipos1)+"_"+save_names_suffix[3];
			this->lri.data_ab_name[Label::ab::a0b0] = "Vs_"+save_names_suffix[1];
			this->lri.data_ab_name[Label::ab::b   ] = "Cs_"+save_names_suffix[0];

			this->lri.cal_loop3(
				{Label::ab_ab::a0b0_a1b1,
				Label::ab_ab::a0b0_a2b1},
				dHs_vec,
				1.0);

			this->lri.data_ab_name[Label::ab::a   ] = "Cs_"+save_names_suffix[0];
			this->lri.data_ab_name[Label::ab::a0b0] = "dVRs_"+std::to_string(ipos0)+"_"+std::to_string(ipos1)+"_"+save_names_suffix[4];

			this->lri.cal_loop3(
				{Label::ab_ab::a0b0_a1b1,
				Label::ab_ab::a0b0_a2b1},
				dHs_vec,
				1.0);

			this->lri.data_ab_name[Label::ab::a0b0] = "Vs_"+save_names_suffix[1];
			this->lri.data_ab_name[Label::ab::b   ] = "dCRs_"+std::to_string(ipos0)+"_"+std::to_string(ipos1)+"_"+save_names_suffix[3];

			this->lri.cal_loop3(
				{Label::ab_ab::a0b0_a1b1,
				Label::ab_ab::a0b0_a2b1},
				dHs_vec,
				1.0);

				this->stress(ipos0,ipos1) = post_2D.cal_energy(
					this->post_2D.saves["Ds_"+save_names_suffix[2]],
					this->post_2D.set_tensors_map2(dHs_vec[0]));
		}
}

}