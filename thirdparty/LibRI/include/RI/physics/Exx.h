// ===================
//  Author: Peize Lin
//  date: 2022.06.02
// ===================

#pragma once

#include "Exx_Post_2D.h"
#include "../global/Global_Func-2.h"
#include "../global/Tensor.h"
#include "../ri/LRI.h"

#include <mpi.h>
#include <array>
#include <map>
#include <set>

namespace RI
{

template<typename TA, typename Tcell, std::size_t Ndim, typename Tdata>
class Exx
{
public:
	using TC = std::array<Tcell,Ndim>;
	using TAC = std::pair<TA,TC>;
	using Tdata_real = Global_Func::To_Real_t<Tdata>;
	using Tpos = double;							// tmp
	constexpr static std::size_t Npos = Ndim;		// tmp
	using Tatom_pos = std::array<Tpos,Npos>;		// tmp

	void set_parallel(
		const MPI_Comm &mpi_comm,
		const std::map<TA,Tatom_pos> &atoms_pos,
		const std::array<Tatom_pos,Ndim> &latvec,
		const std::array<Tcell,Ndim> &period);

	void set_symmetry(
		const bool flag_symmetry,
		const std::map<std::pair<TA,TA>, std::set<TC>> &irreducible_sector);

	void set_Cs(
		const std::map<TA, std::map<TAC, Tensor<Tdata>>> &Cs,
		const Tdata_real &threshold,
		const std::string &save_name_suffix="");
	void set_Vs(
		const std::map<TA, std::map<TAC, Tensor<Tdata>>> &Vs,
		const Tdata_real &threshold,
		const std::string &save_name_suffix="");
	void set_Ds(
		const std::map<TA, std::map<TAC, Tensor<Tdata>>> &Ds,
		const Tdata_real &threshold,
		const std::string &save_name_suffix="");
	void set_Ds_delta(
		const std::map<TA, std::map<TAC, Tensor<Tdata>>> &Ds,
		const Tdata_real &threshold,
		const std::string &save_name_suffix="");
	void set_dCs(
		const std::array<std::map<TA, std::map<TAC, Tensor<Tdata>>>,Npos> &dCs,
		const Tdata_real &threshold,
		const std::string &save_name_suffix="");
	void set_dVs(
		const std::array<std::map<TA, std::map<TAC, Tensor<Tdata>>>,Npos> &dVs,
		const Tdata_real &threshold,
		const std::string &save_name_suffix="");
	void set_dCRs(
		const std::array<std::array<std::map<TA, std::map<TAC, Tensor<Tdata>>>,Npos>,Npos> &dCRs,
		const Tdata_real &threshold,
		const std::string &save_name_suffix="");
	void set_dVRs(
		const std::array<std::array<std::map<TA, std::map<TAC, Tensor<Tdata>>>,Npos>,Npos> &dVRs,
		const Tdata_real &threshold,
		const std::string &save_name_suffix="");

	void cal_Hs(
		const std::array<std::string,3> &save_names_suffix={"","",""});		// "Cs","Vs","Ds"
	void cal_force(
		const std::array<std::string,5> &save_names_suffix={"","","","",""});	// "Cs","Vs","Ds","dCs","dVs"
	void cal_stress(
		const std::array<std::string,5> &save_names_suffix={"","","","",""});	// "Cs","Vs","Ds","dCRs","dVRs"

	std::map<TA, std::map<TAC, Tensor<Tdata>>> Hs;
	std::array<std::array< std::map<TA, std::map<TAC, Tensor<Tdata>>> ,2>,Npos> dHs;
	std::array<std::array< std::map<TA, std::map<TAC, Tensor<Tdata>>> ,Npos>,Npos> dHRs;
	Tdata energy = 0;
	std::array<std::map<TA,Tdata>,Ndim> force;
	Tensor<Tdata> stress = Tensor<Tdata>({Npos, Npos});

	Exx_Post_2D<TA,TC,Tdata> post_2D;

	void free_Cs(const std::string &save_name_suffix="");
	void free_Vs(const std::string &save_name_suffix="");
	void free_Ds(const std::string &save_name_suffix="");
	void free_Ds_delta(const std::string &save_name_suffix="");
	void free_dCs(const std::string &save_name_suffix="");
	void free_dVs(const std::string &save_name_suffix="");
	void free_dCRs(const std::string &save_name_suffix="");
	void free_dVRs(const std::string &save_name_suffix="");

public:
	LRI<TA,Tcell,Ndim,Tdata> lri;

	struct Flag_Finish
	{
		bool stru=false;
		bool Cs=false;
		bool Vs=false;
		bool Ds=false;
		bool Ds_delta=false;
		bool dCs=false;
		bool dVs=false;
		bool dCRs=false;
		bool dVRs=false;
	};
	Flag_Finish flag_finish;

	struct Flag_Save_Result
	{
		bool Hs=true;
		bool dHs=false;
		bool dHRs=false;
	};
	Flag_Save_Result flag_save_result;

};

}

#include "Exx.hpp"