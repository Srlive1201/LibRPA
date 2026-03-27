// ==========================================
//  Author: Peize Lin, Rong Shi, Minye Zhang
//  Date:   2022.07.25
// ==========================================

#pragma once

#include "../ri/LRI.h"
#include "../global/Tensor.h"
#include "../global/Global_Func-2.h"

#include <mpi.h>
#include <array>
#include <map>
#include <set>

namespace RI
{

template<typename TA, typename Tcell, std::size_t Ndim, typename Tdata>
class RPA
{
public:
	using TC = std::array<Tcell,Ndim>;
	using TAC = std::pair<TA,TC>;
	using Tdata_real = Global_Func::To_Real_t<Tdata>;
	using Tatom_pos = std::array<double,Ndim>;		// tmp

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

	void set_Gs_pos(
		const std::map<TA, std::map<TAC, Tensor<Tdata>>> &Gs_pos,
		const Tdata_real &threshold,
		const std::string &save_name_suffix="");

	void set_Gs_neg(
		const std::map<TA, std::map<TAC, Tensor<Tdata>>> &Gs_neg,
		const Tdata_real &threshold,
		const std::string &save_name_suffix="");

	void cal_chi0s(
		const std::array<std::string,3> &save_names_suffix={"","",""});						// "Cs","Gs_pos","Gs_neg"

	std::map<TA, std::map<TAC, Tensor<Tdata>>> chi0s;

	void free_Cs(const std::string &save_name_suffix="");
	void free_Gs_pos(const std::string &save_name_suffix="");
	void free_Gs_neg(const std::string &save_name_suffix="");

public:
	LRI<TA,Tcell,Ndim,Tdata> lri;

	struct Flag_Finish
	{
		bool stru=false;
		bool Cs=false;
		bool Gs_pos=false;
		bool Gs_neg=false;
	};
	Flag_Finish flag_finish;
};

}

#include "RPA.hpp"