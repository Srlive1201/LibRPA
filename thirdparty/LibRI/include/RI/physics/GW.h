// ===================
//  Author: Minye Zhang, almost completely copied from Exx.h
//  date: 2022.12.18
// ===================

#pragma once

// #include "Exx_Post_2D.h"
#include "../global/Global_Func-2.h"
#include "../global/Tensor.h"
#include "../ri/LRI.h"

#include <mpi.h>
#include <array>
#include <map>
#include <set>

namespace RI
{

//! class to compute the correlation self-energy in GW approximation by space-time method
template<typename TA, typename Tcell, std::size_t Ndim, typename Tdata>
class GW
{
public:
	using TC = std::array<Tcell,Ndim>;
	using TAC = std::pair<TA,TC>;
	using Tdata_real = Global_Func::To_Real_t<Tdata>;
	constexpr static std::size_t Npos = Ndim;		// tmp
	using Tatom_pos = std::array<double,Npos>;		// tmp

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

	// setup screened Coulomb interaction
	void set_Ws(
		const std::map<TA, std::map<TAC, Tensor<Tdata>>> &Ws,
		const Tdata_real &threshold,
		const std::string &save_name_suffix="");

	// setup Green's function
	void set_Gs(
		const std::map<TA, std::map<TAC, Tensor<Tdata>>> &Gs,
		const Tdata_real &threshold,
		const std::string &save_name_suffix="");

	void cal_Sigmas(
		const std::array<std::string,3> &save_names_suffix={"","",""});		// "Cs","Ws","Gs"

	std::map<TA, std::map<TAC, Tensor<Tdata>>> Sigmas;

	void free_Cs(const std::string &save_name_suffix="");
	void free_Ws(const std::string &save_name_suffix="");
	void free_Gs(const std::string &save_name_suffix="");

public:
	LRI<TA,Tcell,Ndim,Tdata> lri;

	struct Flag_Finish
	{
		bool stru=false;
		bool Cs=false;
		bool Ws=false;
		bool Gs=false;
	};
	Flag_Finish flag_finish;
};

}

#include "GW.hpp"