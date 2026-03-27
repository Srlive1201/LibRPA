// ===================
//  Author: Peize Lin
//  date: 2022.06.02
// ===================

#pragma once

#include "LRI_Cal_Tools.h"
#include "Label.h"
#include "../global/Tensor.h"
#include "Data_Pack.h"
#include "../parallel/Parallel_LRI_Equally.h"
#include "RI_Tools.h"
#include "../global/Global_Func-2.h"
#include "Filter_Atom.h"

#include <mpi.h>
#include <array>
#include <vector>
#include <map>
#include <set>
#include <unordered_map>
#include <memory>
#include <functional>

namespace RI
{

template<typename TA, typename Tcell, std::size_t Ndim, typename Tdata>
class LRI
{
public:
	using TC = std::array<Tcell,Ndim>;
	using TAC = std::pair<TA,TC>;
	using Tdata_real = Global_Func::To_Real_t<Tdata>;
	using Tatom_pos = std::array<double,Ndim>;		// tmp

	LRI();

	//template<typename Tatom_pos>
	void set_parallel(
		const MPI_Comm &mpi_comm_in,
		const std::map<TA,Tatom_pos> &atoms_pos,
		const std::array<Tatom_pos,Ndim> &latvec,
		const std::array<Tcell,Ndim> &period_in,
		const std::vector<Label::ab_ab> &labels_all_in);

	void set_tensors_map2(
		const std::map<TA, std::map<TAC, Tensor<Tdata>>> &Ds_local,
		const std::vector<Label::ab> &label_list,
		const std::map<std::string, double> &para_in = {},
		const std::string &save_name_in = "default");
			// para:
			//     "flag_period",      true
			//     "flag_comm",        true
			//     "flag_filter",      true
			//     "threshold_filter", 0.0
			// save_name:              Label_Tools::get_name(label)
	void free_tensors_map2(
		const std::string &save_name);

	void cal_loop3(
		const std::vector<Label::ab_ab> &labels,
		std::map<TA, std::map<TAC, Tensor<Tdata>>> &Ds_result,
		const double fac_add_Ds = 1.0);

public:
	std::shared_ptr<Parallel_LRI<TA,Tcell,Ndim,Tdata>>
		parallel = std::make_shared<Parallel_LRI_Equally<TA,Tcell,Ndim,Tdata>>();
	std::unordered_map< Label::ab, RI_Tools::T_filter_func<Tdata> >
		filter_funcs;
	std::shared_ptr<Filter_Atom<TA,TAC>> filter_atom = std::make_shared<Filter_Atom<TA,TAC>>();

public:		// private:
	TC period;
	MPI_Comm mpi_comm;
	std::map<std::string, Data_Pack<TA,TC,Tdata>> data_pool;
	std::unordered_map<Label::ab, std::string> data_ab_name;

public:		// private:
	using T_cal_func = std::function<void(
		const Label::ab_ab &label,
		const TA &Aa01, const TAC &Aa2, const TAC &Ab01, const TAC &Ab2,
		const Tensor<Tdata> &D_a, const Tensor<Tdata> &D_b,
		const Tensor<Tdata> &D_a_transpose, const Tensor<Tdata> &D_b_transpose,
		std::unordered_map<Label::ab_ab, Tensor<Tdata>> &Ds_b01,
		std::unordered_map<Label::ab_ab, Tdata_real> &Ds_b01_csm,
		LRI_Cal_Tools<TA,TC,Tdata> &tools)>;
	std::unordered_map<Label::ab_ab, T_cal_func> cal_funcs;
};

}

#include "LRI.hpp"
#include "LRI-set.hpp"
#include "LRI-cal_loop3.hpp"
