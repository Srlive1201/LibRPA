// ===================
//  Author: Peize Lin
//  date: 2022.07.23
// ===================

#pragma once

#include "Parallel_LRI.h"

namespace RI
{

template<typename TA, typename Tcell, std::size_t Ndim, typename Tdata>
class Parallel_LRI_Equally: public Parallel_LRI<TA,Tcell,Ndim,Tdata>
{
  public:
	using TC = std::array<Tcell,Ndim>;
	using TAC = std::pair<TA,TC>;
	using Tatom_pos = std::array<double,Ndim>;		// tmp

	void set_parallel(
		const MPI_Comm &mpi_comm_in,
		const std::map<TA,Tatom_pos> &atoms_pos,
		const std::array<Tatom_pos,Ndim> &latvec,
		const std::array<Tcell,Ndim> &period_in,
		const std::set<Label::Aab_Aab> &labels) override;

//	std::map<TA,std::map<TAC,Tensor<Tdata>>> comm_tensors_map2(
//		const Label::ab &label,
//		const std::map<TA,std::map<TAC,Tensor<Tdata>>> &Ds) const override;
	std::map<TA,std::map<TAC,Tensor<Tdata>>> comm_tensors_map2(
		const std::vector<Label::ab> &label,
		const std::map<TA,std::map<TAC,Tensor<Tdata>>> &Ds) const override;

	const std::vector<TA >& get_list_Aa01() const override { return this->list_Aa01; }
	const std::vector<TAC>& get_list_Aa2 (const TA &Aa01) const override { return this->list_Aa2;  }
	const std::vector<TAC>& get_list_Ab01(const TA &Aa01, const TAC &Aa2) const override { return this->list_Ab01; }
	const std::vector<TAC>& get_list_Ab2 (const TA &Aa01, const TAC &Aa2, const TAC &Ab01) const override { return this->list_Ab2;  }

  public:	// private:
	std::vector<TA>  list_Aa01;
	std::vector<TAC> list_Aa2;
	std::vector<TAC> list_Ab01;
	std::vector<TAC> list_Ab2;

	MPI_Comm mpi_comm;
	std::array<Tcell,Ndim> period;

  public:	// private:
	void set_parallel_loop4(
		const std::vector<TA> &atoms_vec);
	void set_parallel_loop3(
		const std::vector<TA> &atoms_vec,
		const std::set<Label::Aab_Aab> &labels);
};

}

#include "Parallel_LRI_Equally.hpp"