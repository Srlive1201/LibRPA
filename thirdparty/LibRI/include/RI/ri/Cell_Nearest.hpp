// ===================
//  Author: Peize Lin
//  date: 2022.12.30
// ===================

#pragma once

#include "Cell_Nearest.h"

#include "../global/Tensor.h"
#include "../global/Array_Operator.h"
#include "../global/Tensor_Algorithm.h"

#include <cmath>

namespace RI
{

template<typename TA, typename Tcell, std::size_t Ndim, typename Tpos, std::size_t Npos>
void Cell_Nearest<TA,Tcell,Ndim,Tpos,Npos>::init(
	const std::map<TA,Tatom_pos> &atoms_pos,
	const std::array<Tatom_pos,Ndim> &latvec_in,
	const std::array<Tcell,Ndim> &period_in)
{
	using namespace Array_Operator;
	this->period = period_in;

	Tensor<Tpos> latvec({Ndim, Npos});
	for(std::size_t idim=0; idim<Ndim; ++idim)
		for(std::size_t ipos=0; ipos<Npos; ++ipos)
			latvec(idim,ipos) = latvec_in[idim][ipos];
	const Tensor<Tpos> least_square_tmp				// shape:{Ndim,Npos}
		= - Tensor_Algorithm::inverse_matrix_heev(latvec * latvec.transpose()) * latvec;

	for(const auto &atoms_pos_x : atoms_pos)
	{
		const TA &Ax = atoms_pos_x.first;
		for(const auto &atoms_pos_y : atoms_pos)
		{
			const TA &Ay = atoms_pos_y.first;
			const Tensor<Tpos> delta_pos = to_Tensor(atoms_pos_y.second - atoms_pos_x.second);
			this->cells_nearest_continuous[Ax][Ay] = to_array<Tpos,Ndim>(least_square_tmp * delta_pos);
		}
	}
}

template<typename TA, typename Tcell, std::size_t Ndim, typename Tpos, std::size_t Npos>
auto Cell_Nearest<TA,Tcell,Ndim,Tpos,Npos>::get_cell_nearest_discrete(
	const TA &Ax, const TA &Ay, const TC &cell) const
-> TC
{
	const std::array<Tpos,Ndim> &cell_nearest_continuous = this->cells_nearest_continuous.at(Ax).at(Ay);
	TC cell_nearest_discrete;
	for(std::size_t idim=0; idim<Ndim; ++idim)
		cell_nearest_discrete[idim]
			=  std::round( (cell_nearest_continuous[idim]-cell[idim]) / this->period[idim] )
				* this->period[idim]
				+ cell[idim];
	return cell_nearest_discrete;
}

}
