// ===================
//  Author: Peize Lin
//  date: 2022.12.30
// ===================

#pragma once

#include <array>
#include <map>

namespace RI
{

template<typename TA, typename Tcell, std::size_t Ndim, typename Tpos, std::size_t Npos>
class Cell_Nearest
{
public:
	using Tatom_pos = std::array<Tpos,Npos>;		// tmp
	using TC = std::array<Tcell,Ndim>;

	void init(
		const std::map<TA,Tatom_pos> &atoms_pos,
		const std::array<Tatom_pos,Ndim> &latvec_in,
		const std::array<Tcell,Ndim> &period_in);

	TC get_cell_nearest_discrete(const TA &Ax, const TA &Ay, const TC &cell) const;

public:		//private:
	TC period;
	std::map<TA,std::map<TA,std::array<Tpos,Ndim>>> cells_nearest_continuous;
};

}

#include "Cell_Nearest.hpp"