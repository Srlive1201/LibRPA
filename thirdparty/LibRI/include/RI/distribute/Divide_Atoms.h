// ===================
//  Author: Peize Lin
//  date: 2022.07.13
// ===================

#pragma once

#include <vector>
#include <array>
#include <utility>

namespace RI
{

namespace Divide_Atoms
{
	// equally divide atoms:
	// 	[0,1,2]  [3,4,5]  [6,7]
	template<typename TA>
	extern std::vector<TA> divide_atoms(
		const std::size_t group_rank,
		const std::size_t group_size,
		const std::vector<TA> &atoms);

	// equally divide atoms:
	// 	[0,1]  [2]
	// with all period (period=2)
	// 	[{0,0},{0,1},{1,0},{1,1}]  [{2,0},{2,1}]
	template<typename TA, typename Tcell, std::size_t Ndim>
	extern std::vector<std::pair<TA,std::array<Tcell,Ndim>>> divide_atoms(
		const std::size_t group_rank,
		const std::size_t group_size,
		const std::vector<TA> &atoms,
		const std::array<Tcell,Ndim> &period);

	// equally divide atoms and periods (period=2):
	// 	[{0,0},{0,1},{1,0}]  [{1,1},{2,0},{2,1}]  [{3,0},{3,1}]
	template<typename TA, typename Tcell, std::size_t Ndim>
	std::vector<std::pair<TA,std::array<Tcell,Ndim>>> divide_atoms_periods(
		const std::size_t group_rank,
		const std::size_t group_size,
		const std::vector<TA> &atoms,
		const std::array<Tcell,Ndim> &period);
}

}

#include "Divide_Atoms.hpp"