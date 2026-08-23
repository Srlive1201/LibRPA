// ===================
//  Author: maki49
//  date: 2026.07.10
// ===================

#pragma once

#include <vector>
#include <array>
#include <utility>
#include <map>
#include <cstddef>

namespace RI
{

namespace Divide_Atoms_Weighted
{
	// weights[atom] is the block size of atom (e.g. its number of atomic orbitals).
	// A missing or zero weight is clamped to 1.
	template<typename TA>
	extern std::size_t atom_weight(
		const std::map<TA,std::size_t> &weights,
		const TA &atom);

	// Greedy LPT (longest-processing-time): walk the items from heaviest to lightest,
	// giving each to the currently lightest group. Returns owners[i] = the group of item i.
	// max(load) is within 4/3 of the optimum.
	// Every process calls this with identical arguments and no communication, so the result
	// must be a pure function of them: ties are broken by the smallest original index
	// (stable_sort) and then by the smallest group index (the {load,group} min-heap).
	// group_size > item_weights.size() leaves the surplus groups empty, as the unweighted
	// divide_atoms below already does.
	inline std::vector<std::size_t> assign_owners_weighted(
		const std::size_t group_size,
		const std::vector<std::size_t> &item_weights);

	// divide atoms balancing sum(weights) instead of the atom count.
	// weights.empty() falls back to the unweighted overload above.
	// The returned atoms keep their order in `atoms`, but are no longer a contiguous slice of it.
	template<typename TA>
	extern std::vector<TA> divide_atoms(
		const std::size_t group_rank,
		const std::size_t group_size,
		const std::vector<TA> &atoms,
		const std::map<TA,std::size_t> &weights);

	// divide atoms balancing sum(weights), then expand every atom over all periods.
	// An atom's periods all stay in the same group, so weighting by weights[atom] and by
	// weights[atom]*n_period are the same partition.
	template<typename TA, typename Tcell, std::size_t Ndim>
	extern std::vector<std::pair<TA,std::array<Tcell,Ndim>>> divide_atoms(
		const std::size_t group_rank,
		const std::size_t group_size,
		const std::vector<TA> &atoms,
		const std::array<Tcell,Ndim> &period,
		const std::map<TA,std::size_t> &weights);

	// divide {atom,period} pairs balancing sum(weights) instead of the pair count.
	// Each {atom,period} is an item of weight weights[atom], so an atom's periods may end up
	// in different groups — as they already may in the unweighted overload above.
	template<typename TA, typename Tcell, std::size_t Ndim>
	std::vector<std::pair<TA,std::array<Tcell,Ndim>>> divide_atoms_periods(
		const std::size_t group_rank,
		const std::size_t group_size,
		const std::vector<TA> &atoms,
		const std::array<Tcell,Ndim> &period,
		const std::map<TA,std::size_t> &weights);
}

}

#include "Divide_Atoms_Weighted.hpp"