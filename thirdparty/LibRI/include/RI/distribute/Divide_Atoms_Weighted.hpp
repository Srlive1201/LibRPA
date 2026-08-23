// ===================
//  Author: maki49
//  date: 2026.07.10
// ===================

#pragma once

#include "Divide_Atoms_Weighted.h"
#include "Divide_Atoms.h"
#include "../global/Global_Func-3.h"

#include <numeric>
#include <stdexcept>
#include <string>
#include <algorithm>
#include <functional>
#include <queue>
#include <cassert>

namespace RI
{

namespace Divide_Atoms_Weighted
{
	template<typename TA>
	std::size_t atom_weight(
		const std::map<TA,std::size_t> &weights,
		const TA &atom)
	{
		const auto ptr = weights.find(atom);
		return (ptr==weights.end()) ? 1 : ptr->second;
	}

	// Core distribution function. 
	// Assign the heaviest item to the currently lightest group.
	inline std::vector<std::size_t> assign_owners_weighted(
		const std::size_t group_size,
		const std::vector<std::size_t> &item_weights)
	{
		assert(group_size>0);

		std::vector<std::size_t> order(item_weights.size());
		std::iota(order.begin(), order.end(), std::size_t(0));
		std::stable_sort(order.begin(), order.end(),
			[&item_weights](const std::size_t i, const std::size_t j) -> bool
			{ return item_weights[i] > item_weights[j]; });

		using Tload_group = std::pair<std::size_t,std::size_t>;	// weight, atom-index
		std::priority_queue<Tload_group, std::vector<Tload_group>, std::greater<Tload_group>> groups;	//min heap
		for(std::size_t group=0; group<group_size; ++group)
			groups.push(std::make_pair(std::size_t(0), group));

		std::vector<std::size_t> owners(item_weights.size());
		for(const std::size_t i : order)
		{
			const Tload_group lightest = groups.top();
			groups.pop();
			owners[i] = lightest.second;
			groups.push(std::make_pair(lightest.first+item_weights[i], lightest.second));
		}
		return owners;
	}

	template<typename TA>
	std::vector<TA> divide_atoms(
		const std::size_t group_rank,
		const std::size_t group_size,
		const std::vector<TA> &atoms,
		const std::map<TA,std::size_t> &weights)
	{
		if(weights.empty())
			return Divide_Atoms::divide_atoms(group_rank, group_size, atoms);

		std::vector<std::size_t> item_weights;
		item_weights.reserve(atoms.size());
		for(const TA &atom : atoms)
			item_weights.push_back(atom_weight(weights, atom));

		const std::vector<std::size_t> owners = assign_owners_weighted(group_size, item_weights);

		// get the atoms distributed onto the current rank
		std::vector<TA> atoms_divide;
		for(std::size_t i=0; i<atoms.size(); ++i)
			if(owners[i]==group_rank)
				atoms_divide.push_back(atoms[i]);
		return atoms_divide;
	}

	template<typename TA, typename Tcell, std::size_t Ndim>
	std::vector<std::pair<TA,std::array<Tcell,Ndim>>> divide_atoms(
		const std::size_t group_rank,
		const std::size_t group_size,
		const std::vector<TA> &atoms,
		const std::array<Tcell,Ndim> &period,
		const std::map<TA,std::size_t> &weights)
	{
		const std::vector<TA> atoms_divide = divide_atoms(group_rank, group_size, atoms, weights);
		return Divide_Atoms::traversal_atom_period(atoms_divide, period);
	}

	template<typename TA, typename Tcell, std::size_t Ndim>
	std::vector<std::pair<TA,std::array<Tcell,Ndim>>> divide_atoms_periods(
		const std::size_t group_rank,
		const std::size_t group_size,
		const std::vector<TA> &atoms,
		const std::array<Tcell,Ndim> &period,
		const std::map<TA,std::size_t> &weights)
	{
		if(weights.empty())
			return Divide_Atoms::divide_atoms_periods(group_rank, group_size, atoms, period);

		using TC = std::array<Tcell,Ndim>;
		using TAC = std::pair<TA,TC>;
		const std::vector<TC> cells_origin = Divide_Atoms::traversal_period(period);
		const std::vector<TC> cells = Global_Func::mod_period(cells_origin, period);

		std::vector<TAC> atoms_periods;
		std::vector<std::size_t> item_weights;
		atoms_periods.reserve( atoms.size() * cells.size() );
		item_weights.reserve( atoms.size() * cells.size() );
		for(const TA &atom : atoms)
		{
			const std::size_t weight = atom_weight(weights, atom);
			for(const TC &cell : cells)
			{
				atoms_periods.push_back(std::make_pair(atom,cell));
				item_weights.push_back(weight);
			}
		}

		const std::vector<std::size_t> owners = assign_owners_weighted(group_size, item_weights);

		std::vector<TAC> atoms_periods_divide;
		for(std::size_t i=0; i<atoms_periods.size(); ++i)
			if(owners[i]==group_rank)
				atoms_periods_divide.push_back(atoms_periods[i]);
		return atoms_periods_divide;
	}
}

}