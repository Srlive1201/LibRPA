// ===================
//  Author: Peize Lin
//  date: 2022.07.11
// ===================

#pragma once

#include "Split_Processes.h"
#include "../global/MPI_Wrapper.h"

#include <numeric>
#include <cassert>
#include <string>
#include <stdexcept>

#define MPI_CHECK(x) if((x)!=MPI_SUCCESS)	throw std::runtime_error(std::string(__FILE__)+" line "+std::to_string(__LINE__));

namespace RI
{

namespace Split_Processes
{
	// comm_color
	static std::tuple<MPI_Wrapper::mpi_comm,std::size_t>
	split(
		const MPI_Comm &mc,
		const std::size_t &group_size)
	{
		assert(group_size>0);
		const std::size_t rank_mine = static_cast<std::size_t>(MPI_Wrapper::mpi_get_rank(mc));
		const std::size_t rank_size = static_cast<std::size_t>(MPI_Wrapper::mpi_get_size(mc));
		assert(rank_size>=group_size);

		std::vector<std::size_t> num(group_size);			// sum(num) = rank_size
		const std::size_t mod = rank_size % group_size;
		for(std::size_t i=0; i<mod; ++i)
			num[i] = rank_size/group_size+1;
		for(std::size_t i=mod; i<group_size; ++i)
			num[i] = rank_size/group_size;

		std::vector<std::size_t> slice(group_size);			// slice.back() = rank_size
		slice[0] = num[0];
		for(std::size_t i=1; i<group_size; ++i)
			slice[i] = slice[i-1] + num[i];

		const std::size_t color_group = [&]() -> std::size_t		// which group should rank_mine in
		{
			for(std::size_t i=0; i<group_size; ++i)
				if(rank_mine < slice[i])
					return i;
			throw std::range_error(std::string(__FILE__)+" line "+std::to_string(__LINE__));
		}();

		MPI_Wrapper::mpi_comm mc_split;
		MPI_CHECK( MPI_Comm_split(
			mc, static_cast<int>(color_group), static_cast<int>(rank_mine), &mc_split() ) );
		mc_split.flag_allocate = true;

		return std::forward_as_tuple(std::move(mc_split), color_group);
	}

	// comm_color_size
	static std::tuple<MPI_Wrapper::mpi_comm, std::size_t, std::size_t>
	split_first(
		const MPI_Comm &mc,
		const std::vector<std::size_t> &task_sizes)
	{
		assert(task_sizes.size()>=1);
		const std::size_t rank_size = static_cast<std::size_t>(MPI_Wrapper::mpi_get_size(mc));
		const std::size_t task_product = std::accumulate(
			task_sizes.begin(), task_sizes.end(), std::size_t(1), std::multiplies<std::size_t>() );		// double for numerical range
		const double num_average =
			task_product < rank_size
			? 1.0		// if task_product<rank_size, then num_average<1, then group_size>task_sizes[0]. Set group_size=task_sizes[0]
			: std::pow(static_cast<double>(task_product)/rank_size, 1.0/task_sizes.size());
		const std::size_t group_size =
			task_sizes[0] < num_average
			? 1			// if task_sizes[0]<<task_sizes[1:], then group_size<0.5. Set group_size=1
			: static_cast<std::size_t>(std::round(task_sizes[0]/num_average));
		std::tuple<MPI_Wrapper::mpi_comm, std::size_t>
			comm_color = split(mc, group_size);
		return std::make_tuple(std::move(std::get<0>(comm_color)), std::get<1>(comm_color), group_size);
	}

	// vector<comm_color_size>
	static std::vector<std::tuple<MPI_Wrapper::mpi_comm, std::size_t, std::size_t>>
	split_all(
		const MPI_Comm &mc,
		const std::vector<std::size_t> &task_sizes)
	{
		std::vector<std::tuple<MPI_Wrapper::mpi_comm, std::size_t,std::size_t>>
			comm_color_sizes(task_sizes.size()+1);
		comm_color_sizes[0] = std::forward_as_tuple(
			MPI_Wrapper::mpi_comm(mc,false),
			0,
			1);
		for(std::size_t m=0; m<task_sizes.size(); ++m)
		{
			comm_color_sizes[m+1] = split_first(
				std::get<0>(comm_color_sizes[m])(),
				{task_sizes.begin()+m, task_sizes.end()});
		}
		return comm_color_sizes;
	}
}

}

#undef MPI_CHECK