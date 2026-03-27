// ===================
//  Author: Peize Lin
//  date: 2022.07.15
// ===================

#pragma once

#include "Distribute_Equally.hpp"

#include <vector>
#include <array>
#include <utility>
#include <mpi.h>

namespace RI
{

namespace Distribute_Equally
{
	// num_index 维张量，第0维为atoms，剩余维为{atom,cell}。
	// 返回值为每维上本进程被分配到的值列表，first为第0维，second[i-1]为第i维。
	// 当 任务数<进程数 时，部分进程会被分配到重复任务。if(!flag_task_repeatable)，重复进程只保留一个，其他进程将任务删空。

	// 全部维按照atoms，尽可能均分
	template<typename TA, typename Tcell, std::size_t Ndim>
	extern std::pair<std::vector<TA>,
	                 std::vector<std::vector<std::pair<TA,std::array<Tcell,Ndim>>>>>
	distribute_atoms(
		const MPI_Comm &mpi_comm,
		const std::vector<TA> &atoms,
		const std::array<Tcell,Ndim> &period,
		const std::size_t num_index,
		const bool flag_task_repeatable);

	// 第0维按照atoms、剩余维按照{atom,period}，尽可能均分
	template<typename TA, typename Tcell, std::size_t Ndim>
	extern std::pair<std::vector<TA>,
	                 std::vector<std::vector<std::pair<TA,std::array<Tcell,Ndim>>>>>
	distribute_atoms_periods(
		const MPI_Comm &mpi_comm,
		const std::vector<TA> &atoms,
		const std::array<Tcell,Ndim> &period,
		const std::size_t num_index,
		const bool flag_task_repeatable);



	// num_index 维张量，每维皆为{atom,cell}。
	// 返回值为每维上本进程被分配到的值列表，[i]为第i维。
	// 当 任务数<进程数 时，部分进程会被分配到重复任务。if(!flag_task_repeatable)，重复进程只保留一个，其他进程将任务删空。

	// 全部维按照{atom,period}，尽可能均分
	template<typename TA, typename Tcell, std::size_t Ndim>
	extern std::vector<std::vector<std::pair<TA,std::array<Tcell,Ndim>>>>
	distribute_periods(
		const MPI_Comm &mpi_comm,
		const std::vector<TA> &atoms,
		const std::array<Tcell,Ndim> &period,
		const std::size_t num_index,
		const bool flag_task_repeatable);
}

}