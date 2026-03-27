// ===================
//  Author: Peize Lin
//  date: 2022.07.11
// ===================

#pragma once

#include "../global/MPI_Wrapper.h"

#include <mpi.h>
#include <tuple>
#include <vector>

namespace RI
{

namespace Split_Processes
{
	// 将所有进程按划分为 group_size 组，每组进程数尽量相同
	// 返回 {本进程所在组mpi_comm, 本进程属于第几组}
	static std::tuple<MPI_Wrapper::mpi_comm, std::size_t> split(
		const MPI_Comm &mc,
		const std::size_t &group_size);

	// 任务数多维，所有进程多维划分组，每维分得任务数尽量相同。
	// 返回按第0维划分结果，{本进程所在组mpi_comm, 本进程属于第几组, 总组数}
	static std::tuple<MPI_Wrapper::mpi_comm, std::size_t, std::size_t> split_first(
		const MPI_Comm &mc,
		const std::vector<std::size_t> &task_sizes);

	// 任务数多维，所有进程多维划分，每维分得任务数尽量相同。
	// 返回按所有维划分结果，返回[0]为 {mc,0,1}，返回[i+1]为 按第i维划分结果
	static std::vector<std::tuple<MPI_Wrapper::mpi_comm, std::size_t, std::size_t>> split_all(
		const MPI_Comm &mc,
		const std::vector<std::size_t> &task_sizes);
}

}

#include "Split_Processes.hpp"