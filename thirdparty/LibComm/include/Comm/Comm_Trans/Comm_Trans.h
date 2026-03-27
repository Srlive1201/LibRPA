//=======================
// AUTHOR : Peize Lin
// DATE :   2022-01-05
//=======================

#pragma once

#include "../Comm_Tools.h"
#include "../global/Global_Func.h"
#include <mpi.h>
#include <functional>
#include <future>
#include <sstream>

namespace Comm
{

template<typename Tkey, typename Tvalue, typename Tdatas_isend, typename Tdatas_recv>
class Comm_Trans
{
public:
	std::function<
		void(
			const Tdatas_isend &datas_isend,
			const int rank_isend,
			std::function<void(const Tkey&, const Tvalue&)> &func)>
		traverse_isend;
	std::function<
		void(
			Tkey &&key,
			Tvalue &&value,
			Tdatas_recv &datas_recv)>
		set_value_recv;

	Comm_Tools::Lock_Type flag_lock_set_value;
	std::function<
		Tdatas_recv(
			const int rank_recv)>
		init_datas_local;
	std::function<
		void(
			Tdatas_recv &&datas_local,
			Tdatas_recv &datas_recv)>
		add_datas;

public:
	Comm_Trans(const MPI_Comm &mpi_comm_in);
//	Comm_Trans(const Comm_Trans &com);
	void communicate(
		const Tdatas_isend &datas_isend,
		Tdatas_recv &datas_recv);

private:
	void isend_data (const int rank_isend, const Tdatas_isend &datas_isend, std::string &str_isend, MPI_Request &request_isend, std::atomic<std::size_t> &memory_max_isend);
	void recv_data (Tdatas_recv &datas_recv, const MPI_Status status_recv, MPI_Message message_recv, std::atomic_flag &lock_set_value, std::atomic<std::size_t> &memory_max_isend);
	void post_process(
		std::vector<MPI_Request> &requests_isend,
		std::vector<std::string> &strs_isend,
		std::vector<std::future<void>> &futures_isend,
		std::vector<std::future<void>> &futures_recv) const;
	bool memory_enough(const std::atomic<size_t> &memory_max) const { return Global_Func::memory_available() > memory_max.load() * 2; }

public:
	const MPI_Comm &mpi_comm;
	int rank_mine = 0;
	int comm_size = 1;

private:
	const int tag_data = 0;
	Comm::Cereal_Func cereal_func;
};

}

#include "Comm_Trans.hpp"