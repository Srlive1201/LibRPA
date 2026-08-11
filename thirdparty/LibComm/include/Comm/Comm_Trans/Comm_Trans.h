//=======================
// AUTHOR : Peize Lin
// DATE :   2022-01-05
//=======================

#pragma once

#include "../Comm_Tools.h"
#include "../global/Cereal_Func.h"
#include "State.h"
#include <mpi.h>
#include <vector>
#include <functional>
#include <future>
#include <sstream>

namespace Comm
{
	class Memory_Check;

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
	void communicate_nonblocking(
		const Tdatas_isend &datas_isend,
		Tdatas_recv &datas_recv);
	void communicate_sendrecv_ring(
		const Tdatas_isend &datas_isend,
		Tdatas_recv &datas_recv);
	std::vector<char> serialize_for_rank(
		const int rank_isend,
		const Tdatas_isend &datas_isend) const;
	void deserialize_from_rank(
		const int rank_recv,
		const std::vector<char> &buffer_recv,
		Tdatas_recv &datas_recv) const;
	std::size_t oar_data(
		const int rank_isend,
		const Tdatas_isend &datas_isend,
		std::string &buffer_isend,
		std::atomic<State_Send> &state_send,
		Memory_Check &memory);
	void isend_data(
		const int rank_isend,
		const std::size_t exponent_align,
		std::string &buffer_isend,
		MPI_Request &request_isend,
		std::atomic<State_Send> &state_send);
	void recv_data (
		const MPI_Status status_recv,
		const MPI_Message message_recv,
		Memory_Check &memory,
		std::vector<char> &buffer_recv,
		std::atomic<State_Recv> &state_recv);
	void iar_data (
		const int rank_recv,
		std::vector<char> &buffer_recv,
		std::atomic_flag &lock_set_value,
		Tdatas_recv &datas_recv,
		std::atomic<State_Recv> &state_recv) const;
	bool check_finish(
		const std::vector<std::atomic<State_Send>> &states_send,
		const std::vector<std::atomic<State_Recv>> &states_recv) const;

public:
	MPI_Comm mpi_comm;
	int rank_mine = 0;
	int comm_size = 1;

private:
	const int tag_data = 0;
	Comm::Cereal_Func cereal_func;
};

}

#include "Comm_Trans.hpp"
