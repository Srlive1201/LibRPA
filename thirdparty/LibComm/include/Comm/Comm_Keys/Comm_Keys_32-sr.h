// ===================
//  Author: Peize Lin
//  date: 2022.06.22
// ===================

#pragma once

#include "../global/Cereal_Func.h"

#include <vector>
#include <functional>
#include <mpi.h>
#include <shared_mutex>

namespace Comm
{

template<typename Tkey, typename Tkeys_provide, typename Tkeys_require>
class Comm_Keys_32
{
public:
	Comm_Keys_32(const MPI_Comm &mpi_comm_in);

	std::vector<std::vector<Tkey>> trans(
		const Tkeys_provide &keys_provide_mine,
		const Tkeys_require &keys_require_mine);

protected:
	void send_keys_require_mine(
		const Tkeys_require &keys_require_mine);

	void recv_require_intersection(
		std::vector<Tkey> &keys_provide_mine,
		std::vector<std::vector<Tkey>> &keys_trans_list);

	void intersection(
		std::vector<Tkey> &keys_provide_mine,
		const Tkeys_require &keys_require,
		std::vector<Tkey> &keys_trans);

	virtual std::vector<Tkey> change_keys_provide_mine(
		const Tkeys_provide &keys_provide_mine)=0;

	MPI_Comm mpi_comm;
	int rank_mine;
	int rank_size;

	const int tag_keys = 88;
	Comm::Cereal_Func cereal_func;
	std::shared_mutex lock_provide;
};

template<typename Tkey, typename Tkeys_provide, typename Tkeys_require>
class Comm_Keys_32_SenderTraversal: public Comm_Keys_32<Tkey, Tkeys_provide, Tkeys_require>
{
public:
	Comm_Keys_32_SenderTraversal(const MPI_Comm &mpi_comm);

	std::function<
		void(
			const Tkeys_provide &keys_provide_mine,
			std::function<void(const Tkey&)> &func )>
		traverse_keys_provide;

private:
	std::vector<Tkey> change_keys_provide_mine(
		const Tkeys_provide &keys_provide_mine);
};

template<typename Tkey, typename Tkeys_provide, typename Tkeys_require>
class Comm_Keys_32_SenderJudge: public Comm_Keys_32<Tkey, Tkeys_provide, Tkeys_require>
{
public:
	Comm_Keys_32_SenderJudge(const MPI_Comm &mpi_comm);

	std::function<
		void(
			std::function<void(const Tkey&)> &func )>
		traverse_keys_all;

private:
	std::vector<Tkey> change_keys_provide_mine(
		const Tkeys_provide &keys_provide_mine);
};

}

#include "Comm_Keys_32-sr.hpp"