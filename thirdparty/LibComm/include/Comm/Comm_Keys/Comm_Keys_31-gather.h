// ===================
//  Author: Peize Lin
//  date: 2022.06.22
// ===================

#pragma once

#include <vector>
#include <functional>
#include <mpi.h>

namespace Comm
{

template<typename Tkey, typename Tkeys_provide, typename Tkeys_require>
class Comm_Keys_31
{
public:
	Comm_Keys_31(const MPI_Comm &mpi_comm_in);

	std::vector<std::vector<Tkey>> trans(
		const Tkeys_provide &keys_provide_mine,
		const Tkeys_require &keys_require_mine);

protected:
#if MPI_VERSION>=4
	void send_keys_require_mine(
		const Tkeys_require &keys_require_mine,
		std::vector<MPI_Count> &sss_size,
		std::vector<MPI_Aint> &sss_displs,
		std::vector<char> &buffer_recv);
#else
	void send_keys_require_mine(
		const Tkeys_require &keys_require_mine,
		std::vector<int> &sss_size,
		std::vector<int> &sss_displs,
		std::vector<char> &buffer_recv);
#endif

	//void recv_require_intersection(
	//	const Tkeys_provide &keys_provide_mine,
	//	std::vector<std::vector<Tkey>> &keys_trans_list);

	virtual void intersection(
		const Tkeys_provide &keys_provide_mine,
		const Tkeys_require &keys_require,
		std::vector<Tkey> &keys_trans)=0;

	MPI_Comm mpi_comm;
	int rank_mine;
	int rank_size;

	const int tag_keys = 99;
};

template<typename Tkey, typename Tkeys_provide, typename Tkeys_require>
class Comm_Keys_31_SenderTraversal: public Comm_Keys_31<Tkey, Tkeys_provide, Tkeys_require>
{
public:
	Comm_Keys_31_SenderTraversal(const MPI_Comm &mpi_comm);

	std::function<
		void(
			const Tkeys_provide &keys_provide_mine,
			std::function<void(const Tkey&)> &func )>
		traverse_keys_provide;

private:
	void intersection(
		const Tkeys_provide &keys_provide_mine,
		const Tkeys_require &keys_require,
		std::vector<Tkey> &keys_trans);
};

template<typename Tkey, typename Tkeys_provide, typename Tkeys_require>
class Comm_Keys_31_SenderJudge: public Comm_Keys_31<Tkey, Tkeys_provide, Tkeys_require>
{
public:
	Comm_Keys_31_SenderJudge(const MPI_Comm &mpi_comm);

	std::function<
		void(
			std::function<void(const Tkey&)> &func )>
		traverse_keys_all;

private:
	void intersection(
		const Tkeys_provide &keys_provide_mine,
		const Tkeys_require &keys_require,
		std::vector<Tkey> &keys_trans);
};

}

#include "Comm_Keys_31-gather.hpp"