// ===================
//  Author: Peize Lin
//  date: 2022.06.22
// ===================

#pragma once

#include "Comm_Keys_31-sr.h"

#include <mpi.h>
#include <thread>
#include <sstream>
#include <cereal/archives/binary.hpp>
#include <string>
#include <cassert>

#define MPI_CHECK(x) if((x)!=MPI_SUCCESS)	throw std::runtime_error(std::string(__FILE__)+" line "+std::to_string(__LINE__));

namespace Comm
{

template<typename Tkey, typename Tkeys_provide, typename Tkeys_require>
Comm_Keys_31<Tkey,Tkeys_provide,Tkeys_require>::Comm_Keys_31(
	const MPI_Comm &mpi_comm_in)
	:mpi_comm(mpi_comm_in)
{
	MPI_CHECK( MPI_Comm_size( this->mpi_comm, &this->rank_size ) );
	MPI_CHECK( MPI_Comm_rank( this->mpi_comm, &this->rank_mine ) );
}

template<typename Tkey, typename Tkeys_provide, typename Tkeys_require>
Comm_Keys_31_SenderTraversal<Tkey,Tkeys_provide,Tkeys_require>::Comm_Keys_31_SenderTraversal(
	const MPI_Comm &mpi_comm)
	:Comm_Keys_31<Tkey,Tkeys_provide,Tkeys_require>(mpi_comm)
{
	this->traverse_keys_provide =
		[](const Tkeys_provide &keys_provide, std::function<void(const Tkey&)> &func)
		{ throw std::logic_error("Function traverse not set."); };
}

template<typename Tkey, typename Tkeys_provide, typename Tkeys_require>
Comm_Keys_31_SenderJudge<Tkey,Tkeys_provide,Tkeys_require>::Comm_Keys_31_SenderJudge(
	const MPI_Comm &mpi_comm)
	:Comm_Keys_31<Tkey,Tkeys_provide,Tkeys_require>(mpi_comm)
{
	this->traverse_keys_all =
		[](std::function<void(const Tkey&)> &func)
		{ throw std::logic_error("Function traverse not set."); };
}


template<typename Tkey, typename Tkeys_provide, typename Tkeys_require>
std::vector<std::vector<Tkey>> Comm_Keys_31<Tkey,Tkeys_provide,Tkeys_require>::trans(
	const Tkeys_provide &keys_provide_mine,
	const Tkeys_require &keys_require_mine)
{
	std::vector<std::vector<Tkey>> keys_trans_list(this->rank_size);
	int unfinish_rank = this->rank_size;
	std::vector<std::thread> threads;
	threads.reserve(this->rank_size+2);
	MPI_Status status_working;
	status_working.MPI_SOURCE=-1;

	threads.emplace_back(
		&Comm_Keys_31::send_keys_require_mine, this,
			std::cref(keys_require_mine) );

	threads.emplace_back(
		&Comm_Keys_31::intersection, this,
			std::cref(keys_provide_mine),
			std::cref(keys_require_mine),
			std::ref(keys_trans_list[rank_mine]));
	--unfinish_rank;

	while(unfinish_rank)
	{
		int flag_iprobe = false;
		MPI_Status status;
		MPI_CHECK( MPI_Iprobe( MPI_ANY_SOURCE, MPI_ANY_TAG, this->mpi_comm, &flag_iprobe, &status ) );
		if(flag_iprobe
			&& (status_working.MPI_SOURCE!=status.MPI_SOURCE))
		{
			threads.emplace_back(
				&Comm_Keys_31::recv_require_intersection, this,
					std::cref(keys_provide_mine), std::ref(keys_trans_list) );
			--unfinish_rank;
			status_working = status;
		}
		else
		{
			std::this_thread::yield();
		}
	}

	for(std::thread &t : threads)
		t.join();

	return keys_trans_list;
}


template<typename Tkey, typename Tkeys_provide, typename Tkeys_require>
void Comm_Keys_31<Tkey,Tkeys_provide,Tkeys_require>::send_keys_require_mine(
	const Tkeys_require &keys_require_mine)
{
	std::stringstream ss_isend;
	{
		cereal::BinaryOutputArchive ar(ss_isend);
		ar(keys_require_mine);
	}
	const std::size_t exponent_align = this->cereal_func.align_stringstream(ss_isend);
	const std::string str_isend = ss_isend.str();

	std::vector<MPI_Request> requests_isend(this->rank_size);
	for(int rank_recv_tmp=1; rank_recv_tmp<this->rank_size; ++rank_recv_tmp)
	{
		const int rank_recv = (this->rank_mine + rank_recv_tmp) % this->rank_size;
		this->cereal_func.mpi_isend(str_isend, exponent_align, rank_recv, this->tag_keys, this->mpi_comm, requests_isend[rank_recv]);
		std::this_thread::yield();
	}

	for(int rank_recv_tmp=1; rank_recv_tmp<this->rank_size; ++rank_recv_tmp)
	{
		const int rank_recv = (this->rank_mine + rank_recv_tmp) % this->rank_size;
		while(true)
		{
			int flag_finish = false;
			MPI_CHECK( MPI_Test(&requests_isend[rank_recv], &flag_finish, MPI_STATUS_IGNORE) );
			if(flag_finish)	break;
			std::this_thread::yield();
		}
	}
}


template<typename Tkey, typename Tkeys_provide, typename Tkeys_require>
void Comm_Keys_31<Tkey,Tkeys_provide,Tkeys_require>::recv_require_intersection(
	const Tkeys_provide &keys_provide_mine,
	std::vector<std::vector<Tkey>> &keys_trans_list)
{
	Tkeys_require keys_require;
	const MPI_Status status_recv = this->cereal_func.mpi_recv( this->mpi_comm,
		keys_require);
	const int rank_require = status_recv.MPI_SOURCE;
	assert(this->tag_keys==status_recv.MPI_TAG);

	this->intersection(keys_provide_mine, keys_require, keys_trans_list[rank_require]);
}


template<typename Tkey, typename Tkeys_provide, typename Tkeys_require>
void Comm_Keys_31_SenderTraversal<Tkey,Tkeys_provide,Tkeys_require>::intersection(
	const Tkeys_provide &keys_provide_mine,
	const Tkeys_require &keys_require,
	std::vector<Tkey> &keys_trans)
{
	std::function<void(const Tkey&)> inter = [&](const Tkey &key)
	{
		if(keys_require.judge(key))
			keys_trans.push_back(key);
	};
	this->traverse_keys_provide( keys_provide_mine, inter );
}


template<typename Tkey, typename Tkeys_provide, typename Tkeys_require>
void Comm_Keys_31_SenderJudge<Tkey,Tkeys_provide,Tkeys_require>::intersection(
	const Tkeys_provide &keys_provide_mine,
	const Tkeys_require &keys_require,
	std::vector<Tkey> &keys_trans)
{
	std::function<void(const Tkey&)> inter = [&](const Tkey &key)
	{
		if(keys_provide_mine.judge(key) && keys_require.judge(key))
			keys_trans.push_back(key);
	};
	this->traverse_keys_all(inter);
}

}

#undef MPI_CHECK