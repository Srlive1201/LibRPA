// ===================
//  Author: Peize Lin
//  date: 2022.06.22
// ===================

#pragma once

#include "Comm_Keys_31-gather.h"
#include "../global/Cereal_Func.h"

#include <mpi.h>
#include <omp.h>
#include <sstream>
#include <cereal/archives/binary.hpp>
#include <string>

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

#if MPI_VERSION>=4
	std::vector<MPI_Count> sss_size;
	std::vector<MPI_Aint> sss_displs;
#else
	std::vector<int> sss_size;
	std::vector<int> sss_displs;
#endif
	std::vector<char> buffer_recv;
	this->send_keys_require_mine( keys_require_mine,
		sss_size, sss_displs, buffer_recv);

	#pragma omp parallel for
	for(int rank_require=0; rank_require<this->rank_size; ++rank_require)
	{
		Tkeys_require keys_require;
		std::stringstream ss_recv;
		ss_recv.rdbuf()->pubsetbuf(buffer_recv.data()+sss_displs[rank_require], sss_size[rank_require]);
		{
			cereal::BinaryInputArchive ar(ss_recv);
			ar(keys_require);
		}

		this->intersection(keys_provide_mine, keys_require, keys_trans_list[rank_require]);
	}

	return keys_trans_list;
}


#if MPI_VERSION>=4
template<typename Tkey, typename Tkeys_provide, typename Tkeys_require>
void Comm_Keys_31<Tkey,Tkeys_provide,Tkeys_require>::send_keys_require_mine(
	const Tkeys_require &keys_require_mine,
	std::vector<MPI_Count> &sss_size,
	std::vector<MPI_Aint> &sss_displs,
	std::vector<char> &buffer_recv)
#else
template<typename Tkey, typename Tkeys_provide, typename Tkeys_require>
void Comm_Keys_31<Tkey,Tkeys_provide,Tkeys_require>::send_keys_require_mine(
	const Tkeys_require &keys_require_mine,
	std::vector<int> &sss_size,
	std::vector<int> &sss_displs,
	std::vector<char> &buffer_recv)
#endif
{
	std::stringstream ss_send;
	{
		cereal::BinaryOutputArchive ar(ss_send);
		ar(keys_require_mine);
	}

	sss_size.resize(this->rank_size);
#if MPI_VERSION>=4
	const MPI_Count ss_size = ss_send.str().size();
	MPI_CHECK( MPI_Allgather( &ss_size, 1, MPI_COUNT, sss_size.data(), 1, MPI_COUNT, this->mpi_comm ) );
#else
	const int ss_size = ss_send.str().size();
	MPI_CHECK( MPI_Allgather( &ss_size, 1, MPI_INT, sss_size.data(), 1, MPI_INT, this->mpi_comm ) );
#endif

	sss_displs.resize(this->rank_size);
	sss_displs[0] = 0;
	for(int i=1; i<this->rank_size; ++i)
		sss_displs[i] = sss_displs[i-1] + sss_size[i-1];

	buffer_recv.resize(sss_displs.back() + sss_size.back());
#if MPI_VERSION>=4
	MPI_CHECK( MPI_Allgatherv_c( ss_send.str().c_str(), ss_send.str().size(), MPI_CHAR, buffer_recv.data(), sss_size.data(), sss_displs.data(), MPI_CHAR, this->mpi_comm ) );
#else
	MPI_CHECK( MPI_Allgatherv  ( ss_send.str().c_str(), ss_send.str().size(), MPI_CHAR, buffer_recv.data(), sss_size.data(), sss_displs.data(), MPI_CHAR, this->mpi_comm ) );
#endif
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