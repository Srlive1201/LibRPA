//=======================
// AUTHOR : Peize Lin
// DATE :   2022-01-05
//=======================

#pragma once

#include "Comm_Trans.h"
#include "../global/Cereal_Func.h"

#include <vector>
#include <string>
#include <stdexcept>

#include <cereal/archives/binary.hpp>
#include <cereal/types/tuple.hpp>
#include <cereal/types/map.hpp>

#define MPI_CHECK(x) if((x)!=MPI_SUCCESS)	throw std::runtime_error(std::string(__FILE__)+" line "+std::to_string(__LINE__));

namespace Comm
{

template<typename Tkey, typename Tvalue, typename Tdatas_isend, typename Tdatas_recv>
Comm_Trans<Tkey,Tvalue,Tdatas_isend,Tdatas_recv>::Comm_Trans(const MPI_Comm &mpi_comm_in)
	:mpi_comm(mpi_comm_in)
{
	MPI_CHECK (MPI_Comm_size (this->mpi_comm, &this->comm_size));
	MPI_CHECK (MPI_Comm_rank (this->mpi_comm, &this->rank_mine));

	this->set_value_recv
		= [](Tkey &&key, Tvalue &&value, Tdatas_recv &datas_recv)
		{ throw std::logic_error("Function set_value not set."); };
	this->traverse_isend
		= [](const Tdatas_isend &datas_isend, const int rank_isend, std::function<void(const Tkey&, const Tvalue&)> &func)
		{ throw std::logic_error("Function traverse not set."); };
	this->init_datas_local
		= [](const int rank_recv) -> Tdatas_recv
		{ throw std::logic_error("Function init_datas_local not set."); };
	this->add_datas
		= [](Tdatas_recv &&datas_local, Tdatas_recv &datas_recv)
		{ throw std::logic_error("Function add_datas not set."); };
}

/*
template<typename Tkey, typename Tvalue, typename Tdatas_isend, typename Tdatas_recv>
Comm_Trans<Tkey,Tvalue,Tdatas_isend,Tdatas_recv>::Comm_Trans(const Comm_Trans &com)
	:mpi_comm(com.mpi_comm)
{
	//ofs<<"C"<<" ";
	MPI_CHECK (MPI_Comm_size (this->mpi_comm, &this->comm_size));
	MPI_CHECK (MPI_Comm_rank (this->mpi_comm, &this->rank_mine));
	this->set_value_recv = com.set_value_recv;
	this->traverse_isend = com.traverse_isend;
	this->flag_lock_set_value = com.flag_lock_set_value;
	this->init_datas_local = com.init_datas_local;
	this->add_datas = com.add_datas;
}
*/

template<typename Tkey, typename Tvalue, typename Tdatas_isend, typename Tdatas_recv>
void Comm_Trans<Tkey,Tvalue,Tdatas_isend,Tdatas_recv>::communicate(
	const Tdatas_isend &datas_isend,
	Tdatas_recv &datas_recv)
{
	// initialization
	int rank_isend_tmp = 0;
	int rank_recv_working = -1;

	std::vector<MPI_Request> requests_isend(comm_size);
	std::vector<std::string> strs_isend(comm_size);
	std::vector<std::future<void>> futures_isend(comm_size);
	std::vector<std::future<void>> futures_recv(comm_size);
	std::atomic_flag lock_set_value = ATOMIC_FLAG_INIT;
	std::atomic<std::size_t> memory_max_isend(0);
	std::atomic<std::size_t> memory_max_recv(0);

	std::future<void> future_post_process = std::async (std::launch::async,
		&Comm_Trans::post_process, this,
		std::ref(requests_isend), std::ref(strs_isend), std::ref(futures_isend), std::ref(futures_recv));

	while (future_post_process.wait_for(std::chrono::seconds(0)) != std::future_status::ready)
	{
		int flag_iprobe=0;
		MPI_Status status_recv;
		MPI_Message message_recv;
		MPI_CHECK (MPI_Improbe(MPI_ANY_SOURCE, this->tag_data, this->mpi_comm, &flag_iprobe, &message_recv, &status_recv));
		if (flag_iprobe && rank_recv_working!=status_recv.MPI_SOURCE && memory_enough(memory_max_recv))
		{
			futures_recv[status_recv.MPI_SOURCE] = std::async (std::launch::async,
				&Comm_Trans::recv_data, this,
				std::ref(datas_recv), status_recv, message_recv, std::ref(lock_set_value), std::ref(memory_max_recv));
			rank_recv_working = status_recv.MPI_SOURCE;
		}

		if (rank_isend_tmp < this->comm_size && memory_enough(memory_max_isend))
		{
			const int rank_isend = (rank_isend_tmp + this->rank_mine) % this->comm_size;
			futures_isend[rank_isend] = std::async (std::launch::async,
				&Comm_Trans::isend_data, this,
				rank_isend, std::cref(datas_isend), std::ref(strs_isend[rank_isend]), std::ref(requests_isend[rank_isend]), std::ref(memory_max_isend));
			++rank_isend_tmp;
		}
	}
	future_post_process.get();
}


template<typename Tkey, typename Tvalue, typename Tdatas_isend, typename Tdatas_recv>
void Comm_Trans<Tkey,Tvalue,Tdatas_isend,Tdatas_recv>::isend_data(
	const int rank_isend,
	const Tdatas_isend &datas_isend,
	std::string &str_isend,
	MPI_Request &request_isend,
	std::atomic<std::size_t> &memory_max_isend)
{
	std::stringstream ss_isend;
	{
		cereal::BinaryOutputArchive oar(ss_isend);

		size_t size_item = 0;
		oar(size_item);					// 占位

		std::function<void(const Tkey&, const Tvalue&)> archive_data = [&oar, &size_item](
			const Tkey &key, const Tvalue &value)
		{
			oar(key, value);
			++size_item;
		};
		this->traverse_isend(datas_isend, rank_isend, archive_data);

		ss_isend.rdbuf()->pubseekpos(0);		// 返回size_item的占位，序列化真正的size_item值
		oar(size_item);
	} // end cereal::BinaryOutputArchive
	const std::size_t exponent_align = this->cereal_func.align_stringstream(ss_isend);
	str_isend = ss_isend.str();
	memory_max_isend.store( std::max(str_isend.size()*sizeof(char), memory_max_isend.load()) );
	this->cereal_func.mpi_isend(str_isend, exponent_align, rank_isend, this->tag_data, this->mpi_comm, request_isend);
}



template<typename Tkey, typename Tvalue, typename Tdatas_isend, typename Tdatas_recv>
void Comm_Trans<Tkey,Tvalue,Tdatas_isend,Tdatas_recv>::recv_data (
	Tdatas_recv &datas_recv,
	const MPI_Status status_recv,
	MPI_Message message_recv,
	std::atomic_flag &lock_set_value,
	std::atomic<std::size_t> &memory_max_recv)
{
	std::vector<char> buffer_recv = this->cereal_func.mpi_mrecv(message_recv, status_recv);

	std::stringstream ss_recv;
	ss_recv.rdbuf()->pubsetbuf(buffer_recv.data(), buffer_recv.size());
	memory_max_recv.store( std::max(buffer_recv.size()*sizeof(char), memory_max_recv.load()) );

	{
		cereal::BinaryInputArchive iar(ss_recv);
		size_t size_item;	iar(size_item);

		if (this->flag_lock_set_value==Comm_Tools::Lock_Type::Lock_free)
		{
			for (size_t i=0; i<size_item; ++i)
			{
				Tkey key;
				Tvalue value;
				iar(key, value);

				this->set_value_recv(std::move(key), std::move(value), datas_recv);
			}
		}
		else if (this->flag_lock_set_value==Comm_Tools::Lock_Type::Lock_item)
		{
			for (size_t i=0; i<size_item; ++i)
			{
				Tkey key;
				Tvalue value;
				iar(key, value);

				while (lock_set_value.test_and_set(std::memory_order_seq_cst)) std::this_thread::yield();
				this->set_value_recv(std::move(key), std::move(value), datas_recv);
				lock_set_value.clear(std::memory_order_seq_cst);
			}
		}
		else if (this->flag_lock_set_value==Comm_Tools::Lock_Type::Lock_Process)
		{
			while (lock_set_value.test_and_set(std::memory_order_seq_cst)) std::this_thread::yield();
			for (size_t i=0; i<size_item; ++i)
			{
				Tkey key;
				Tvalue value;
				iar(key, value);

				this->set_value_recv(std::move(key), std::move(value), datas_recv);
			}
			lock_set_value.clear(std::memory_order_seq_cst);
		}
		else if (this->flag_lock_set_value==Comm_Tools::Lock_Type::Copy_merge)
		{
			Tdatas_recv datas_local = this->init_datas_local (status_recv.MPI_SOURCE);
			for (size_t i=0; i<size_item; ++i)
			{
				Tkey key;
				Tvalue value;
				iar(key, value);

				this->set_value_recv (std::move(key), std::move(value), datas_local);
			}
			while (lock_set_value.test_and_set(std::memory_order_seq_cst)) std::this_thread::yield();
			this->add_datas (std::move(datas_local), datas_recv);
			lock_set_value.clear(std::memory_order_seq_cst);
		}
		else
		{
			throw std::invalid_argument(
				+" file "+std::string(__FILE__)
				+" line "+std::to_string(__LINE__)
				+" rank_mine "+std::to_string(this->rank_mine)
				+" rank_recv "+std::to_string(status_recv.MPI_SOURCE));
		}
	} // end cereal::BinaryInputArchive
}


template<typename Tkey, typename Tvalue, typename Tdatas_isend, typename Tdatas_recv>
void Comm_Trans<Tkey,Tvalue,Tdatas_isend,Tdatas_recv>::post_process(
	std::vector<MPI_Request> &requests_isend,
	std::vector<std::string> &strs_isend,
	std::vector<std::future<void>> &futures_isend,
	std::vector<std::future<void>> &futures_recv) const
{
	int rank_isend_free_tmp = 0;
	int rank_recv_free_tmp = 0;
	while (rank_isend_free_tmp < this->comm_size
		|| rank_recv_free_tmp < this->comm_size)
	{
		while (rank_isend_free_tmp < this->comm_size)
		{
			const int rank_isend_free = (this->rank_mine+rank_isend_free_tmp)%this->comm_size;
			if (futures_isend[rank_isend_free].valid()
				&& futures_isend[rank_isend_free].wait_for(std::chrono::seconds(0)) == std::future_status::ready)
			{
				int flag_finish=0;
				MPI_CHECK (MPI_Test (&(requests_isend[rank_isend_free]), &flag_finish, MPI_STATUS_IGNORE));
				if (flag_finish)
				{
	//				MPI_CHECK (MPI_Request_free (&requests_isend[rank_isend_free]));
					futures_isend[rank_isend_free].get();
					strs_isend[rank_isend_free].clear();
					++rank_isend_free_tmp;
				}
				else{ break; }
			}
			else{ break; }
		}

		while (rank_recv_free_tmp < this->comm_size)
		{
			const int rank_recv_free = (this->rank_mine+rank_recv_free_tmp)%this->comm_size;
			if (futures_recv[rank_recv_free].valid()
				&& futures_recv[rank_recv_free].wait_for(std::chrono::seconds(0)) == std::future_status::ready)
			{
				futures_recv[rank_recv_free].get();
				++rank_recv_free_tmp;
			}
			else{ break; }
		}

		std::this_thread::yield();
	}
}

}

#undef MPI_CHECK

/*
get_send_keys()
{

	if(unique)
	{
		for(irank in all)
			send(irank_send, atom_pairs_remove);
	}
}
*/