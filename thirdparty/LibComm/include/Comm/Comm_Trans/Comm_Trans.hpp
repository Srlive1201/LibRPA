//=======================
// AUTHOR : Peize Lin
// DATE :   2022-01-05
//=======================

#pragma once

#include "Comm_Trans.h"
#include "Memory_Check.h"

#include <vector>
#include <queue>
#include <string>
#include <stdexcept>
#include <cassert>
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <istream>
#include <limits>
#include <streambuf>

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
	const char *mode_env = std::getenv("LIBCOMM_TRANS_MODE");
	const std::string mode = mode_env == nullptr || mode_env[0] == '\0'
		? "sendrecv_ring" : mode_env;
	const int mode_code = mode == "nonblocking" ? 0 : mode == "sendrecv_ring" ? 1 : -1;
	int mode_min = 0;
	int mode_max = 0;
	MPI_CHECK(MPI_Allreduce(&mode_code, &mode_min, 1, MPI_INT, MPI_MIN, this->mpi_comm));
	MPI_CHECK(MPI_Allreduce(&mode_code, &mode_max, 1, MPI_INT, MPI_MAX, this->mpi_comm));
	if(mode_min != mode_max)
	{
		throw std::invalid_argument("Inconsistent LIBCOMM_TRANS_MODE across MPI ranks");
	}
	if(mode_code < 0)
	{
		throw std::invalid_argument("Unknown LIBCOMM_TRANS_MODE: " + mode);
	}

	if(mode_code == 0)
	{
		this->communicate_nonblocking(datas_isend, datas_recv);
	}
	else
	{
		this->communicate_sendrecv_ring(datas_isend, datas_recv);
	}
}


template<typename Tkey, typename Tvalue, typename Tdatas_isend, typename Tdatas_recv>
void Comm_Trans<Tkey,Tvalue,Tdatas_isend,Tdatas_recv>::communicate_nonblocking(
	const Tdatas_isend &datas_isend,
	Tdatas_recv &datas_recv)
{
	std::vector<MPI_Request> requests_isend(this->comm_size);
	std::vector<std::string> buffers_isend(this->comm_size);
	std::vector<std::future<std::size_t>> futures_oar(this->comm_size);
	std::vector<std::vector<char>> buffers_recv(this->comm_size);
	std::queue<std::pair<MPI_Status,MPI_Message>> status_message_s_recv;
	std::vector<std::atomic<State_Send>> states_send(this->comm_size);
	for(std::atomic<State_Send> & state_send : states_send)
		state_send = State_Send::unstart;
	std::vector<std::atomic<State_Recv>> states_recv(this->comm_size);
	for(std::atomic<State_Recv> & state_recv : states_recv)
		state_recv = State_Recv::unstart;
	std::atomic_flag lock_set_value = ATOMIC_FLAG_INIT;
	Memory_Check memory(states_send, states_recv);

	// initialization
	int rank_isend_tmp1=1, rank_isend_tmp2=1, rank_isend_tmp3=1;

	while (!check_finish(states_send, states_recv))
	{
		while(true)
		{
			int flag_iprobe=0;
			MPI_Status status_recv;
			MPI_Message message_recv;
			MPI_CHECK (MPI_Improbe(MPI_ANY_SOURCE, this->tag_data, this->mpi_comm, &flag_iprobe, &message_recv, &status_recv));
			if (flag_iprobe)
				status_message_s_recv.emplace(status_recv, message_recv);
			else
				break;
		}

		if (!status_message_s_recv.empty() && memory.enough_recv())
		{
			const MPI_Status status_recv = status_message_s_recv.front().first;
			const MPI_Message message_recv = status_message_s_recv.front().second;
			const int rank_recv = status_recv.MPI_SOURCE;
			status_message_s_recv.pop();

			this->recv_data(
				status_recv,
				message_recv,
				memory,
				buffers_recv[rank_recv],
				states_recv[rank_recv]);

			std::async (std::launch::async,
				&Comm_Trans::iar_data, this,
					rank_recv,
					std::ref(buffers_recv[rank_recv]),
					std::ref(lock_set_value),
					std::ref(datas_recv),
					std::ref(states_recv[rank_recv])).wait();
		}

		if (rank_isend_tmp1<this->comm_size+1 && memory.enough_send())
		{
			const int rank_isend = (rank_isend_tmp1 + this->rank_mine) % this->comm_size;
			futures_oar[rank_isend] = std::async (std::launch::async,
				&Comm_Trans::oar_data, this,
					rank_isend,
					std::cref(datas_isend),
					std::ref(buffers_isend[rank_isend]),
					std::ref(states_send[rank_isend]),
					std::ref(memory));
			++rank_isend_tmp1;
		}

		while(rank_isend_tmp2<rank_isend_tmp1)
		{
			const int rank_isend = (rank_isend_tmp2 + this->rank_mine) % this->comm_size;
			if(futures_oar[rank_isend].valid())
			{
				const std::size_t exponent_align = futures_oar[rank_isend].get();
				this->isend_data(
					rank_isend,
					exponent_align,
					buffers_isend[rank_isend],
					requests_isend[rank_isend],
					states_send[rank_isend]);
				++rank_isend_tmp2;
			}
		}

		while(rank_isend_tmp3<rank_isend_tmp2)
		{
			const int rank_isend = (rank_isend_tmp3 + this->rank_mine) % this->comm_size;
			if(states_send[rank_isend] == State_Send::begin_isend)
			{
				int flag_finish=0;
				MPI_CHECK (MPI_Test (&(requests_isend[rank_isend]), &flag_finish, MPI_STATUS_IGNORE));
				if (flag_finish)
				{
					//MPI_CHECK (MPI_Request_free (&requests_isend[rank_isend]));
					buffers_isend[rank_isend].clear();
					buffers_isend[rank_isend].shrink_to_fit();
					states_send[rank_isend] = State_Send::finish_isend;
					++rank_isend_tmp3;
				}
				else {break;}
			}
			else {break;}
		}
	}
}


template<typename Tkey, typename Tvalue, typename Tdatas_isend, typename Tdatas_recv>
std::vector<char> Comm_Trans<Tkey,Tvalue,Tdatas_isend,Tdatas_recv>::serialize_for_rank(
	const int rank_isend,
	const Tdatas_isend &datas_isend) const
{
	class VectorOutputBuffer : public std::streambuf
	{
	public:
		std::vector<char> release()
		{
			return std::move(this->buffer);
		}

	protected:
		std::streamsize xsputn(const char *data, const std::streamsize size) override
		{
			if(size <= 0) return 0;
			const std::size_t count = static_cast<std::size_t>(size);
			if(count > std::numeric_limits<std::size_t>::max() - this->position)
				throw std::length_error("LibComm ring send payload exceeds size_t");
			const std::size_t end = this->position + count;
			if(end > this->buffer.size()) this->buffer.resize(end);
			std::copy(data, data + count, this->buffer.begin() + this->position);
			this->position = end;
			return size;
		}

		int_type overflow(const int_type value) override
		{
			if(traits_type::eq_int_type(value, traits_type::eof()))
				return traits_type::not_eof(value);
			const char data = traits_type::to_char_type(value);
			this->xsputn(&data, 1);
			return value;
		}

		pos_type seekoff(const off_type offset, const std::ios_base::seekdir direction,
		                 const std::ios_base::openmode mode) override
		{
			if((mode & std::ios_base::out) == 0) return pos_type(off_type(-1));
			std::size_t base = 0;
			if(direction == std::ios_base::cur) base = this->position;
			else if(direction == std::ios_base::end) base = this->buffer.size();
			else if(direction != std::ios_base::beg) return pos_type(off_type(-1));
			if(offset < 0)
			{
				const unsigned long long magnitude =
					static_cast<unsigned long long>(-(offset + 1)) + 1;
				if(magnitude > base) return pos_type(off_type(-1));
				this->position = base - static_cast<std::size_t>(magnitude);
			}
			else
			{
				if(static_cast<unsigned long long>(offset)
					> std::numeric_limits<std::size_t>::max() - base)
					return pos_type(off_type(-1));
				this->position = base + static_cast<std::size_t>(offset);
			}
			return pos_type(this->position);
		}

		pos_type seekpos(const pos_type position, const std::ios_base::openmode mode) override
		{
			return this->seekoff(static_cast<off_type>(position), std::ios_base::beg, mode);
		}

	private:
		std::vector<char> buffer;
		std::size_t position = 0;
	};

	VectorOutputBuffer stream_buffer;
	std::ostream stream(&stream_buffer);
	{
		cereal::BinaryOutputArchive archive(stream);
		std::size_t size_item = 0;
		archive(size_item);
		std::function<void(const Tkey&, const Tvalue&)> archive_data = [&archive, &size_item](
			const Tkey &key, const Tvalue &value)
		{
			archive(key, value);
			++size_item;
		};
		this->traverse_isend(datas_isend, rank_isend, archive_data);
		stream.seekp(0);
		archive(size_item);
	}
	if(!stream)
		throw std::runtime_error("Failed to serialize LibComm ring payload");
	return stream_buffer.release();
}


template<typename Tkey, typename Tvalue, typename Tdatas_isend, typename Tdatas_recv>
void Comm_Trans<Tkey,Tvalue,Tdatas_isend,Tdatas_recv>::deserialize_from_rank(
	const int rank_recv,
	const std::vector<char> &buffer_recv,
	Tdatas_recv &datas_recv) const
{
	class ReadOnlyVectorBuffer : public std::streambuf
	{
	public:
		explicit ReadOnlyVectorBuffer(const std::vector<char> &buffer)
		{
			if(buffer.empty())
			{
				this->setg(nullptr, nullptr, nullptr);
				return;
			}
			char *begin = const_cast<char *>(buffer.data());
			this->setg(begin, begin, begin + buffer.size());
		}
	};
	ReadOnlyVectorBuffer stream_buffer(buffer_recv);
	std::istream stream(&stream_buffer);
	cereal::BinaryInputArchive archive(stream);
	std::size_t size_item = 0;
	archive(size_item);

	if(this->flag_lock_set_value == Comm_Tools::Lock_Type::Copy_merge)
	{
		Tdatas_recv datas_local = this->init_datas_local(rank_recv);
		for(std::size_t i = 0; i != size_item; ++i)
		{
			Tkey key;
			Tvalue value;
			archive(key, value);
			this->set_value_recv(std::move(key), std::move(value), datas_local);
		}
		this->add_datas(std::move(datas_local), datas_recv);
		return;
	}

	if(this->flag_lock_set_value != Comm_Tools::Lock_Type::Lock_free &&
	   this->flag_lock_set_value != Comm_Tools::Lock_Type::Lock_item &&
	   this->flag_lock_set_value != Comm_Tools::Lock_Type::Lock_Process)
	{
		throw std::invalid_argument(
			"Unknown flag_lock_set_value on rank " + std::to_string(this->rank_mine));
	}

	// The ring backend deserializes one peer at a time on the calling thread.
	for(std::size_t i = 0; i != size_item; ++i)
	{
		Tkey key;
		Tvalue value;
		archive(key, value);
		this->set_value_recv(std::move(key), std::move(value), datas_recv);
	}
}


template<typename Tkey, typename Tvalue, typename Tdatas_isend, typename Tdatas_recv>
void Comm_Trans<Tkey,Tvalue,Tdatas_isend,Tdatas_recv>::communicate_sendrecv_ring(
	const Tdatas_isend &datas_isend,
	Tdatas_recv &datas_recv)
{
	// Keep only one peer's serialized send/receive buffers live at a time.  The
	// duplicated communicator isolates the ring's size and payload messages.
	unsigned long long chunk_value = 64ULL << 20;
	int chunk_valid = 1;
	if(const char *chunk_env = std::getenv("LIBCOMM_TRANS_CHUNK_BYTES"))
	{
		try
		{
			std::size_t parsed = 0;
			chunk_value = std::stoull(chunk_env, &parsed);
			if(parsed != std::string(chunk_env).size() || chunk_value == 0 || chunk_env[0] == '-')
				chunk_valid = 0;
		}
		catch(...)
		{
			chunk_valid = 0;
		}
	}
	int all_chunks_valid = 0;
	MPI_CHECK(MPI_Allreduce(
		&chunk_valid, &all_chunks_valid, 1, MPI_INT, MPI_MIN, this->mpi_comm));
	if(!all_chunks_valid)
		throw std::invalid_argument("Invalid LIBCOMM_TRANS_CHUNK_BYTES on an MPI rank");
	chunk_value = std::min(
		chunk_value, static_cast<unsigned long long>(std::numeric_limits<int>::max()));
	unsigned long long chunk_min = 0;
	unsigned long long chunk_max = 0;
	MPI_CHECK(MPI_Allreduce(
		&chunk_value, &chunk_min, 1, MPI_UNSIGNED_LONG_LONG, MPI_MIN, this->mpi_comm));
	MPI_CHECK(MPI_Allreduce(
		&chunk_value, &chunk_max, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, this->mpi_comm));
	if(chunk_min != chunk_max)
		throw std::invalid_argument("Inconsistent LIBCOMM_TRANS_CHUNK_BYTES across MPI ranks");
	const std::size_t chunk_bytes = static_cast<std::size_t>(chunk_value);
	const bool trace = []()
	{
		const char *trace_env = std::getenv("LIBCOMM_TRANS_TRACE");
		return trace_env != nullptr && std::string(trace_env) != "0";
	}();

	struct RingCommGuard
	{
		MPI_Comm value = MPI_COMM_NULL;
		~RingCommGuard()
		{
			if(this->value != MPI_COMM_NULL)
				MPI_Comm_free(&this->value);
		}
	};
	RingCommGuard ring_comm;
	MPI_CHECK(MPI_Comm_dup(this->mpi_comm, &ring_comm.value));
	auto throw_if_rank_failed = [&](const int local_failed, const char *phase)
	{
		const int local_rank = local_failed ? this->rank_mine : this->comm_size;
		int first_failed_rank = this->comm_size;
		MPI_CHECK(MPI_Allreduce(
			&local_rank, &first_failed_rank, 1, MPI_INT, MPI_MIN, ring_comm.value));
		if(first_failed_rank != this->comm_size)
		{
			throw std::runtime_error(
				std::string("LibComm ring ") + phase + " failed on MPI rank "
				+ std::to_string(first_failed_rank));
		}
	};
	constexpr int tag_size = 0;
	constexpr int tag_payload = 1;

	for(int step = 0; step != this->comm_size; ++step)
	{
		const int rank_isend = (this->rank_mine + step) % this->comm_size;
		const int rank_recv = (this->rank_mine - step + this->comm_size) % this->comm_size;
		std::vector<char> buffer_isend;
		int serialization_failed = 0;
		try
		{
			buffer_isend = this->serialize_for_rank(rank_isend, datas_isend);
		}
		catch(...)
		{
			serialization_failed = 1;
		}
		throw_if_rank_failed(serialization_failed, "serialization");
		const unsigned long long size_isend = buffer_isend.size();
		unsigned long long size_recv = 0;
		MPI_CHECK(MPI_Sendrecv(
			&size_isend, 1, MPI_UNSIGNED_LONG_LONG, rank_isend, tag_size,
			&size_recv, 1, MPI_UNSIGNED_LONG_LONG, rank_recv, tag_size,
			ring_comm.value, MPI_STATUS_IGNORE));

		std::vector<char> buffer_recv;
		int allocation_failed = 0;
		try
		{
			if(size_recv > static_cast<unsigned long long>(
				std::numeric_limits<std::size_t>::max()))
			{
				throw std::length_error("LibComm ring receive payload exceeds size_t");
			}
			buffer_recv.resize(static_cast<std::size_t>(size_recv));
		}
		catch(...)
		{
			allocation_failed = 1;
		}
		throw_if_rank_failed(allocation_failed, "receive allocation");
		const std::size_t send_rounds =
			(buffer_isend.size() + chunk_bytes - 1) / chunk_bytes;
		const std::size_t recv_rounds =
			(buffer_recv.size() + chunk_bytes - 1) / chunk_bytes;
		unsigned long long local_rounds = std::max(send_rounds, recv_rounds);
		unsigned long long global_rounds = 0;
		MPI_CHECK(MPI_Allreduce(
			&local_rounds, &global_rounds, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, ring_comm.value));

		if(trace)
		{
			std::cerr << "LIBCOMM_TRANS_RING rank=" << this->rank_mine
			          << " step=" << step
			          << " send_rank=" << rank_isend
			          << " recv_rank=" << rank_recv
			          << " send_bytes=" << size_isend
			          << " recv_bytes=" << size_recv
			          << " rounds=" << global_rounds << std::endl;
		}

		for(std::size_t round = 0; round != static_cast<std::size_t>(global_rounds); ++round)
		{
			const std::size_t send_offset = std::min(round * chunk_bytes, buffer_isend.size());
			const std::size_t recv_offset = std::min(round * chunk_bytes, buffer_recv.size());
			const int send_count = static_cast<int>(std::min(
				chunk_bytes, buffer_isend.size() - send_offset));
			const int recv_count = static_cast<int>(std::min(
				chunk_bytes, buffer_recv.size() - recv_offset));
			const char *send_ptr = send_count == 0 ? nullptr : buffer_isend.data() + send_offset;
			char *recv_ptr = recv_count == 0 ? nullptr : buffer_recv.data() + recv_offset;
			MPI_CHECK(MPI_Sendrecv(
				send_ptr, send_count, MPI_BYTE, rank_isend, tag_payload,
				recv_ptr, recv_count, MPI_BYTE, rank_recv, tag_payload,
				ring_comm.value, MPI_STATUS_IGNORE));
		}

		int deserialization_failed = 0;
		try
		{
			this->deserialize_from_rank(rank_recv, buffer_recv, datas_recv);
		}
		catch(...)
		{
			deserialization_failed = 1;
		}
		throw_if_rank_failed(deserialization_failed, "deserialization");
	}

	MPI_CHECK(MPI_Comm_free(&ring_comm.value));
}


template<typename Tkey, typename Tvalue, typename Tdatas_isend, typename Tdatas_recv>
std::size_t Comm_Trans<Tkey,Tvalue,Tdatas_isend,Tdatas_recv>::oar_data(
	const int rank_isend,
	const Tdatas_isend &datas_isend,
	std::string &buffer_isend,
	std::atomic<State_Send> &state_send,
	Memory_Check &memory)
{
	assert(state_send == State_Send::unstart);
	state_send = State_Send::begin_oar;
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
	buffer_isend = std::move(ss_isend.str());
	memory.set_max_used_send(buffer_isend.size()*sizeof(char));
	state_send = State_Send::finish_oar;
	return exponent_align;
}


template<typename Tkey, typename Tvalue, typename Tdatas_isend, typename Tdatas_recv>
void Comm_Trans<Tkey,Tvalue,Tdatas_isend,Tdatas_recv>::isend_data(
	const int rank_isend,
	const std::size_t exponent_align,
	std::string &buffer_isend,
	MPI_Request &request_isend,
	std::atomic<State_Send> &state_send)
{
	assert(state_send == State_Send::finish_oar);
	this->cereal_func.mpi_isend(buffer_isend, exponent_align, rank_isend, this->tag_data, this->mpi_comm, request_isend);
	state_send = State_Send::begin_isend;
}


template<typename Tkey, typename Tvalue, typename Tdatas_isend, typename Tdatas_recv>
void Comm_Trans<Tkey,Tvalue,Tdatas_isend,Tdatas_recv>::recv_data (
	const MPI_Status status_recv,
	const MPI_Message message_recv,
	Memory_Check &memory,
	std::vector<char> &buffer_recv,
	std::atomic<State_Recv> &state_recv)
{
	assert(state_recv == State_Recv::unstart);
	state_recv = State_Recv::begin_recv;
	MPI_Message message_recv_tmp = message_recv;
	buffer_recv = this->cereal_func.mpi_mrecv(message_recv_tmp, status_recv);
	assert(message_recv_tmp == MPI_MESSAGE_NULL);
	memory.set_max_used_recv(buffer_recv.size()*sizeof(char));
	state_recv = State_Recv::finish_recv;
}


template<typename Tkey, typename Tvalue, typename Tdatas_isend, typename Tdatas_recv>
void Comm_Trans<Tkey,Tvalue,Tdatas_isend,Tdatas_recv>::iar_data (
	const int rank_recv,
	std::vector<char> &buffer_recv,
	std::atomic_flag &lock_set_value,
	Tdatas_recv &datas_recv,
	std::atomic<State_Recv> &state_recv) const
{
	assert(state_recv == State_Recv::finish_recv);
	state_recv = State_Recv::begin_iar;
	std::stringstream ss_recv;
	ss_recv.rdbuf()->pubsetbuf(buffer_recv.data(), buffer_recv.size());
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
			Tdatas_recv datas_local = this->init_datas_local (rank_recv);
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
				+" rank_recv "+std::to_string(rank_recv));
		}
	} // end cereal::BinaryInputArchive
	buffer_recv.clear();
	buffer_recv.shrink_to_fit();
	state_recv = State_Recv::finish_iar;
}


template<typename Tkey, typename Tvalue, typename Tdatas_isend, typename Tdatas_recv>
bool Comm_Trans<Tkey,Tvalue,Tdatas_isend,Tdatas_recv>::check_finish(
	const std::vector<std::atomic<State_Send>> &states_send,
	const std::vector<std::atomic<State_Recv>> &states_recv) const
{
	for(int rank_isend_tmp=this->comm_size; rank_isend_tmp>0; --rank_isend_tmp)
	{
		const int rank_isend = (rank_isend_tmp + this->rank_mine) % this->comm_size;
		if(states_send[rank_isend] != State_Send::finish_isend)
			return false;
	}
	for(int rank_recv_tmp=0; rank_recv_tmp<this->comm_size; ++rank_recv_tmp)
	{
		const int rank_recv = (rank_recv_tmp + this->rank_mine) % this->comm_size;
		if(states_recv[rank_recv] != State_Recv::finish_iar)
			return false;
	}
	return true;
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
