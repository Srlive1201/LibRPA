// ===================
//  Author: Peize Lin
//  date: 2022.05.01
// ===================

#pragma once

#include "Cereal_Func.h"
#include "Cereal_Types.h"
#include "MPI_Wrapper.h"
#include "Global_Func.h"

#include <cereal/archives/binary.hpp>
#include <sstream>
#include <mpi.h>
#include <cmath>
#include <limits>

#define MPI_CHECK(x) if((x)!=MPI_SUCCESS)	throw std::runtime_error(std::string(__FILE__)+" line "+std::to_string(__LINE__));

namespace Comm
{

inline Cereal_Func::Cereal_Func()
{
  #if MPI_VERSION>=4
	using int_type = MPI_Count;
  #else
	using int_type = int;
  #endif

	// assuming MPI communication is less than max(1TB, memory availablle now) for most case,
	// so initialize here to avoid thread conflict in the future.
	const std::size_t TB = std::size_t(1)<<12;
	const std::size_t memory_max = std::max(TB, Global_Func::memory_available());
	const std::size_t times = std::ceil( double(memory_max) / double(std::numeric_limits<int_type>::max()) );
	const std::size_t exponent_align = std::ceil( std::log(times) / std::log(2) );
	this->char_contiguous.resize(exponent_align);
}	

// every 2^exponent_align char concatenate to 1 word
	// <<exponent_align means *2^exponent_align
	// >>exponent_align means /2^exponent_align
inline std::size_t Cereal_Func::align_stringstream(std::stringstream &ss)
{
  #if MPI_VERSION>=4
	using int_type = MPI_Count;
  #else
	using int_type = int;
  #endif

	const std::size_t size_old = ss.str().size();				// Inefficient, should be optimized
	const std::size_t times = std::ceil( double(size_old) / double(std::numeric_limits<int_type>::max()) );
	const std::size_t exponent_align = std::ceil( std::log(times) / std::log(2) );
	this->char_contiguous.resize(exponent_align);
	constexpr char c0 = 0;
	const std::size_t size_align = 1<<exponent_align;
	if(size_old%size_align)
		for(std::size_t i=size_old%size_align; i<size_align; ++i)
			ss<<c0;
	return exponent_align;
}


// Send str
inline void Cereal_Func::mpi_send(const std::string &str, const std::size_t exponent_align, const int rank_recv, const int tag, const MPI_Comm &mpi_comm)
{
	this->char_contiguous.resize(exponent_align);
  #if MPI_VERSION>=4
	MPI_CHECK( MPI_Send_c( str.c_str(), str.size()>>exponent_align, this->char_contiguous(exponent_align), rank_recv, tag, mpi_comm ) );
  #else
	MPI_CHECK( MPI_Send  ( str.c_str(), str.size()>>exponent_align, this->char_contiguous(exponent_align), rank_recv, tag, mpi_comm ) );
  #endif
}

// Send data
template<typename... Ts>
void Cereal_Func::mpi_send(const int rank_recv, const int tag, const MPI_Comm &mpi_comm,
	const Ts&... data)
{
	std::stringstream ss;
	{
		cereal::BinaryOutputArchive ar(ss);
		ar(data...);
	}
	const std::size_t exponent_align = align_stringstream(ss);
	const std::string &str = ss.str();
	mpi_send(str, exponent_align, rank_recv, tag, mpi_comm);
}


// Isend str
inline void Cereal_Func::mpi_isend(const std::string &str, const std::size_t exponent_align, const int rank_recv, const int tag, const MPI_Comm &mpi_comm, MPI_Request &request)
{
	this->char_contiguous.resize(exponent_align);
  #if MPI_VERSION>=4
	MPI_CHECK( MPI_Isend_c( str.c_str(), str.size()>>exponent_align, this->char_contiguous(exponent_align), rank_recv, tag, mpi_comm, &request ) );
  #else
	MPI_CHECK( MPI_Isend  ( str.c_str(), str.size()>>exponent_align, this->char_contiguous(exponent_align), rank_recv, tag, mpi_comm, &request ) );
  #endif
}

// Isend data using temporary memory str
template<typename... Ts>
void Cereal_Func::mpi_isend(const int rank_recv, const int tag, const MPI_Comm &mpi_comm,
	std::string &str, MPI_Request &request,
	const Ts&... data)
{
	std::stringstream ss;
	{
		cereal::BinaryOutputArchive ar(ss);
		ar(data...);
	}
	const std::size_t exponent_align = align_stringstream(ss);
	str = ss.str();
	mpi_isend(str, exponent_align, rank_recv, tag, mpi_comm, request);
}


// Recv to return
inline std::vector<char> Cereal_Func::mpi_recv(const MPI_Comm &mpi_comm, MPI_Status &status)
{
	for(std::size_t exponent_align=0; ; ++exponent_align)
	{
		this->char_contiguous.resize(exponent_align);
		const MPI_Datatype mpi_type = this->char_contiguous(exponent_align);
	  #if MPI_VERSION>=4
		MPI_Count size;		MPI_CHECK( MPI_Get_count_c( &status, mpi_type, &size ) );
	  #else
		int size;			MPI_CHECK( MPI_Get_count  ( &status, mpi_type, &size ) );
	  #endif
		if(size!=MPI_UNDEFINED)
		{
			std::vector<char> c(std::size_t(size)<<exponent_align);
		  #if MPI_VERSION>=4
			MPI_CHECK( MPI_Recv_c( c.data(), size, mpi_type, status.MPI_SOURCE, status.MPI_TAG, mpi_comm, MPI_STATUS_IGNORE ) );
		  #else
			MPI_CHECK( MPI_Recv  ( c.data(), size, mpi_type, status.MPI_SOURCE, status.MPI_TAG, mpi_comm, MPI_STATUS_IGNORE ) );
		  #endif
			return c;
		}
	}
}

// Recv to data
template<typename... Ts>
MPI_Status Cereal_Func::mpi_recv(const MPI_Comm &mpi_comm,
	Ts&... data)
{
	MPI_Status status;
	MPI_CHECK( MPI_Probe( MPI_ANY_SOURCE, MPI_ANY_TAG, mpi_comm, &status ) );

	std::vector<char> c = mpi_recv(mpi_comm, status);

	std::stringstream ss;
	ss.rdbuf()->pubsetbuf(c.data(), c.size());
	{
		cereal::BinaryInputArchive ar(ss);
		ar(data...);
	}
	return status;
}


// Mrecv to return
inline std::vector<char> Cereal_Func::mpi_mrecv(MPI_Message &message_recv, const MPI_Status &status)
{
	for(std::size_t exponent_align=0; ; ++exponent_align)
	{
		this->char_contiguous.resize(exponent_align);
		const MPI_Datatype mpi_type = this->char_contiguous(exponent_align);
	  #if MPI_VERSION>=4
		MPI_Count size;		MPI_CHECK( MPI_Get_count_c( &status, mpi_type, &size ) );
	  #else
		int size;			MPI_CHECK( MPI_Get_count  ( &status, mpi_type, &size ) );
	  #endif
		if(size!=MPI_UNDEFINED)
		{
			std::vector<char> c(std::size_t(size)<<exponent_align);
		  #if MPI_VERSION>=4
			MPI_CHECK( MPI_Mrecv_c( c.data(), size, mpi_type, &message_recv, MPI_STATUS_IGNORE ) );
		  #else
			MPI_CHECK( MPI_Mrecv  ( c.data(), size, mpi_type, &message_recv, MPI_STATUS_IGNORE ) );
		  #endif
			return c;
		}
	}
}

}

#undef MPI_CHECK