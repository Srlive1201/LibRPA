// ===================
//  Author: Peize Lin
//  date: 2022.05.01
// ===================

#pragma once

#include "MPI_Wrapper.h"

#include <mpi.h>
#include <sstream>
#include <vector>

namespace Comm
{

class Cereal_Func
{
  public:

	inline Cereal_Func();

	// every 2^exponent_align char concatenate to 1 word
	inline std::size_t align_stringstream(std::stringstream &ss);

	// Send str
	inline void mpi_send(const std::string &str, const std::size_t exponent_align, const int rank_recv, const int tag, const MPI_Comm &mpi_comm);

	// Send data
	template<typename... Ts>
	void mpi_send(const int rank_recv, const int tag, const MPI_Comm &mpi_comm,
		const Ts&... data);

	// Isend str
	inline void mpi_isend(const std::string &str, const std::size_t exponent_align, const int rank_recv, const int tag, const MPI_Comm &mpi_comm, MPI_Request &request);

	// Isend data using temporary memory str
	template<typename... Ts>
	void mpi_isend(const int rank_recv, const int tag, const MPI_Comm &mpi_comm,
		std::string &str, MPI_Request &request,
		const Ts&... data);

	// Recv to return
	inline std::vector<char> mpi_recv(const MPI_Comm &mpi_comm, MPI_Status &status);

	// Recv to data
	template<typename... Ts>
	MPI_Status mpi_recv(const MPI_Comm &mpi_comm,
		Ts&... data);

	// Mrecv to return
	inline std::vector<char> mpi_mrecv(MPI_Message &message_recv, const MPI_Status &status);

  private:

	MPI_Wrapper::MPI_Type_Contiguous_Pool char_contiguous{MPI_CHAR};
};

}

#include "Cereal_Func.hpp"