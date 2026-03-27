// ===================
//  Author: Peize Lin
//  date: 2023.06.08
// ===================

#include <mpi.h>

#pragma once

namespace RI
{

namespace MPI_Wrapper
{
	extern inline int mpi_get_rank(const MPI_Comm &mpi_comm);
	extern inline int mpi_get_size(const MPI_Comm &mpi_comm);

	template<typename T> extern inline MPI_Datatype mpi_get_datatype(const T&v);

	template<typename T> extern inline void mpi_reduce(T &data, const MPI_Op &op, const int &root, const MPI_Comm &mpi_comm);
	template<typename T> extern inline void mpi_allreduce(T &data, const MPI_Op &op, const MPI_Comm &mpi_comm);
	template<typename T> extern inline void mpi_reduce(T*const ptr, const int &count, const MPI_Op &op, const int &root, const MPI_Comm &mpi_comm);
	template<typename T> extern inline void mpi_allreduce(T*const ptr, const int &count, const MPI_Op &op, const MPI_Comm &mpi_comm);

	class mpi_comm;
}

}

#include "MPI_Wrapper-func.h"
#include "MPI_Wrapper-class.h"