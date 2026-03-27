// ===================
//  Author: Peize Lin
//  date: 2022.06.02
// ===================

#pragma once

#include <mpi.h>
#include <complex>
#include <vector>
#include <string>
#include <stdexcept>

#define MPI_CHECK(x) if((x)!=MPI_SUCCESS)	throw std::runtime_error(std::string(__FILE__)+" line "+std::to_string(__LINE__));

namespace RI
{

namespace MPI_Wrapper
{
	inline int mpi_get_rank(const MPI_Comm &mpi_comm)
	{
		int rank_mine;
		MPI_CHECK( MPI_Comm_rank (mpi_comm, &rank_mine) );
		return rank_mine;
	}

	inline int mpi_get_size(const MPI_Comm &mpi_comm)
	{
		int rank_size;
		MPI_CHECK( MPI_Comm_size (mpi_comm, &rank_size) );
		return rank_size;
	}

	inline MPI_Datatype mpi_get_datatype(const char							&v) { return MPI_CHAR; }
	inline MPI_Datatype mpi_get_datatype(const short						&v) { return MPI_SHORT; }
	inline MPI_Datatype mpi_get_datatype(const int							&v) { return MPI_INT; }
	inline MPI_Datatype mpi_get_datatype(const long							&v) { return MPI_LONG; }
	inline MPI_Datatype mpi_get_datatype(const long long					&v) { return MPI_LONG_LONG; }
	inline MPI_Datatype mpi_get_datatype(const unsigned char				&v) { return MPI_UNSIGNED_CHAR; }
	inline MPI_Datatype mpi_get_datatype(const unsigned short				&v) { return MPI_UNSIGNED_SHORT; }
	inline MPI_Datatype mpi_get_datatype(const unsigned int					&v) { return MPI_UNSIGNED; }
	inline MPI_Datatype mpi_get_datatype(const unsigned long				&v) { return MPI_UNSIGNED_LONG; }
	inline MPI_Datatype mpi_get_datatype(const unsigned long long			&v) { return MPI_UNSIGNED_LONG_LONG; }
	inline MPI_Datatype mpi_get_datatype(const float						&v) { return MPI_FLOAT; }
	inline MPI_Datatype mpi_get_datatype(const double						&v) { return MPI_DOUBLE; }
	inline MPI_Datatype mpi_get_datatype(const long double					&v) { return MPI_LONG_DOUBLE; }
	inline MPI_Datatype mpi_get_datatype(const bool							&v) { return MPI_CXX_BOOL; }
	inline MPI_Datatype mpi_get_datatype(const std::complex<float>			&v) { return MPI_CXX_FLOAT_COMPLEX; }
	inline MPI_Datatype mpi_get_datatype(const std::complex<double>			&v) { return MPI_CXX_DOUBLE_COMPLEX; }
	inline MPI_Datatype mpi_get_datatype(const std::complex<long double>	&v) { return MPI_CXX_LONG_DOUBLE_COMPLEX; }

	//inline int mpi_get_count(const MPI_Status &status, const MPI_Datatype &datatype)
	//{
	//	int count;
	//	MPI_CHECK( MPI_Get_count(&status, datatype, &count) );
	//	return count;
	//}

	template<typename T>
	inline void mpi_reduce(T &data, const MPI_Op &op, const int &root, const MPI_Comm &mpi_comm)
	{
		T data_out;
		MPI_CHECK( MPI_Reduce(&data, &data_out, 1, mpi_get_datatype(data), op, root, mpi_comm) );
		if(mpi_get_rank(mpi_comm)==root)
			data = data_out;
	}
	template<typename T>
	inline void mpi_allreduce(T &data, const MPI_Op &op, const MPI_Comm &mpi_comm)
	{
		T data_out;
		MPI_CHECK( MPI_Allreduce(&data, &data_out, 1, mpi_get_datatype(data), op, mpi_comm) );
		data = data_out;
	}

	template<typename T>
	inline void mpi_reduce(T*const ptr, const int &count, const MPI_Op &op, const int &root, const MPI_Comm &mpi_comm)
	{
		std::vector<T> ptr_out(count);
		MPI_CHECK( MPI_Reduce(ptr, ptr_out.data(), count, mpi_get_datatype(*ptr), op, root, mpi_comm) );
		if(mpi_get_rank(mpi_comm)==root)
			for(std::size_t i=0; i<count; ++i)
				ptr[i] = ptr_out[i];
	}
	template<typename T>
	inline void mpi_allreduce(T*const ptr, const int &count, const MPI_Op &op, const MPI_Comm &mpi_comm)
	{
		std::vector<T> ptr_out(count);
		MPI_CHECK( MPI_Allreduce(ptr, ptr_out.data(), count, mpi_get_datatype(*ptr), op, mpi_comm) );
		for(std::size_t i=0; i<count; ++i)
			ptr[i] = ptr_out[i];
	}
}

}

#undef MPI_CHECK