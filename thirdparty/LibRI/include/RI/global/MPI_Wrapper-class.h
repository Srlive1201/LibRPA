// ===================
//  Author: Peize Lin
//  date: 2023.06.08
// ===================

#pragma once

#include <mpi.h>
#include <stdexcept>
#include <string>

#define MPI_CHECK(x) if((x)!=MPI_SUCCESS)	throw std::runtime_error(std::string(__FILE__)+" line "+std::to_string(__LINE__));

namespace RI
{

namespace MPI_Wrapper
{
	class mpi_comm
	{
	  public:
		MPI_Comm comm;
		bool flag_allocate = false;								// flag_allocate=true is controled by user

		mpi_comm() = default;
		mpi_comm(const MPI_Comm &comm_in, const bool &flag_allocate_in): comm(comm_in), flag_allocate(flag_allocate_in) {}
		mpi_comm(const mpi_comm &mc_in) = delete;
		mpi_comm(mpi_comm &mc_in) = delete;
		mpi_comm(mpi_comm &&mc_in)
		{
			this->free();
			this->comm = mc_in.comm;
			this->flag_allocate = mc_in.flag_allocate;
			mc_in.flag_allocate = false;
		}
		mpi_comm &operator=(const mpi_comm &mc_in) = delete;
		mpi_comm &operator=(mpi_comm &mc_in) = delete;
		mpi_comm &operator=(mpi_comm &&mc_in)
		{
			this->free();
			this->comm = mc_in.comm;
			this->flag_allocate = mc_in.flag_allocate;
			mc_in.flag_allocate = false;
			return *this;
		}

		~mpi_comm() { this->free(); }

		MPI_Comm &operator()(){ return this->comm; }
		const MPI_Comm &operator()()const{ return this->comm; }

		void free()
		{
			if(this->flag_allocate)
			{
				MPI_CHECK( MPI_Comm_free( &this->comm ) );
				this->flag_allocate = false;
			}
		}
	};
}

}

#undef MPI_CHECK