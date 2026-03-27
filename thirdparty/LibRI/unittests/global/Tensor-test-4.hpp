// ===================
//  Author: Peize Lin
//  date: 2023.02.24
// ===================

#include "Tensor-test.h"
#include "Tensor-test-2.hpp"
#include "../include/RI/global/Cereal_Types.h"

#include <Comm/global/Cereal_Func.h>
#include <Comm/global/MPI_Wrapper.h>

#include <mpi.h>
#include <iostream>

namespace Tensor_Test
{
	static void test_cereal(int argc, char *argv[])
	{
		int mpi_init_provide;
		MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &mpi_init_provide);

		assert(Comm::MPI_Wrapper::mpi_get_size(MPI_COMM_WORLD)>=2);
		const int rank_mine = Comm::MPI_Wrapper::mpi_get_rank(MPI_COMM_WORLD);

		if(rank_mine==0)
		{
			const RI::Tensor<double> m = Tensor_Test::init_real_1<double>();
			Comm::Cereal_Func::mpi_send(1, 0, MPI_COMM_WORLD, m);
			std::cout<<m<<std::endl;
		}
		else if(rank_mine==1)
		{
			RI::Tensor<double> m;
			Comm::Cereal_Func::mpi_recv(MPI_COMM_WORLD, m);
			std::cout<<m<<std::endl;
		}

		MPI_Finalize();
	}
}