// ===================
//  Author: Minye Zhang, almost copied from Exx-test.hpp
//  date: 2022.06.02
// ===================

#pragma once

#include "RI/physics/GW.h"
// #include <complex>

namespace GW_Test
{
	template<typename Tdata>
	void main(int argc, char *argv[])
	{
		int mpi_init_provide;
		MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &mpi_init_provide);

		RI::GW<int,int,1,Tdata> gw;
		gw.set_parallel(MPI_COMM_WORLD, {{1,{0}},{2,{4}}}, {}, {1});
		gw.set_Cs({}, 1E-4);
		gw.set_Ws({}, 1E-4);
		gw.set_Gs({}, 1E-4);
		gw.cal_Sigmas();

		gw.free_Cs();
		gw.free_Ws();
		gw.free_Gs();

		MPI_Finalize();
	}
}