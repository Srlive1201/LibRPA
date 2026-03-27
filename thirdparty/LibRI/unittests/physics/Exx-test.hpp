// ===================
//  Author: Peize Lin
//  date: 2022.06.02
// ===================

#pragma once

#include "RI/physics/Exx.h"
#include <complex>

namespace Exx_Test
{
	template<typename Tdata>
	void main(int argc, char *argv[])
	{
		int mpi_init_provide;
		MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &mpi_init_provide);

		RI::Exx<int,int,1,Tdata> exx;
		exx.set_parallel(MPI_COMM_WORLD, {{1,{0}},{2,{4}}}, {}, {1});
		exx.set_symmetry(false, {});
	
		exx.set_Cs({}, 1E-4);
		exx.set_Vs({}, 1E-4);
		exx.set_Ds({}, 1E-4);
		exx.cal_Hs();

		exx.set_dCs({}, 1E-4);
		exx.set_dVs({}, 1E-4);
		exx.cal_force();

		exx.set_dCRs({}, 1E-4);
		exx.set_dVRs({}, 1E-4);
		exx.cal_stress();

		exx.set_Ds_delta({}, 1E-4);
		exx.cal_Hs();

		exx.free_Cs();
		exx.free_Vs();
		exx.free_Ds();
		exx.free_Ds_delta();
		exx.free_dCs();
		exx.free_dVs();
		exx.free_dCRs();
		exx.free_dVRs();

		MPI_Finalize();
	}
}