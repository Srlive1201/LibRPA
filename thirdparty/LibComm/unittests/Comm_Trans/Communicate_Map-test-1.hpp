//=======================
// AUTHOR : Peize Lin
// DATE :   2022-01-05
//=======================

#pragma once

#include <mpi.h>
#include <map>
#include <tuple>
#include <fstream>
#include "Comm/Comm_Trans/Comm_Trans.h"
#include "Comm/example/Communicate_Map-1.h"
#include "Comm/global/MPI_Wrapper.h"
#include "unittests/print_stl.h"

namespace Communicate_Map_Test
{
	static void test_transmission(int argc, char *argv[])
	{
		int provided;
		MPI_Init_thread( &argc, &argv, MPI_THREAD_MULTIPLE, &provided );

		std::map<int,std::map<int,double>> m_in;
		std::map<int,std::map<int,double>> m_out;

		if(Comm::MPI_Wrapper::mpi_get_rank(MPI_COMM_WORLD)==0){ m_in[0][0]=1.2;	m_in[2][1]=3.4;}
		else if(Comm::MPI_Wrapper::mpi_get_rank(MPI_COMM_WORLD)==1){ m_in[2][0]=5.6; }

		Comm::Comm_Trans<
			std::tuple<int,int>,
			double,
			std::map<int,std::map<int,double>>,
			std::map<int,std::map<int,double>>
		> com(MPI_COMM_WORLD);

//		com.set_value_recv = Comm::Communicate_Map::set_value_assignment<int,int,double>;
		com.set_value_recv = Comm::Communicate_Map::set_value_add<int,int,double>;

//		com.traverse_isend = Communicate_Map::traverse_datas_all<int,int,double>;
		com.traverse_isend = std::bind(
			Comm::Communicate_Map::traverse_datas_mask<int,int,double>,
			std::placeholders::_1, std::placeholders::_2, std::placeholders::_3,
			[&com](const int rank_isend, const int &key0) ->bool { return (key0%com.comm_size==rank_isend) ? true : false; },
			[](const int rank_isend, const int &key1) ->bool { return true; });

		com.flag_lock_set_value = Comm::Comm_Tools::Lock_Type::Copy_merge;
		com.init_datas_local = Comm::Communicate_Map::init_datas_local<int,int,double>;
		com.add_datas = Comm::Communicate_Map::add_datas<int,int,double>;

		com.communicate(m_in,m_out);

		std::ofstream ofs_in("in_"+std::to_string(Comm::MPI_Wrapper::mpi_get_rank(MPI_COMM_WORLD)));
		ofs_in<<m_in<<std::endl;
		std::ofstream ofs_out("out_"+std::to_string(Comm::MPI_Wrapper::mpi_get_rank(MPI_COMM_WORLD)));
		ofs_out<<m_out<<std::endl;

		MPI_Finalize();
	}
	/*
		mpirun -n 2
		rank 0:
		{
			0: {0: 1.2},
			2: {0: 5.6,
				1: 3.4}
		}
		rank 1: {}

		mpirun -n 3
		rank 0:
		{
			0: {0: 1.2}
		}
		rank 1: {}
		rank 2:
		{
			2: {0: 5.6,
				1: 3.4}
		}
	*/
}