//=======================
// AUTHOR : Peize Lin
// DATE :   2022-07-06
//=======================

#pragma once

#include "Comm/Comm_Assemble/Comm_Assemble.h"
#include "Comm/example/Communicate_Map-1.h"
#include "Comm/example/Communicate_Map-2.h"
#include "Comm/global/MPI_Wrapper.h"
#include "unittests/print_stl.h"

#include <mpi.h>
#include <map>
#include <tuple>
#include <fstream>
#include <cassert>
#include <stdexcept>
#include <string>

#define MPI_CHECK(x) if((x)!=MPI_SUCCESS)	throw std::runtime_error(std::string(__FILE__)+" line "+std::to_string(__LINE__));

namespace Communicate_Map_Test
{
	static void test_assemble(int argc, char *argv[])
	{
		int provided;
		MPI_CHECK( MPI_Init_thread( &argc, &argv, MPI_THREAD_MULTIPLE, &provided ) );
		assert(Comm::MPI_Wrapper::mpi_get_size(MPI_COMM_WORLD)==6);
		const int rank_mine = Comm::MPI_Wrapper::mpi_get_rank(MPI_COMM_WORLD);

		std::map<int,std::map<int,double>> m_in;
		std::map<int,std::map<int,double>> m_out;

		if(Comm::MPI_Wrapper::mpi_get_rank(MPI_COMM_WORLD)==0){ m_in[0][0]=1.2;	m_in[2][1]=3.4;}
		else if(Comm::MPI_Wrapper::mpi_get_rank(MPI_COMM_WORLD)==1){ m_in[2][0]=5.6; }
		else if(Comm::MPI_Wrapper::mpi_get_rank(MPI_COMM_WORLD)==3){ m_in[2][1]=7.8;	m_in[1][10]=9.0; }

		Comm::Communicate_Map::Judge_Map2<int,int> judge;
		judge.s0 = {rank_mine%3};
		judge.s1 = {rank_mine/3};
		/*
			s0\s1	0	1
			0		0	3
			1		1	4
			2		2	5
		*/
		//judge.s1 = {0,1,2};
		/*
			s0
			0	03
			1	14
			2	25
		*/

		Comm::Comm_Assemble<
			std::tuple<int,int>,
			double,
			std::map<int,std::map<int,double>>,
			Comm::Communicate_Map::Judge_Map2<int,int>,
			std::map<int,std::map<int,double>>
		> com(MPI_COMM_WORLD);

		com.traverse_keys_provide = Comm::Communicate_Map::traverse_keys<int,int,double>;
		com.get_value_provide = Comm::Communicate_Map::get_value<int,int,double>;
		com.set_value_require = Comm::Communicate_Map::set_value_add<int,int,double>;
		com.flag_lock_set_value = Comm::Comm_Tools::Lock_Type::Copy_merge;
		com.init_datas_local = Comm::Communicate_Map::init_datas_local<int,int,double>;
		com.add_datas = Comm::Communicate_Map::add_datas<int,int,double>;

		com.communicate( m_in, judge, m_out );

		std::ofstream ofs_in("in_"+std::to_string(Comm::MPI_Wrapper::mpi_get_rank(MPI_COMM_WORLD)));
		ofs_in<<m_in<<std::endl;
		std::ofstream ofs_out("out_"+std::to_string(Comm::MPI_Wrapper::mpi_get_rank(MPI_COMM_WORLD)));
		ofs_out<<m_out<<std::endl;

		MPI_CHECK( MPI_Finalize() );
	}
	/*
	judge.s1 = {rank_mine/3};
	mpirun -n 6
		rank 0:	[0][0] = 1.2
		rank 1:
		rank 2:	[2][0] = 5.6
		rank 3:
		rank 4:
		rank 5:	[2][1] = 11.2
	*/
	/*
	judge.s1 = {0,1,2};
	mpirun -n 6
		rank 0:	[0][0]=1.2
		rank 1:
		rank 2:	[2][0]=5.6, [2][1]=11.2
		rank 3:	[0][0]=1.2
		rank 4:
		rank 5:	[2][0]=5.6, [2][1]=11.2
	*/
}

#undef MPI_CHECK
