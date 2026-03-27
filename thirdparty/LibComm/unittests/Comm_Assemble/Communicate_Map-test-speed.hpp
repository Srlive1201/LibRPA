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
#include <valarray>
#include <map>
#include <tuple>
#include <fstream>
#include <stdexcept>
#include <string>

#include <omp.h>
#include <sys/time.h>
#include <cereal/types/valarray.hpp>

#define MPI_CHECK(x) if((x)!=MPI_SUCCESS)	throw std::runtime_error(std::string(__FILE__)+" line "+std::to_string(__LINE__));

namespace Communicate_Map_Test
{
	static void test_speed(int argc, char *argv[])
	{
		int provided;
		MPI_CHECK( MPI_Init_thread( &argc, &argv, MPI_THREAD_MULTIPLE, &provided ) );
		const int rank_mine = Comm::MPI_Wrapper::mpi_get_rank(MPI_COMM_WORLD);
		const int rank_size = Comm::MPI_Wrapper::mpi_get_size(MPI_COMM_WORLD);

		std::map<int,std::valarray<double>> m_in;
		std::map<int,std::valarray<double>> m_out;

		const int N = std::stoi(argv[1]);
		const int M = std::stoi(argv[2]);

		const int misson_size = std::floor(N/rank_size);
		for(int i=0; i<misson_size; ++i)
			m_in[rank_mine*misson_size+i] = std::valarray<double>(M);

		Comm::Communicate_Map::Judge_Map<int> judge;
		for(int i=0; i<N; ++i)
			if(i%rank_size==rank_mine)
				judge.s.insert(i);

		Comm::Comm_Assemble<
			int,
			std::valarray<double>,
			std::map<int,std::valarray<double>>,
			Comm::Communicate_Map::Judge_Map<int>,
			std::map<int,std::valarray<double>>
		> com(MPI_COMM_WORLD);

		com.traverse_keys_provide = Comm::Communicate_Map::traverse_keys<int,std::valarray<double>>;
		com.get_value_provide = Comm::Communicate_Map::get_value<int,std::valarray<double>>;
		com.set_value_require = Comm::Communicate_Map::set_value_add<int,std::valarray<double>>;
		com.flag_lock_set_value = Comm::Comm_Tools::Lock_Type::Copy_merge;
		com.init_datas_local = Comm::Communicate_Map::init_datas_local<int,std::valarray<double>>;
		com.add_datas = Comm::Communicate_Map::add_datas<int,std::valarray<double>>;

			MPI_Barrier(MPI_COMM_WORLD);
			timeval t_begin;	gettimeofday( &t_begin, NULL);
		com.communicate( m_in, judge, m_out );
			MPI_Barrier(MPI_COMM_WORLD);
			timeval t_end;		gettimeofday( &t_end, NULL);
			if(rank_mine==0)
				std::cout<<N<<"\t"<<M<<"\t"<<rank_size<<"\t"<<omp_get_max_threads()<<"\t"
			             <<(double)(t_end.tv_sec-t_begin.tv_sec) + (double)(t_end.tv_usec-t_begin.tv_usec)/1000000.0<<std::endl;

		MPI_CHECK( MPI_Finalize() );
	}
}

#undef MPI_CHECK
