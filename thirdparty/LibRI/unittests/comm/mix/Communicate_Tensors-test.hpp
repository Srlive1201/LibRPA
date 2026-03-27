//=======================
// AUTHOR : Peize Lin
// DATE :   2022-07-20
//=======================

#pragma once

#include "RI/comm/mix/Communicate_Tensors_Map_Judge.h"
#include "RI/global/MPI_Wrapper.h"
#include "unittests/global/Tensor-test.h"
#include "unittests/print_stl.h"

#include <mpi.h>
#include <map>
#include <set>
#include <fstream>

namespace Communicate_Tensors_Test
{
	inline RI::Tensor<double> init_Tensor(const double d)
	{
		RI::Tensor<double> t({1});
		t(0) = d;
		return t;
	}

	/*
	mpirun -n 5
		rank 0:
			2	22
			3	303
			5	5
		rank 1:
			3	303
			4	44
			5	5
			6	660
		rank 2:
		rank 3:
			1	1
			3	303
			4	44
		rank 4:
	*/

	static void test_comm_judge_map(int argc, char *argv[])
	{
		int mpi_init_provide;	MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &mpi_init_provide);
		const int rank_mine = RI::MPI_Wrapper::mpi_get_rank(MPI_COMM_WORLD);
		std::map<int,RI::Tensor<double>> Ds_in;
		std::set<int> s;
		if(rank_mine==0)
		{
			for(int i=0; i<6; ++i)
				Ds_in[i]=init_Tensor(i);
			s={2,3,5,7,11,13};
		}
		else if(rank_mine==1)
		{
			for(int i=0; i<10; i+=2)
				Ds_in[i]=init_Tensor(10*i);
			s={7,6,5,4,3};
		}
		else if(rank_mine==2)
		{
			for(int i=0; i<10; i+=3)
				Ds_in[i]=init_Tensor(100*i);
		}
		else if(rank_mine==3)
		{
			s={3,1,4};
		}
		std::map<int,RI::Tensor<double>> Ds_out = RI::Communicate_Tensors_Map_Judge::comm_map(MPI_COMM_WORLD, Ds_in, s);
		std::ofstream ofs("out."+std::to_string(rank_mine));
		ofs<<Ds_out<<std::endl;
		MPI_Finalize();
	}

	static void test_comm_judge_map2(int argc, char *argv[])
	{
		int mpi_init_provide;	MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &mpi_init_provide);
		const int rank_mine = RI::MPI_Wrapper::mpi_get_rank(MPI_COMM_WORLD);
		using TAC = std::pair<int,std::array<int,1>>;
		std::map<int,std::map<TAC,RI::Tensor<double>>> Ds_in;
		std::set<int> s;
		std::set<TAC> s1;
		for(int i=0; i<10; ++i)
			s1.insert({i,{i}});
		if(rank_mine==0)
		{
			for(int i=0; i<6; ++i)
				Ds_in[i][{i,{i}}]=init_Tensor(i);
			s={2,3,5,7,11,13};
		}
		else if(rank_mine==1)
		{
			for(int i=0; i<10; i+=2)
				Ds_in[i][{i,{i}}]=init_Tensor(10*i);
			s={7,6,5,4,3};
		}
		else if(rank_mine==2)
		{
			for(int i=0; i<10; i+=3)
				Ds_in[i][{i,{i}}]=init_Tensor(100*i);
		}
		else if(rank_mine==3)
		{
			s={3,1,4};
		}
		std::map<int,std::map<TAC,RI::Tensor<double>>> Ds_out = RI::Communicate_Tensors_Map_Judge::comm_map2(MPI_COMM_WORLD, Ds_in, s, s1);
		std::ofstream ofs("out."+std::to_string(rank_mine));
		ofs<<Ds_out<<std::endl;
		MPI_Finalize();
	}

	static void test_comm_judge_map3(int argc, char *argv[])
	{
		int mpi_init_provide;	MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &mpi_init_provide);
		const int rank_mine = RI::MPI_Wrapper::mpi_get_rank(MPI_COMM_WORLD);
		using TAC = std::pair<int,std::array<int,1>>;
		std::map<int,std::map<TAC,std::map<TAC,RI::Tensor<double>>>> Ds_in;
		std::set<int> s;
		std::set<TAC> s1;
		for(int i=0; i<10; ++i)
			s1.insert({i,{i}});
		std::set<TAC> s2={{0,{0}}};
		if(rank_mine==0)
		{
			for(int i=0; i<6; ++i)
				Ds_in[i][{i,{i}}][{0,{0}}]=init_Tensor(i);
			s={2,3,5,7,11,13};
		}
		else if(rank_mine==1)
		{
			for(int i=0; i<10; i+=2)
				Ds_in[i][{i,{i}}][{0,{0}}]=init_Tensor(10*i);
			s={7,6,5,4,3};
		}
		else if(rank_mine==2)
		{
			for(int i=0; i<10; i+=3)
				Ds_in[i][{i,{i}}][{0,{0}}]=init_Tensor(100*i);
		}
		else if(rank_mine==3)
		{
			s={3,1,4};
		}
		std::map<int,std::map<TAC,std::map<TAC,RI::Tensor<double>>>> Ds_out = RI::Communicate_Tensors_Map_Judge::comm_map3(MPI_COMM_WORLD, Ds_in, s, s1, s2);
		std::ofstream ofs("out."+std::to_string(rank_mine));
		ofs<<Ds_out<<std::endl;
		MPI_Finalize();
	}

	static void test_comm_judge_map2_first(int argc, char *argv[])
	{
		int mpi_init_provide;	MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &mpi_init_provide);
		const int rank_mine = RI::MPI_Wrapper::mpi_get_rank(MPI_COMM_WORLD);
		using TAC = std::pair<int,std::array<int,1>>;
		std::map<int,std::map<TAC,RI::Tensor<double>>> Ds_in;
		std::set<int> s;
		std::set<int> s1={0};
		if(rank_mine==0)
		{
			for(int i=0; i<6; ++i)
				Ds_in[i][{0,{i}}]=init_Tensor(i);
			s={2,3,5,7,11,13};
		}
		else if(rank_mine==1)
		{
			for(int i=0; i<10; i+=2)
				Ds_in[i][{0,{i}}]=init_Tensor(10*i);
			s={7,6,5,4,3};
		}
		else if(rank_mine==2)
		{
			for(int i=0; i<10; i+=3)
				Ds_in[i][{0,{i}}]=init_Tensor(100*i);
		}
		else if(rank_mine==3)
		{
			s={3,1,4};
		}
		std::map<int,std::map<TAC,RI::Tensor<double>>> Ds_out = RI::Communicate_Tensors_Map_Judge::comm_map2_first(MPI_COMM_WORLD, Ds_in, s, s1);
		std::ofstream ofs("out."+std::to_string(rank_mine));
		ofs<<Ds_out<<std::endl;
		MPI_Finalize();
	}

	static void test_comm_judge_map3_first(int argc, char *argv[])
	{
		int mpi_init_provide;	MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &mpi_init_provide);
		const int rank_mine = RI::MPI_Wrapper::mpi_get_rank(MPI_COMM_WORLD);
		using TAC = std::pair<int,std::array<int,1>>;
		std::map<int,std::map<TAC,std::map<TAC,RI::Tensor<double>>>> Ds_in;
		std::set<int> s;
		std::set<int> s1={0};
		std::set<int> s2={0};
		if(rank_mine==0)
		{
			for(int i=0; i<6; ++i)
				Ds_in[i][{0,{i}}][{0,{i}}]=init_Tensor(i);
			s={2,3,5,7,11,13};
		}
		else if(rank_mine==1)
		{
			for(int i=0; i<10; i+=2)
				Ds_in[i][{0,{i}}][{0,{i}}]=init_Tensor(10*i);
			s={7,6,5,4,3};
		}
		else if(rank_mine==2)
		{
			for(int i=0; i<10; i+=3)
				Ds_in[i][{0,{i}}][{0,{i}}]=init_Tensor(100*i);
		}
		else if(rank_mine==3)
		{
			s={3,1,4};
		}
		std::map<int,std::map<TAC,std::map<TAC,RI::Tensor<double>>>> Ds_out = RI::Communicate_Tensors_Map_Judge::comm_map3_first(MPI_COMM_WORLD, Ds_in, s, s1, s2);
		std::ofstream ofs("out."+std::to_string(rank_mine));
		ofs<<Ds_out<<std::endl;
		MPI_Finalize();
	}

	static void test_comm_judge_map2_period(int argc, char *argv[])
	{
		int mpi_init_provide;	MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &mpi_init_provide);
		const int rank_mine = RI::MPI_Wrapper::mpi_get_rank(MPI_COMM_WORLD);
		using TAC = std::pair<int,std::array<int,1>>;
		std::map<int,std::map<TAC,RI::Tensor<double>>> Ds_in;
		std::set<int> s;
		std::set<TAC> s0;
		std::set<TAC> s1;
		for(int i=0; i<10; ++i)
			s1.insert({i,{i}});
		if(rank_mine==0)
		{
			for(int i=0; i<6; ++i)
				Ds_in[i][{i,{0}}]=init_Tensor(i);
			s={2,3,5,7,11,13};
		}
		else if(rank_mine==1)
		{
			for(int i=0; i<10; i+=2)
				Ds_in[i][{i,{0}}]=init_Tensor(10*i);
			s={7,6,5,4,3};
		}
		else if(rank_mine==2)
		{
			for(int i=0; i<10; i+=3)
				Ds_in[i][{i,{0}}]=init_Tensor(100*i);
		}
		else if(rank_mine==3)
		{
			s={3,1,4};
		}
		for(const int is : s)
			s0.insert({is,{0}});
		std::map<int,std::map<TAC,RI::Tensor<double>>> Ds_out = RI::Communicate_Tensors_Map_Judge::comm_map2_period(MPI_COMM_WORLD, Ds_in, s0, s1, {1});
		std::ofstream ofs("out."+std::to_string(rank_mine));
		ofs<<Ds_out<<std::endl;
		MPI_Finalize();
	}

	static void test_comm_judge_map3_period(int argc, char *argv[])
	{
		int mpi_init_provide;	MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &mpi_init_provide);
		const int rank_mine = RI::MPI_Wrapper::mpi_get_rank(MPI_COMM_WORLD);
		using TAC = std::pair<int,std::array<int,1>>;
		std::map<int,std::map<TAC,std::map<TAC,RI::Tensor<double>>>> Ds_in;
		std::set<int> s;
		std::set<TAC> s0;
		std::set<TAC> s1;
		for(int i=0; i<10; ++i)
			s1.insert({i,{i}});
		std::set<TAC> s2={{0,{0}}};
		if(rank_mine==0)
		{
			for(int i=0; i<6; ++i)
				Ds_in[i][{i,{0}}][{0,{0}}]=init_Tensor(i);
			s={2,3,5,7,11,13};
		}
		else if(rank_mine==1)
		{
			for(int i=0; i<10; i+=2)
				Ds_in[i][{i,{0}}][{0,{0}}]=init_Tensor(10*i);
			s={7,6,5,4,3};
		}
		else if(rank_mine==2)
		{
			for(int i=0; i<10; i+=3)
				Ds_in[i][{i,{0}}][{0,{0}}]=init_Tensor(100*i);
		}
		else if(rank_mine==3)
		{
			s={3,1,4};
		}
		for(const int is : s)
			s0.insert({is,{0}});
		std::map<int,std::map<TAC,std::map<TAC,RI::Tensor<double>>>> Ds_out = RI::Communicate_Tensors_Map_Judge::comm_map3_period(MPI_COMM_WORLD, Ds_in, s0, s1, s2, {1});
		std::ofstream ofs("out."+std::to_string(rank_mine));
		ofs<<Ds_out<<std::endl;
		MPI_Finalize();
	}

}