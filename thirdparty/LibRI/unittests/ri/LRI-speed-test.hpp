// ===================
//  Author: Peize Lin
//  date: 2022.09.21
// ===================

#pragma once

#include"RI/ri/Label.h"
#include"RI/ri/Label_Tools.h"
#include"RI/ri/LRI.h"
#include"RI/global/MPI_Wrapper.h"
#include"RI/global/Global_Func-1.h"

#include<array>
#include<map>
#include<unordered_map>
#include<iostream>
#include<mpi.h>
#include<sys/time.h>

namespace LRI_Speed_Test
{
	template<typename Tdata>
	RI::Tensor<Tdata> init_tensor(const RI::Shape_Vector &shape)
	{
		RI::Tensor<Tdata> D(shape);
		for(std::size_t i=0; i<D.data->size(); ++i)
			(*D.data)[i] = i;
		return D;
	}

	inline double time_during(const timeval &t_begin)
	{
		timeval t_end;
		gettimeofday(&t_end, NULL);
		return (double)(t_end.tv_sec-t_begin.tv_sec) + (double)(t_end.tv_usec-t_begin.tv_usec)/1000000.0;
	}


	template<typename Tdata>
	void test_speed(int argc, char *argv[], const int NA, const std::size_t Ni)
	{

		int mpi_init_provide;
		MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &mpi_init_provide);

		constexpr std::size_t Ndim = 1;
		//const std::size_t Na0=20, Nb0=30, Na1=40, Nb1=50, Na2=60, Nb2=70;
		const std::size_t Na0=Ni, Nb0=Ni, Na1=Ni, Nb1=Ni, Na2=Ni, Nb2=Ni;

		using T_Ds = std::map<int, std::map<std::pair<int,std::array<int,Ndim>>, RI::Tensor<Tdata>>>;
		std::unordered_map<RI::Label::ab, T_Ds> Ds_ab;
		Ds_ab.reserve(11);
		auto init_Ds = [&Ds_ab, NA](const RI::Label::ab &label, const RI::Shape_Vector &shape)
		{
			const int rank_mine = RI::MPI_Wrapper::mpi_get_rank(MPI_COMM_WORLD);
			const int rank_size = RI::MPI_Wrapper::mpi_get_size(MPI_COMM_WORLD);
			const RI::Tensor<Tdata> D = init_tensor<Tdata>(shape);
			for(int iAx=0; iAx<NA; ++iAx)
			{
				if(iAx%rank_size!=rank_mine)	continue;
				for(int iAy=0; iAy<NA; ++iAy)
					Ds_ab[label][iAx][{iAy,{0}}] = D;
			}
		};
		init_Ds(RI::Label::ab::a, {Na0,Na1,Na2});
		init_Ds(RI::Label::ab::b, {Nb0,Nb1,Nb2});
		init_Ds(RI::Label::ab::a0b0, {Na0,Nb0});
		init_Ds(RI::Label::ab::a0b1, {Na0,Nb1});
		init_Ds(RI::Label::ab::a0b2, {Na0,Nb2});
		init_Ds(RI::Label::ab::a1b0, {Na1,Nb0});
		init_Ds(RI::Label::ab::a1b1, {Na1,Nb1});
		init_Ds(RI::Label::ab::a1b2, {Na1,Nb2});
		init_Ds(RI::Label::ab::a2b0, {Na2,Nb0});
		init_Ds(RI::Label::ab::a2b1, {Na2,Nb1});
		init_Ds(RI::Label::ab::a2b2, {Na2,Nb2});

		std::map<int,std::array<double,1>> atoms_pos;
		for(int iA=0; iA<NA; ++iA)
			atoms_pos[iA] = {0};

		RI::LRI<int,int,Ndim,Tdata> lri;
		lri.set_parallel(MPI_COMM_WORLD, atoms_pos, {}, {1}, RI::Global_Func::to_vector(RI::Label::array_ab_ab));

		for(const RI::Label::ab &label : RI::Label::array_ab)
			lri.set_tensors_map2(Ds_ab[label], {label});

		timeval t_begin;
		gettimeofday(&t_begin, NULL);

		T_Ds Ds_result;
		lri.cal_loop3(RI::Global_Func::to_vector(RI::Label::array_ab_ab), Ds_result);

		std::cout<<time_during(t_begin)<<std::endl;

		MPI_Finalize();
	}

}