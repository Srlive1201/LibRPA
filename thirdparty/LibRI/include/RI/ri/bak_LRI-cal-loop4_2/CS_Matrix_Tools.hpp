// ===================
//  Author: Peize Lin
//  date: 2022.06.06
// ===================

#pragma once

#include "CS_Matrix_Tools.h"
#include "../global/Blas_Interface-Tensor.h"
#include "../global/Tensor_Wrapper.h"
#include <stdexcept>
#include <memory.h>

namespace RI
{

namespace CS_Matrix_Tools
{
	// {Dfirst label, Dsecond label}。顺序必须与D4=DDDD的顺序相同，与实现顺序直接耦合
	std::pair<Label::ab, Label::ab> split_label(const Label::ab_ab &label)
	{
		switch(label)
		{
			case Label::ab_ab::a0b0_a1b1:	return {Label::ab::a0b0, Label::ab::a1b1};
			case Label::ab_ab::a0b1_a1b0:	return {Label::ab::a1b0, Label::ab::a0b1};
			case Label::ab_ab::a0b0_a2b1:	return {Label::ab::a0b0, Label::ab::a2b1};
			case Label::ab_ab::a0b1_a2b0:	return {Label::ab::a2b0, Label::ab::a0b1};
			case Label::ab_ab::a1b0_a2b1:	return {Label::ab::a2b1, Label::ab::a1b0};
			case Label::ab_ab::a1b1_a2b0:	return {Label::ab::a2b0, Label::ab::a1b1};
			case Label::ab_ab::a0b0_a1b2:	return {Label::ab::a0b0, Label::ab::a1b2};
			case Label::ab_ab::a0b1_a1b2:	return {Label::ab::a0b1, Label::ab::a1b2};
			case Label::ab_ab::a0b2_a1b0:	return {Label::ab::a1b0, Label::ab::a0b2};
			case Label::ab_ab::a0b2_a1b1:	return {Label::ab::a1b1, Label::ab::a0b2};
			case Label::ab_ab::a0b0_a2b2:	return {Label::ab::a0b0, Label::ab::a2b2};
			case Label::ab_ab::a0b1_a2b2:	return {Label::ab::a0b1, Label::ab::a2b2};
			case Label::ab_ab::a0b2_a2b0:	return {Label::ab::a2b0, Label::ab::a0b2};
			case Label::ab_ab::a0b2_a2b1:	return {Label::ab::a2b1, Label::ab::a0b2};
			case Label::ab_ab::a1b0_a2b2:	return {Label::ab::a1b0, Label::ab::a2b2};
			case Label::ab_ab::a1b1_a2b2:	return {Label::ab::a1b1, Label::ab::a2b2};
			case Label::ab_ab::a1b2_a2b0:	return {Label::ab::a2b0, Label::ab::a1b2};
			case Label::ab_ab::a1b2_a2b1:	return {Label::ab::a2b1, Label::ab::a1b2};
			default:	throw std::invalid_argument("CS_Matrix_Tools::split_label");
		}
	}

	template<typename Tdata>
	Global_Func::To_Real_t<Tdata> cal_uplimit(
		const Uplimit_Type &uplimit_type,
		const Tensor<Tdata> &D)
	{
		using Tlim = Global_Func::To_Real_t<Tdata>;

		auto three_0 = [&D](const std::function<Tlim(Tensor<Tdata>)> &func) -> Tlim
		{
			std::valarray<Tlim> uplimits(D.shape[0]);
			Tensor<Tdata> D_sub({D.shape[1], D.shape[2]});
			for(std::size_t i0=0; i0<D.shape[0]; ++i0)
			{
				memcpy(
					D_sub.ptr(),
					D.ptr()+i0*D.shape[1]*D.shape[2],
					sizeof(Tdata)*D.shape[1]*D.shape[2]);
				uplimits[i0] = func(D_sub);
			}
			return uplimits.max();
		};

		auto three_1 = [&D](const std::function<Tlim(Tensor<Tdata>)> &func) -> Tlim
		{
			std::valarray<Tlim> uplimits(D.shape[1]);
			Tensor<Tdata> D_sub({D.shape[0], D.shape[2]});
			for(std::size_t i1=0; i1<D.shape[1]; ++i1)
			{
				for(std::size_t i0=0; i0<D.shape[0]; ++i0)
					memcpy(
						D_sub.ptr()+i0*D.shape[2],
						D.ptr()+(i0*D.shape[1]+i1)*D.shape[2],
						sizeof(Tdata)*D.shape[2]);
				uplimits[i1] = func(D_sub);
			}
			return uplimits.max();
		};

		auto three_2 = [&D](const std::function<Tlim(Tensor<Tdata>)> &func) -> Tlim
		{
			std::vector<Tensor<Tdata>> Ds_sub;
			Ds_sub.reserve(D.shape[2]);
			for(std::size_t i2=0; i2<D.shape[2]; ++i2)
				Ds_sub.emplace_back(Shape_Vector{D.shape[0],D.shape[1]});

			const Tdata* D_ptr = D.ptr();
			std::vector<Tdata*> Ds_sub_ptr(D.shape[2]);
			for(std::size_t i2=0; i2<D.shape[2]; ++i2)
				Ds_sub_ptr[i2] = Ds_sub[i2].ptr();

			for(std::size_t i0=0; i0<D.shape[0]; ++i0)
				for(std::size_t i1=0; i1<D.shape[1]; ++i1)
					for(std::size_t i2=0; i2<D.shape[2]; ++i2)
						*(Ds_sub_ptr[i2]++) = *(D_ptr++);

			std::valarray<Tlim> uplimits(D.shape[2]);
			for(std::size_t i2=0; i2<D.shape[2]; ++i2)
				uplimits[i2] = func(Ds_sub[i2]);

			return uplimits.max();
		};

		auto three_2_norm = [&D]() -> Tlim
		{
#ifdef __MKL_RI
			const Tensor<Tdata> Ds_sub = Blas_Interface::omatcopy(
				'T', Tdata{1.0},
				D.reshape({D.shape[0]*D.shape[1],D.shape[2]}));
#else
			Tensor<Tdata> Ds_sub({D.shape[2], D.shape[0]*D.shape[1]});

			std::vector<Tdata*> Ds_sub_ptr(D.shape[2]);
			for(std::size_t i2=0; i2<D.shape[2]; ++i2)
				Ds_sub_ptr[i2] = Ds_sub.ptr()+i2*Ds_sub.shape[1]-1;

			const Tdata* D_ptr = D.ptr()-1;
			const std::size_t size2 = D.shape[2];
			for(std::size_t i01=0; i01<Ds_sub.shape[1]; ++i01)
				for(std::size_t i2=0; i2<size2; ++i2)
					*(++Ds_sub_ptr[i2]) = *(++D_ptr);
#endif

			Tensor_Wrapper<Tdata> D_sub({D.shape[0],D.shape[1]}, nullptr);
			std::valarray<Global_Func::To_Real_t<Tdata>> uplimits(D.shape[2]);
			for(std::size_t i2=0; i2<D.shape[2]; ++i2)
			{
				D_sub.ptr_ = Ds_sub.ptr()+i2*Ds_sub.shape[1];
				uplimits[i2] = D_sub.norm(2);
			}
			return uplimits.max();
		};

		auto norm = [](const Tensor<Tdata> &D) -> Tlim
		{
			return D.norm(2);
		};

		auto square = [](const Tensor<Tdata> &D) -> Tlim
		{
			return std::sqrt(Blas_Interface::gemm('C', 'N', Tdata(1), D, D).norm(2));
		};

		switch(uplimit_type)
		{
			case Uplimit_Type::norm_two:
				return norm(D);
			case Uplimit_Type::norm_three_0:
				return three_0(norm);
			case Uplimit_Type::norm_three_1:
				return three_1(norm);
			case Uplimit_Type::norm_three_2:
//				return three_2(norm);
				return three_2_norm();
			case Uplimit_Type::square_two:
				return square(D);
			case Uplimit_Type::square_three_0:
				return three_0(square);
			case Uplimit_Type::square_three_1:
				return three_1(square);
			case Uplimit_Type::square_three_2:
				return three_2(square);
			default:
				throw std::invalid_argument("CS_Matrix_Tools::cal_uplimit");
		}
	}

	template<typename Tkey, typename Tvalue>
	auto cal_uplimit(
		const Uplimit_Type &uplimit_type,
		const std::map<Tkey,Tvalue> &Ds)
	-> std::map<Tkey, decltype(cal_uplimit(uplimit_type,Ds.begin()->second))>
	{
		std::map<Tkey, decltype(cal_uplimit(uplimit_type,Ds.begin()->second))> uplimits;
		for(const auto &Ds_tmp : Ds)
			uplimits[Ds_tmp.first] = cal_uplimit(uplimit_type, Ds_tmp.second);
		return uplimits;
	}
}

}