// ===================
//  Author: Peize Lin
//  date: 2021.10.31
// ===================

#pragma once

#include "Tensor-test.h"
#include "RI/global/Tensor.h"
#include "unittests/print_stl.h"
#include <iostream>
#include <complex>

namespace Tensor_Test
{

	template<typename T>
	RI::Tensor<T> init_real_1()
	{
		RI::Tensor<T> t({2,2});
		t(0,0)=1;	t(0,1)=2;	t(1,0)=3; 	t(1,1)=4;
		return t;
	}
	template<typename T>
	RI::Tensor<T> init_real_2()
	{
		RI::Tensor<T> t({2,2});
		t(0,0)=5;	t(0,1)=6;	t(1,0)=7; 	t(1,1)=8;
		return t;
	}
	template<typename T>
	RI::Tensor<std::complex<T>> init_complex_1()
	{
		RI::Tensor<std::complex<T>> t({2,2});
		t(0,0)={1,2};	t(0,1)={3,4};	t(1,0)={5,6}; 	t(1,1)={7,8};
		return t;
	}
	template<typename T>
	RI::Tensor<std::complex<T>> init_complex_2()
	{
		RI::Tensor<std::complex<T>> t({2,2});
		t(0,0)={9,10};	t(0,1)={11,12};	t(1,0)={13,14}; 	t(1,1)={15,16};
		return t;
	}

	template<typename T>
	void test_multiply_real(const char transA, const char transB)
	{
		const RI::Tensor<T> t1=init_real_1<T>(), t2=init_real_2<T>();
		//RI::Tensor<T> t({2,2});
		//RI::Blas_Interface::gemm(transA,transB,2,2,2,1,t1.ptr(),t2.ptr(),0,t.ptr());
		//std::cout<<t<<std::endl;
		std::cout<<RI::Blas_Interface::gemm(transA, transB, T{1}, t1, t2)<<std::endl;
	}

	template<typename T>
	void test_multiply_complex(const char transA, const char transB)
	{
		const RI::Tensor<std::complex<T>> t1=init_complex_1<T>(), t2=init_complex_2<T>();
		//RI::Tensor<std::complex<T>> t({2,2});
		//RI::Blas_Interface::gemm(transA,transB,2,2,2,1,t1.ptr(),t2.ptr(),0,t.ptr());
		//std::cout<<t<<std::endl;
		std::cout<<RI::Blas_Interface::gemm(transA, transB, std::complex<T>{1}, t1, t2)<<std::endl;
	}

	static void test_multiply_2()
	{
		test_multiply_real<float>('N','N');		// [[19,22],[43,50]]
		test_multiply_real<float>('N','T');		// [[17,23],[39,53]]
		test_multiply_real<float>('N','C');		// [[17,23],[39,53]]

		test_multiply_real<double>('N','N');		// [[19,22],[43,50]]
		test_multiply_real<double>('N','T');		// [[17,23],[39,53]]
		test_multiply_real<double>('N','C');		// [[17,23],[39,53]]

		test_multiply_complex<float>('N','N');		// [[-28+122i, -32+142i], [-36+306i, 40+358i]]
		test_multiply_complex<float>('N','T');		// [[-26+108i, -34+148i], [-34+276i, -42+380i]]
		test_multiply_complex<float>('N','C');		// [[110+16i, 150+24i], [278+8i, 382+16i]]

		test_multiply_complex<double>('N','N');		// [[-28+122i, -32+142i], [-36+306i, 40+358i]]
		test_multiply_complex<double>('N','T');		// [[-26+108i, -34+148i], [-34+276i, -42+380i]]
		test_multiply_complex<double>('N','C');		// [[110+16i, 150+24i], [278+8i, 382+16i]]
	}
}