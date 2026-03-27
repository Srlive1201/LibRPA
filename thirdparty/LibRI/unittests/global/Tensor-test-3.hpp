// ===================
//  Author: Peize Lin
//  date: 2021.10.31
// ===================

#pragma once

#include "Tensor-test.h"
#include "unittests/print_stl.h"
#include <iostream>
#include <complex>

namespace Tensor_Test
{
	template<typename Tdata>
	RI::Tensor<Tdata> init_vector(const std::size_t n)
	{
		RI::Tensor<Tdata> t({n});
		for(std::size_t i=0; i<t.shape[0]; ++i)
			t(i) = i;
		return t;
	}

	template<typename Tdata>
	RI::Tensor<Tdata> init_matrix(const std::size_t nr, const std::size_t nc)
	{
		RI::Tensor<Tdata> t({nr,nc});
		for(std::size_t i0=0; i0<t.shape[0]; ++i0)
			for(std::size_t i1=0; i1<t.shape[1]; ++i1)
				t(i0,i1) = i0*10+i1;
		return t;
	}

	template<typename Tdata>
	void test_property()
	{
		const RI::Tensor<Tdata> t = init_matrix<Tdata>(3,2);

		std::cout<<t.shape<<std::endl;
		std::cout<<t.shape.get_shape_all()<<std::endl;

		std::cout<<t<<std::endl;
		std::cout<<t.reshape({2,3})<<std::endl;
		std::cout<<t.transpose()<<std::endl;
		std::cout<<t.reshape({6})<<std::endl;
	}

	template<typename Tdata>
	void test_add()
	{
		std::cout<<init_matrix<Tdata>(2,3)+init_matrix<Tdata>(2,3)<<std::endl;
		// 0  2  4
		// 20 22 24
	}

	template<typename Tdata>
	void test_multiply()
	{
		const RI::Tensor<Tdata> v3 = init_vector<Tdata>(3);
		const RI::Tensor<Tdata> m23 = init_matrix<Tdata>(2,3);
		const RI::Tensor<Tdata> m34 = init_matrix<Tdata>(3,4);

		std::cout<<v3<<std::endl;
		std::cout<<v3*v3<<std::endl<<std::endl;
		// 5

		std::cout<<m23<<std::endl;
		std::cout<<v3<<std::endl;
		std::cout<<m23*v3<<std::endl<<std::endl;
		// 5 35

		std::cout<<v3<<std::endl;
		std::cout<<m34<<std::endl;
		std::cout<<v3*m34<<std::endl<<std::endl;
		// 50 53 56 59

		std::cout<<m23<<std::endl;
		std::cout<<m34<<std::endl;
		std::cout<<m23*m34<<std::endl<<std::endl;
		// 50  53  56  59
		// 350 383 416 449
	}

	template<typename Tdata>
	void test_norm()
	{
		RI::Tensor<double> td({2});
		td(0)=3;	td(1)=-4;
		std::cout<<td.norm(1)<<std::endl;									// 7
		std::cout<<td.norm(2)<<std::endl;									// 5
		std::cout<<td.norm(3)<<std::endl;									// 4.497941445275415
		std::cout<<td.norm(std::numeric_limits<double>::max())<<std::endl;	// 4
	}

	template<typename Tdata>
	void test_norm_complex()
	{
		RI::Tensor<std::complex<Tdata>> tc({2});
		tc(0)={-3,4}; tc(1)={6,-8};
		std::cout<<tc.norm(1)<<std::endl;									// 15
		std::cout<<tc.norm(2)<<std::endl;									// 11.180339887498949
		std::cout<<tc.norm(3)<<std::endl;									// 10.400419115259519
		std::cout<<tc.norm(std::numeric_limits<double>::max())<<std::endl;	// 10
	}

	static void test_operator_all_3()
	{
		test_property<float>();
		test_property<double>();
		test_property<std::complex<float>>();
		test_property<std::complex<double>>();
		test_add<float>();
		test_add<double>();
		test_add<std::complex<float>>();
		test_add<std::complex<double>>();
		test_multiply<float>();
		test_multiply<double>();
		test_multiply<std::complex<float>>();
		test_multiply<std::complex<double>>();
		test_norm<float>();
		test_norm<double>();
		test_norm<std::complex<float>>();
		test_norm<std::complex<double>>();
		test_norm_complex<float>();
		test_norm_complex<double>();
	}
}
