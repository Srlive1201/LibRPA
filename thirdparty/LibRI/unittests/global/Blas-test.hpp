// ===================
//  Author: Peize Lin
//  date: 2021.10.31
// ===================

#pragma once

#include <vector>
#include <iostream>

#include "RI/global/Blas_Interface.h"
#include "RI/global/Blas_Interface-Contiguous.h"
#include "RI/global/Blas_Interface-Tensor.h"
#include "RI/global/Tensor.h"
#include "../global/Tensor-test-3.hpp"

#ifdef __MKL_RI
#include <mkl_blas.h>
#endif

namespace Blas_Test
{
	template<typename Tdata>
	void nrm2()
	{
		const RI::Tensor<Tdata> v3 = Tensor_Test::init_vector<Tdata>(3);
		std::cout<<v3<<std::endl;
		std::cout<<RI::Blas_Interface::nrm2(v3)<<std::endl;
		// 2.23607
	}

	template<typename Tdata>
	void axpy()
	{
		const RI::Tensor<Tdata> vx = Tensor_Test::init_vector<Tdata>(3);
		RI::Tensor<Tdata> vy = Tensor_Test::init_vector<Tdata>(3);
		RI::Blas_Interface::axpy(vx.shape[0], 2, vx.ptr(), 1, vy.ptr(), 1);
		std::cout<<vy<<std::endl;
		// 0, 3, 6
		RI::Tensor<Tdata> vz = RI::Blas_Interface::axpy(Tdata(2), vx);
		std::cout<<vz<<std::endl;
		// 0, 2, 4
	}

	template<typename Tdata>
	void dot_1()
	{
		std::vector<Tdata> a = {1,2,3};
		std::vector<Tdata> b = {4,5,6};
		std::cout<<RI::Blas_Interface::dot(3, a.data(), 1, b.data(), 1)<<std::endl;
		// 32
	}

//#ifdef __MKL_RI
//	void dot_2()
//	{
//		const MKL_INT N=3, INC=1;
//		std::vector<double> a = {1,2,3};
//		std::vector<double> b = {4,5,6};
//		std::cout<<DDOT(&N, a.data(), &INC, b.data(), &INC)<<std::endl;
//		// 32
//	}
//#endif

	template<typename Tdata>
	void dot_complex()
	{
		const std::vector<Tdata> a = {{1,2},{3,4}};
		const std::vector<Tdata> b = {{5,6},{7,8}};
//		std::cout<<a<<std::endl<<b<<std::endl;
		std::cout<<RI::Blas_Interface::dotu(a.size(), a.data(), 1, b.data(), 1)<<std::endl;
		/* -18+68i */
		std::cout<<RI::Blas_Interface::dotc(a.size(), a.data(), 1, b.data(), 1)<<std::endl;
		/* 70-8i */

		// test for different incX and incY
		const std::vector<Tdata> a3 = { {1,2}, {1,4}, {4,2}, {3,4}, {2,8}, {1,8}};
		const std::vector<Tdata> b2 = { {5,6}, {5,7}, {7,8}, {7,5}};
		std::cout << RI::Blas_Interface::dotu(a.size(), a3.data(), 3, b2.data(), 2) << std::endl;
		/* -18+68i */
		std::cout << RI::Blas_Interface::dotc(a.size(), a3.data(), 3, b2.data(), 2) << std::endl;
		/* 70-8i */
	}

	template<typename Tdata>
	void gemv()
	{
		const RI::Tensor<Tdata> m23 = Tensor_Test::init_matrix<Tdata>(2,3);
		const RI::Tensor<Tdata> v2 = Tensor_Test::init_vector<Tdata>(2);
		const RI::Tensor<Tdata> v3 = Tensor_Test::init_vector<Tdata>(3);
		RI::Tensor<Tdata> vr2({2});
		RI::Tensor<Tdata> vr3({3});

		std::cout<<m23<<std::endl;
		std::cout<<v3<<std::endl;
		RI::Blas_Interface::gemv('N', m23.shape[0], m23.shape[1], Tdata{1.0}, m23.ptr(), v3.ptr(), Tdata{0.0}, vr2.ptr());
		std::cout<<vr2<<std::endl;
		RI::Blas_Interface::gemv('N', Tdata{1.0}, m23, v3, Tdata{0.0}, vr2);
		std::cout<<vr2<<std::endl;
		std::cout<<RI::Blas_Interface::gemv('N', Tdata{1.0}, m23, v3)<<std::endl;
		// 5 35

		std::cout<<m23.transpose()<<std::endl;
		std::cout<<v2<<std::endl;
		RI::Blas_Interface::gemv('T', m23.shape[0], m23.shape[1], Tdata{1.0}, m23.ptr(), v2.ptr(), Tdata{0.0}, vr3.ptr());
		std::cout<<vr3<<std::endl;
		RI::Blas_Interface::gemv('T', Tdata{1.0}, m23, v2, Tdata{0.0}, vr3);
		std::cout<<vr3<<std::endl;
		std::cout<<RI::Blas_Interface::gemv('T', Tdata{1.0}, m23, v2)<<std::endl;
		// 10 11 12
	}

	template<typename Tdata>
	void gemv_complex()
	{
		const RI::Tensor<Tdata> m23 = Tensor_Test::init_matrix<Tdata>(2,3);
		const RI::Tensor<Tdata> v3 = Tensor_Test::init_vector<Tdata>(3);
		RI::Tensor<Tdata> vr2 = Tensor_Test::init_vector<Tdata>(2);
		std::cout<<m23<<std::endl;
		std::cout<<v3<<std::endl;
		std::cout<<vr2<<std::endl;
		RI::Blas_Interface::gemv('N', Tdata{1.0,2.0}, m23, v3, Tdata{3,-4}, vr2);
		std::cout<<vr2<<std::endl;
		// (5,10)	(38,66)
	}

	template<typename Tdata>
	void gemm_real()
	{
		const RI::Tensor<Tdata> m23 = Tensor_Test::init_matrix<Tdata>(2,3);
		std::cout<<m23<<std::endl;
		std::cout<<RI::Blas_Interface::gemm('T','N',Tdata(1.0),m23,m23)<<std::endl;
		std::cout<<RI::Blas_Interface::gemm('N','T',Tdata(1.0),m23,m23)<<std::endl;
		// 100     110     120
		// 110     122     134
		// 120     134     148

		// 5       35
		// 35      365
	}

	static void syrk()
	{
		const RI::Tensor<double> m = Tensor_Test::init_matrix<double>(2,3);

		std::cout<<m<<std::endl;
		std::cout<<m.transpose()<<std::endl;
		RI::Tensor<double> m_mT = RI::Tensor<double>({2,2});
		RI::Blas_Interface::syrk('U', 'N', m.shape[0], m.shape[1], 1.0, m.ptr(), m.shape[1], 0.0, m_mT.ptr(), m_mT.shape[0]);
		std::cout<<m_mT<<std::endl;
		RI::Blas_Interface::syrk('U', 'N', m.shape[0], m.shape[1], 1.0, m.ptr(), 0.0, m_mT.ptr());
		std::cout<<m_mT<<std::endl;
		RI::Blas_Interface::syrk('U', 'N', 1.0, m, 0.0, m_mT);
		std::cout<<m_mT<<std::endl;
		std::cout<<RI::Blas_Interface::syrk('U', 'N', 1.0, m)<<std::endl;
		// 5 35
		// 0 365

		std::cout<<m.transpose()<<std::endl;
		std::cout<<m<<std::endl;
		RI::Tensor<double> mT_m = RI::Tensor<double>({3,3});
		RI::Blas_Interface::syrk('U', 'T', m.shape[1], m.shape[0], 1.0, m.ptr(), m.shape[1], 0.0, m_mT.ptr(), mT_m.shape[0]);
		std::cout<<mT_m<<std::endl;
		RI::Blas_Interface::syrk('U', 'T', m.shape[1], m.shape[0], 1.0, m.ptr(), 0.0, m_mT.ptr());
		std::cout<<mT_m<<std::endl;
		RI::Blas_Interface::syrk('U', 'T', 1.0, m, 0.0, mT_m);
		std::cout<<mT_m<<std::endl;
		std::cout<<RI::Blas_Interface::syrk('U', 'T', 1.0, m)<<std::endl;
		// 100 110 120
		// 0   122 134
		// 0   0   148
	}

#ifdef __MKL_RI
	template<typename Tdata>
	void matcopy()
	{
		RI::Tensor<Tdata> t1 = Tensor_Test::init_matrix<Tdata>(2,3);
		RI::Blas_Interface::imatcopy('N',Tdata{1},t1);		std::cout<<t1<<std::endl;
		RI::Blas_Interface::imatcopy('T',Tdata{1},t1);		std::cout<<t1<<std::endl;
		RI::Blas_Interface::imatcopy('C',Tdata{1},t1);		std::cout<<t1<<std::endl;
		RI::Blas_Interface::imatcopy('R',Tdata{1},t1);		std::cout<<t1<<std::endl;
		std::cout<<RI::Blas_Interface::omatcopy('N',Tdata{1},t1)<<std::endl;
		std::cout<<RI::Blas_Interface::omatcopy('T',Tdata{1},t1)<<std::endl;
		std::cout<<RI::Blas_Interface::omatcopy('C',Tdata{1},t1)<<std::endl;
		std::cout<<RI::Blas_Interface::omatcopy('R',Tdata{1},t1)<<std::endl;
	}
#endif

	static void test_all()
	{
		nrm2<float>();
		nrm2<double>();
		nrm2<std::complex<float>>();
		nrm2<std::complex<double>>();
		axpy<float>();
		axpy<double>();
		axpy<std::complex<float>>();
		axpy<std::complex<double>>();
		dot_1<float>();
		dot_1<double>();
//		dot_2();
		dot_complex<std::complex<float>>();
		dot_complex<std::complex<double>>();
		gemv<float>();
		gemv<double>();
		gemv<std::complex<float>>();
		gemv<std::complex<double>>();
		gemv_complex<std::complex<float>>();
		gemv_complex<std::complex<double>>();
		gemm_real<float>();
		gemm_real<double>();
		gemm_real<std::complex<float>>();
		gemm_real<std::complex<double>>();
		syrk();
#ifdef __MKL_RI
		matcopy<float>();
		matcopy<double>();
		matcopy<std::complex<float>>();
		matcopy<std::complex<double>>();
#endif
	}
}