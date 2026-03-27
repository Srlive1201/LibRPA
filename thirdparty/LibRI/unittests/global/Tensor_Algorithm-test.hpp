// ===================
//  Author: Peize Lin
//  date: 2022.12.27
// ===================

#pragma once

#include "RI/global/Tensor_Algorithm.h"

#include "../print_stl.h"
#include "Tensor-test.h"

namespace Tensor_Algorithm_Test
{
	template<typename T>
	void test_inverse_matrix_real()
	{
		const RI::Tensor<T> m({2,2});
		m(0,0)=1; m(0,1)=2; m(1,0)=2, m(1,1)=3;
		const RI::Tensor<T> mI = RI::Tensor_Algorithm::inverse_matrix_heev(m);
		std::cout<<m<<std::endl;
		std::cout<<mI<<std::endl;
		std::cout<<mI*m<<std::endl;
		/* m:   1  2
		        2  3
		   mI: -3  2
		        2 -1
		*/
	}
	template<typename T>
	void test_inverse_matrix_complex()
	{
		const RI::Tensor<std::complex<T>> m({2,2});
		m(0,0)={1,0}; m(0,1)={3,4}; m(1,0)={3,-4}, m(1,1)={5,0};
		const RI::Tensor<std::complex<T>> mI = RI::Tensor_Algorithm::inverse_matrix_heev(m);
		std::cout<<m<<std::endl;
		std::cout<<mI<<std::endl;
		std::cout<<mI*m<<std::endl;
		/* m:	(1, 0)			(3,4)
				(3,-4)			(5,0)
		   mI:	(-0.25, 0.00)	( 0.15, 0.20)
				( 0.15,-0.20)	(-0.05, 0.00)
		*/
	}
}