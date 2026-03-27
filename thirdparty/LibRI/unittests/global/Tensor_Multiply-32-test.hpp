// ===================
//  Author: Peize Lin
//  date: 2023.08.03
// ===================

#pragma once

#include "Tensor_Multiply-test.hpp"

#define FOR_0(i,N)					for(std::size_t i=0; i<N; ++i)

namespace Tensor_Multiply_Test
{
	// Txy(x0,x1,y0) = Tx(x0,x1,a) * Ty(y0,a)
	template<typename Tdata>
	void x0x1y0_x0x1a_y0a_test()
	{
		std::cout<<"Txy(x0,x1,y0) = Tx(x0,x1,a) * Ty(y0,a)"<<std::endl;
		const std::size_t X0=2, X1=3, Y0=4, A=5;
		const RI::Tensor<Tdata> Tx = init_tensor<Tdata>({X0,X1,A});
		const RI::Tensor<Tdata> Ty = init_tensor<Tdata>({Y0,A});
		const RI::Tensor<Tdata> Txy = RI::Tensor_Multiply::x0x1y0_x0x1a_y0a(Tx,Ty);
		FOR_0(x0,X0)
			FOR_0(x1,X1)
				FOR_0(y0,Y0)
				{
					Tdata result = 0;
					FOR_0(a,A)
						result += Tx(x0,x1,a) * Ty(y0,a);
					assert(Txy(x0,x1,y0) == result);
				}
	}

	// Txy(x0,x1,y1) = Tx(x0,x1,a) * Ty(a,y1)
	template<typename Tdata>
	void x0x1y1_x0x1a_ay1_test()
	{
		std::cout<<"Txy(x0,x1,y1) = Tx(x0,x1,a) * Ty(a,y1)"<<std::endl;
		const std::size_t X0=2, X1=3, Y1=4, A=5;
		const RI::Tensor<Tdata> Tx = init_tensor<Tdata>({X0,X1,A});
		const RI::Tensor<Tdata> Ty = init_tensor<Tdata>({A,Y1});
		const RI::Tensor<Tdata> Txy = RI::Tensor_Multiply::x0x1y1_x0x1a_ay1(Tx,Ty);
		FOR_0(x0,X0)
			FOR_0(x1,X1)
				FOR_0(y1,Y1)
				{
					Tdata result = 0;
					FOR_0(a,A)
						result += Tx(x0,x1,a) * Ty(a,y1);
					assert(Txy(x0,x1,y1) == result);
				}
	}

	// Txy(x1,x2,y0) = Tx(a,x1,x2) * Ty(y0,a)
	template<typename Tdata>
	void x1x2y0_ax1x2_y0a_test()
	{
		std::cout<<"Txy(x1,x2,y0) = Tx(a,x1,x2) * Ty(y0,a)"<<std::endl;
		const std::size_t X1=2, X2=3, Y0=4, A=5;
		const RI::Tensor<Tdata> Tx = init_tensor<Tdata>({A,X1,X2});
		const RI::Tensor<Tdata> Ty = init_tensor<Tdata>({Y0,A});
		const RI::Tensor<Tdata> Txy = RI::Tensor_Multiply::x1x2y0_ax1x2_y0a(Tx,Ty);
		FOR_0(x1,X1)
			FOR_0(x2,X2)
				FOR_0(y0,Y0)
				{
					Tdata result = 0;
					FOR_0(a,A)
						result += Tx(a,x1,x2) * Ty(y0,a);
					assert(Txy(x1,x2,y0) == result);
				}
	}

	// Txy(x1,x2,y1) = Tx(a,x1,x2) * Ty(a,y1)
	template<typename Tdata>
	void x1x2y1_ax1x2_ay1_test()
	{
		std::cout<<"Txy(x1,x2,y1) = Tx(a,x1,x2) * Ty(a,y1)"<<std::endl;
		const std::size_t X1=2, X2=3, Y1=4, A=5;
		const RI::Tensor<Tdata> Tx = init_tensor<Tdata>({A,X1,X2});
		const RI::Tensor<Tdata> Ty = init_tensor<Tdata>({A,Y1});
		const RI::Tensor<Tdata> Txy = RI::Tensor_Multiply::x1x2y1_ax1x2_ay1(Tx,Ty);
		FOR_0(x1,X1)
			FOR_0(x2,X2)
				FOR_0(y1,Y1)
				{
					Tdata result = 0;
					FOR_0(a,A)
						result += Tx(a,x1,x2) * Ty(a,y1);
					assert(Txy(x1,x2,y1) == result);
				}
	}

}

#undef FOR_0
