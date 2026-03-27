// ===================
//  Author: Peize Lin
//  date: 2023.08.03
// ===================

#pragma once

#include "Tensor_Multiply-test.hpp"

#define FOR_0(i,N)					for(std::size_t i=0; i<N; ++i)

namespace Tensor_Multiply_Test
{
	// Txy(x0,y0,y1) = Tx(x0,a) * Ty(y0,y1,a)
	template<typename Tdata>
	void x0y0y1_x0a_y0y1a_test()
	{
		std::cout<<"Txy(x0,y0,y1) = Tx(x0,a) * Ty(y0,y1,a)"<<std::endl;
		const std::size_t X0=2, Y0=3, Y1=4, A=5;
		const RI::Tensor<Tdata> Tx = init_tensor<Tdata>({X0,A});
		const RI::Tensor<Tdata> Ty = init_tensor<Tdata>({Y0,Y1,A});
		const RI::Tensor<Tdata> Txy = RI::Tensor_Multiply::x0y0y1_x0a_y0y1a(Tx,Ty);
		FOR_0(x0,X0)
			FOR_0(y0,Y0)
				FOR_0(y1,Y1)
				{
					Tdata result = 0;
					FOR_0(a,A)
						result += Tx(x0,a) * Ty(y0,y1,a);
					assert(Txy(x0,y0,y1) == result);
				}
	}

	// Txy(x0,y1,y2) = Tx(x0,a) * Ty(a,y1,y2)
	template<typename Tdata>
	void x0y1y2_x0a_ay1y2_test()
	{
		std::cout<<"Txy(x0,y1,y2) = Tx(x0,a) * Ty(a,y1,y2)"<<std::endl;
		const std::size_t X0=2, Y1=3, Y2=4, A=5;
		const RI::Tensor<Tdata> Tx = init_tensor<Tdata>({X0,A});
		const RI::Tensor<Tdata> Ty = init_tensor<Tdata>({A,Y1,Y2});
		const RI::Tensor<Tdata> Txy = RI::Tensor_Multiply::x0y1y2_x0a_ay1y2(Tx,Ty);
		FOR_0(x0,X0)
			FOR_0(y1,Y1)
				FOR_0(y2,Y2)
				{
					Tdata result = 0;
					FOR_0(a,A)
						result += Tx(x0,a) * Ty(a,y1,y2);
					assert(Txy(x0,y1,y2) == result);
				}
	}

	// Txy(x1,y0,y1) = Tx(a,x1) * Ty(y0,y1,a)
	template<typename Tdata>
	void x1y0y1_ax1_y0y1a_test()
	{
		std::cout<<"Txy(x1,y0,y1) = Tx(a,x1) * Ty(y0,y1,a)"<<std::endl;
		const std::size_t X1=2, Y0=3, Y1=4, A=5;
		const RI::Tensor<Tdata> Tx = init_tensor<Tdata>({A,X1});
		const RI::Tensor<Tdata> Ty = init_tensor<Tdata>({Y0,Y1,A});
		const RI::Tensor<Tdata> Txy = RI::Tensor_Multiply::x1y0y1_ax1_y0y1a(Tx,Ty);
		FOR_0(x1,X1)
			FOR_0(y0,Y0)
				FOR_0(y1,Y1)
				{
					Tdata result = 0;
					FOR_0(a,A)
						result += Tx(a,x1) * Ty(y0,y1,a);
					assert(Txy(x1,y0,y1) == result);
				}
	}

	// Txy(x1,y1,y2) = Tx(a,x1) * Ty(a,y1,y2)
	template<typename Tdata>
	void x1y1y2_ax1_ay1y2_test()
	{
		std::cout<<"Txy(x1,y1,y2) = Tx(a,x1) * Ty(a,y1,y2)"<<std::endl;
		const std::size_t X1=2, Y1=3, Y2=4, A=5;
		const RI::Tensor<Tdata> Tx = init_tensor<Tdata>({A,X1});
		const RI::Tensor<Tdata> Ty = init_tensor<Tdata>({A,Y1,Y2});
		const RI::Tensor<Tdata> Txy = RI::Tensor_Multiply::x1y1y2_ax1_ay1y2(Tx,Ty);
		FOR_0(x1,X1)
			FOR_0(y1,Y1)
				FOR_0(y2,Y2)
				{
					Tdata result = 0;
					FOR_0(a,A)
						result += Tx(a,x1) * Ty(a,y1,y2);
					assert(Txy(x1,y1,y2) == result);
				}
	}

}

#undef FOR_0
