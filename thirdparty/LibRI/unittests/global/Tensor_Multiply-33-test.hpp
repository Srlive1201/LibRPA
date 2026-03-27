// ===================
//  Author: Peize Lin
//  date: 2023.08.03
// ===================

#pragma once

#include "Tensor_Multiply-test.hpp"

#define FOR_0(i,N)					for(std::size_t i=0; i<N; ++i)

namespace Tensor_Multiply_Test
{
	// Txy(x0,y0) = Tx(x0,a,b) * Ty(y0,a,b)
	template<typename Tdata>
	void x0y0_x0ab_y0ab_test()
	{
		std::cout<<"Txy(x0,y0) = Tx(x0,a,b) * Ty(y0,a,b)"<<std::endl;
		const std::size_t X0=2, Y0=3, A=4, B=5;
		const RI::Tensor<Tdata> Tx = init_tensor<Tdata>({X0,A,B});
		const RI::Tensor<Tdata> Ty = init_tensor<Tdata>({Y0,A,B});
		const RI::Tensor<Tdata> Txy = RI::Tensor_Multiply::x0y0_x0ab_y0ab(Tx,Ty);
		FOR_0(x0,X0)
			FOR_0(y0,Y0)
			{
				Tdata result = 0;
				FOR_0(a,A)
					FOR_0(b,B)
						result += Tx(x0,a,b) * Ty(y0,a,b);
				assert(Txy(x0,y0) == result);
			}
	}

	// Txy(x0,y2) = Tx(x0,a,b) * Ty(a,b,y2)
	template<typename Tdata>
	void x0y2_x0ab_aby2_test()
	{
		std::cout<<"Txy(x0,y2) = Tx(x0,a,b) * Ty(a,b,y2)"<<std::endl;
		const std::size_t X0=2, Y2=3, A=4, B=5;
		const RI::Tensor<Tdata> Tx = init_tensor<Tdata>({X0,A,B});
		const RI::Tensor<Tdata> Ty = init_tensor<Tdata>({A,B,Y2});
		const RI::Tensor<Tdata> Txy = RI::Tensor_Multiply::x0y2_x0ab_aby2(Tx,Ty);
		FOR_0(x0,X0)
			FOR_0(y2,Y2)
			{
				Tdata result = 0;
				FOR_0(a,A)
					FOR_0(b,B)
						result += Tx(x0,a,b) * Ty(a,b,y2);
				assert(Txy(x0,y2) == result);
			}
	}

	// Txy(x2,y0) = Tx(a,b,x2) * Ty(y0,a,b)
	template<typename Tdata>
	void x2y0_abx2_y0ab_test()
	{
		std::cout<<"Txy(x2,y0) = Tx(a,b,x2) * Ty(y0,a,b)"<<std::endl;
		const std::size_t X2=2, Y0=3, A=4, B=5;
		const RI::Tensor<Tdata> Tx = init_tensor<Tdata>({A,B,X2});
		const RI::Tensor<Tdata> Ty = init_tensor<Tdata>({Y0,A,B});
		const RI::Tensor<Tdata> Txy = RI::Tensor_Multiply::x2y0_abx2_y0ab(Tx,Ty);
		FOR_0(x2,X2)
			FOR_0(y0,Y0)
			{
				Tdata result = 0;
				FOR_0(a,A)
					FOR_0(b,B)
						result += Tx(a,b,x2) * Ty(y0,a,b);
				assert(Txy(x2,y0) == result);
			}
	}

	// Txy(x2,y2) = Tx(a,b,x2) * Ty(a,b,y2)
	template<typename Tdata>
	void x2y2_abx2_aby2_test()
	{
		std::cout<<"Txy(x2,y2) = Tx(a,b,x2) * Ty(a,b,y2)"<<std::endl;
		const std::size_t X2=2, Y2=3, A=4, B=5;
		const RI::Tensor<Tdata> Tx = init_tensor<Tdata>({A,B,X2});
		const RI::Tensor<Tdata> Ty = init_tensor<Tdata>({A,B,Y2});
		const RI::Tensor<Tdata> Txy = RI::Tensor_Multiply::x2y2_abx2_aby2(Tx,Ty);
		FOR_0(x2,X2)
			FOR_0(y2,Y2)
			{
				Tdata result = 0;
				FOR_0(a,A)
					FOR_0(b,B)
						result += Tx(a,b,x2) * Ty(a,b,y2);
				assert(Txy(x2,y2) == result);
			}
	}

	// Txy(x0,x1,y0,y1) = Tx(x0,x1,a) * Ty(y0,y1,a)
	template<typename Tdata>
	void x0x1y0y1_x0x1a_y0y1a_test()
	{
		std::cout<<"Txy(x0,x1,y0,y1) = Tx(x0,x1,a) * Ty(y0,y1,a)"<<std::endl;
		const std::size_t X0=2, X1=3, Y0=4, Y1=5, A=6;
		const RI::Tensor<Tdata> Tx = init_tensor<Tdata>({X0,X1,A});
		const RI::Tensor<Tdata> Ty = init_tensor<Tdata>({Y0,Y1,A});
		const RI::Tensor<Tdata> Txy = RI::Tensor_Multiply::x0x1y0y1_x0x1a_y0y1a(Tx,Ty);
		FOR_0(x0,X0)
			FOR_0(x1,X1)
				FOR_0(y0,Y0)
					FOR_0(y1,Y1)
					{
						Tdata result = 0;
						FOR_0(a,A)
							result += Tx(x0,x1,a) * Ty(y0,y1,a);
						assert(Txy(x0,x1,y0,y1) == result);
					}
	}

	// Txy(x0,x1,y1,y2) = Tx(x0,x1,a) * Ty(a,y1,y2)
	template<typename Tdata>
	void x0x1y1y2_x0x1a_ay1y2_test()
	{
		std::cout<<"Txy(x0,x1,y1,y2) = Tx(x0,x1,a) * Ty(a,y1,y2)"<<std::endl;
		const std::size_t X0=2, X1=3, Y1=4, Y2=5, A=6;
		const RI::Tensor<Tdata> Tx = init_tensor<Tdata>({X0,X1,A});
		const RI::Tensor<Tdata> Ty = init_tensor<Tdata>({A,Y1,Y2});
		const RI::Tensor<Tdata> Txy = RI::Tensor_Multiply::x0x1y1y2_x0x1a_ay1y2(Tx,Ty);
		FOR_0(x0,X0)
			FOR_0(x1,X1)
				FOR_0(y1,Y1)
					FOR_0(y2,Y2)
					{
						Tdata result = 0;
						FOR_0(a,A)
							result += Tx(x0,x1,a) * Ty(a,y1,y2);
						assert(Txy(x0,x1,y1,y2) == result);
					}
	}

	// Txy(x1,x2,y0,y1) = Tx(a,x1,x2) * Ty(y0,y1,a)
	template<typename Tdata>
	void x1x2y0y1_ax1x2_y0y1a_test()
	{
		std::cout<<"Txy(x1,x2,y0,y1) = Tx(a,x1,x2) * Ty(y0,y1,a)"<<std::endl;
		const std::size_t X1=2, X2=3, Y0=4, Y1=5, A=6;
		const RI::Tensor<Tdata> Tx = init_tensor<Tdata>({A,X1,X2});
		const RI::Tensor<Tdata> Ty = init_tensor<Tdata>({Y0,Y1,A});
		const RI::Tensor<Tdata> Txy = RI::Tensor_Multiply::x1x2y0y1_ax1x2_y0y1a(Tx,Ty);
		FOR_0(x1,X1)
			FOR_0(x2,X2)
				FOR_0(y0,Y0)
					FOR_0(y1,Y1)
					{
						Tdata result = 0;
						FOR_0(a,A)
							result += Tx(a,x1,x2) * Ty(y0,y1,a);
						assert(Txy(x1,x2,y0,y1) == result);
					}
	}

	// Txy(x1,x2,y1,y2) = Tx(a,x1,x2) * Ty(a,y1,y2)
	template<typename Tdata>
	void x1x2y1y2_ax1x2_ay1y2_test()
	{
		std::cout<<"Txy(x1,x2,y1,y2) = Tx(a,x1,x2) * Ty(a,y1,y2)"<<std::endl;
		const std::size_t X1=2, X2=3, Y1=4, Y2=5, A=6;
		const RI::Tensor<Tdata> Tx = init_tensor<Tdata>({A,X1,X2});
		const RI::Tensor<Tdata> Ty = init_tensor<Tdata>({A,Y1,Y2});
		const RI::Tensor<Tdata> Txy = RI::Tensor_Multiply::x1x2y1y2_ax1x2_ay1y2(Tx,Ty);
		FOR_0(x1,X1)
			FOR_0(x2,X2)
				FOR_0(y1,Y1)
					FOR_0(y2,Y2)
					{
						Tdata result = 0;
						FOR_0(a,A)
							result += Tx(a,x1,x2) * Ty(a,y1,y2);
						assert(Txy(x1,x2,y1,y2) == result);
					}
	}

}

#undef FOR_0
