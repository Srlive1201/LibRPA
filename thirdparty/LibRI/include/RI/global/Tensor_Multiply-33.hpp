// ===================
//  Author: Peize Lin
//  date: 2023.08.03
// ===================

#pragma once

#include "Tensor_Multiply.h"

namespace RI
{

namespace Tensor_Multiply
{
	// Txy(x0,y0) = Tx(x0,a,b) * Ty(y0,a,b)
	template<typename Tdata>
	Tensor<Tdata> x0y0_x0ab_y0ab(const Tensor<Tdata> &Tx, const Tensor<Tdata> &Ty)
	{
		assert(Tx.shape.size()==3);
		assert(Ty.shape.size()==3);
		Tensor<Tdata> Txy({Tx.shape[0], Ty.shape[0]});
		Blas_Interface::gemm(
			'N', 'T',
			Tx.shape[0],
			Ty.shape[0],
			Tx.shape[1] * Tx.shape[2],
			Tdata(1.0), Tx.ptr(), Ty.ptr(),
			Tdata(0.0), Txy.ptr());
		return Txy;
	}

	// Txy(x0,y2) = Tx(x0,a,b) * Ty(a,b,y2)
	template<typename Tdata>
	Tensor<Tdata> x0y2_x0ab_aby2(const Tensor<Tdata> &Tx, const Tensor<Tdata> &Ty)
	{
		assert(Tx.shape.size()==3);
		assert(Ty.shape.size()==3);
		Tensor<Tdata> Txy({Tx.shape[0], Ty.shape[2]});
		Blas_Interface::gemm(
			'N', 'N',
			Tx.shape[0],
			Ty.shape[2],
			Tx.shape[1] * Tx.shape[2],
			Tdata(1.0), Tx.ptr(), Ty.ptr(),
			Tdata(0.0), Txy.ptr());
		return Txy;
	}

	// Txy(x2,y0) = Tx(a,b,x2) * Ty(y0,a,b)
	template<typename Tdata>
	Tensor<Tdata> x2y0_abx2_y0ab(const Tensor<Tdata> &Tx, const Tensor<Tdata> &Ty)
	{
		assert(Tx.shape.size()==3);
		assert(Ty.shape.size()==3);
		Tensor<Tdata> Txy({Tx.shape[2], Ty.shape[0]});
		Blas_Interface::gemm(
			'T', 'T',
			Tx.shape[2],
			Ty.shape[0],
			Tx.shape[0] * Tx.shape[1],
			Tdata(1.0), Tx.ptr(), Ty.ptr(),
			Tdata(0.0), Txy.ptr());
		return Txy;
	}

	// Txy(x2,y2) = Tx(a,b,x2) * Ty(a,b,y2)
	template<typename Tdata>
	Tensor<Tdata> x2y2_abx2_aby2(const Tensor<Tdata> &Tx, const Tensor<Tdata> &Ty)
	{
		assert(Tx.shape.size()==3);
		assert(Ty.shape.size()==3);
		Tensor<Tdata> Txy({Tx.shape[2], Ty.shape[2]});
		Blas_Interface::gemm(
			'T', 'N',
			Tx.shape[2],
			Ty.shape[2],
			Tx.shape[0] * Tx.shape[1],
			Tdata(1.0), Tx.ptr(), Ty.ptr(),
			Tdata(0.0), Txy.ptr());
		return Txy;
	}

	// Txy(x0,x1,y0,y1) = Tx(x0,x1,a) * Ty(y0,y1,a)
	template<typename Tdata>
	Tensor<Tdata> x0x1y0y1_x0x1a_y0y1a(const Tensor<Tdata> &Tx, const Tensor<Tdata> &Ty)
	{
		assert(Tx.shape.size()==3);
		assert(Ty.shape.size()==3);
		Tensor<Tdata> Txy({Tx.shape[0], Tx.shape[1], Ty.shape[0], Ty.shape[1]});
		Blas_Interface::gemm(
			'N', 'T',
			Tx.shape[0] * Tx.shape[1],
			Ty.shape[0] * Ty.shape[1],
			Tx.shape[2],
			Tdata(1.0), Tx.ptr(), Ty.ptr(),
			Tdata(0.0), Txy.ptr());
		return Txy;
	}

	// Txy(x0,x1,y1,y2) = Tx(x0,x1,a) * Ty(a,y1,y2)
	template<typename Tdata>
	Tensor<Tdata> x0x1y1y2_x0x1a_ay1y2(const Tensor<Tdata> &Tx, const Tensor<Tdata> &Ty)
	{
		assert(Tx.shape.size()==3);
		assert(Ty.shape.size()==3);
		Tensor<Tdata> Txy({Tx.shape[0], Tx.shape[1], Ty.shape[1], Ty.shape[2]});
		Blas_Interface::gemm(
			'N', 'N',
			Tx.shape[0] * Tx.shape[1],
			Ty.shape[1] * Ty.shape[2],
			Tx.shape[2],
			Tdata(1.0), Tx.ptr(), Ty.ptr(),
			Tdata(0.0), Txy.ptr());
		return Txy;
	}

	// Txy(x1,x2,y0,y1) = Tx(a,x1,x2) * Ty(y0,y1,a)
	template<typename Tdata>
	Tensor<Tdata> x1x2y0y1_ax1x2_y0y1a(const Tensor<Tdata> &Tx, const Tensor<Tdata> &Ty)
	{
		assert(Tx.shape.size()==3);
		assert(Ty.shape.size()==3);
		Tensor<Tdata> Txy({Tx.shape[1], Tx.shape[2], Ty.shape[0], Ty.shape[1]});
		Blas_Interface::gemm(
			'T', 'T',
			Tx.shape[1] * Tx.shape[2],
			Ty.shape[0] * Ty.shape[1],
			Tx.shape[0],
			Tdata(1.0), Tx.ptr(), Ty.ptr(),
			Tdata(0.0), Txy.ptr());
		return Txy;
	}

	// Txy(x1,x2,y1,y2) = Tx(a,x1,x2) * Ty(a,y1,y2)
	template<typename Tdata>
	Tensor<Tdata> x1x2y1y2_ax1x2_ay1y2(const Tensor<Tdata> &Tx, const Tensor<Tdata> &Ty)
	{
		assert(Tx.shape.size()==3);
		assert(Ty.shape.size()==3);
		Tensor<Tdata> Txy({Tx.shape[1], Tx.shape[2], Ty.shape[1], Ty.shape[2]});
		Blas_Interface::gemm(
			'T', 'N',
			Tx.shape[1] * Tx.shape[2],
			Ty.shape[1] * Ty.shape[2],
			Tx.shape[0],
			Tdata(1.0), Tx.ptr(), Ty.ptr(),
			Tdata(0.0), Txy.ptr());
		return Txy;
	}

}

}
