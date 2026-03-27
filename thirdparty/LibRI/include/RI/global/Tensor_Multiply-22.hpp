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
	// Txy(x0,y0) = Tx(x0,a) * Ty(y0,a)
	template<typename Tdata>
	Tensor<Tdata> x0y0_x0a_y0a(const Tensor<Tdata> &Tx, const Tensor<Tdata> &Ty)
	{
		assert(Tx.shape.size()==2);
		assert(Ty.shape.size()==2);
		Tensor<Tdata> Txy({Tx.shape[0], Ty.shape[0]});
		Blas_Interface::gemm(
			'N', 'T',
			Tx.shape[0],
			Ty.shape[0],
			Tx.shape[1],
			Tdata(1.0), Tx.ptr(), Ty.ptr(),
			Tdata(0.0), Txy.ptr());
		return Txy;
	}

	// Txy(x0,y1) = Tx(x0,a) * Ty(a,y1)
	template<typename Tdata>
	Tensor<Tdata> x0y1_x0a_ay1(const Tensor<Tdata> &Tx, const Tensor<Tdata> &Ty)
	{
		assert(Tx.shape.size()==2);
		assert(Ty.shape.size()==2);
		Tensor<Tdata> Txy({Tx.shape[0], Ty.shape[1]});
		Blas_Interface::gemm(
			'N', 'N',
			Tx.shape[0],
			Ty.shape[1],
			Tx.shape[1],
			Tdata(1.0), Tx.ptr(), Ty.ptr(),
			Tdata(0.0), Txy.ptr());
		return Txy;
	}

	// Txy(x1,y0) = Tx(a,x1) * Ty(y0,a)
	template<typename Tdata>
	Tensor<Tdata> x1y0_ax1_y0a(const Tensor<Tdata> &Tx, const Tensor<Tdata> &Ty)
	{
		assert(Tx.shape.size()==2);
		assert(Ty.shape.size()==2);
		Tensor<Tdata> Txy({Tx.shape[1], Ty.shape[0]});
		Blas_Interface::gemm(
			'T', 'T',
			Tx.shape[1],
			Ty.shape[0],
			Tx.shape[0],
			Tdata(1.0), Tx.ptr(), Ty.ptr(),
			Tdata(0.0), Txy.ptr());
		return Txy;
	}

	// Txy(x1,y1) = Tx(a,x1) * Ty(a,y1)
	template<typename Tdata>
	Tensor<Tdata> x1y1_ax1_ay1(const Tensor<Tdata> &Tx, const Tensor<Tdata> &Ty)
	{
		assert(Tx.shape.size()==2);
		assert(Ty.shape.size()==2);
		Tensor<Tdata> Txy({Tx.shape[1], Ty.shape[1]});
		Blas_Interface::gemm(
			'T', 'N',
			Tx.shape[1],
			Ty.shape[1],
			Tx.shape[0],
			Tdata(1.0), Tx.ptr(), Ty.ptr(),
			Tdata(0.0), Txy.ptr());
		return Txy;
	}

}

}
