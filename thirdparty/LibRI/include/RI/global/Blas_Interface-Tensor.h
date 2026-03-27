// ===================
//  Author: Peize Lin
//  date: 2021.10.31
// ===================

#pragma once

#include "Blas_Interface-Contiguous.h"
#include "Tensor.h"
#include <cassert>
#include <string>

namespace RI
{

namespace Blas_Interface
{
	// nrm2 = ||x||_2
	template<typename T, template<typename> class Tvec>
	inline Global_Func::To_Real_t<T> nrm2(const Tvec<T> &X)
	{
		return nrm2(X.shape.get_shape_all(), X.ptr());
	}

	// Vy = alpha * Vx + Vy
	template<typename T>
	inline void axpy(const T &alpha, const Tensor<T> &X, Tensor<T> &Y)
	{
		assert(X.shape.get_shape_all() == Y.shape.get_shape_all());
		axpy(X.shape.get_shape_all(), alpha, X.ptr(), Y.ptr());
	}
	template<typename T>
	inline Tensor<T> axpy(const T &alpha, const Tensor<T> &X)
	{
		Tensor<T> Y(X.shape);
		axpy(alpha, X, Y);
		return Y;
	}

	// d = Vx * Vy
	template<typename T>
	inline T dot(const Tensor<T> &X, const Tensor<T> &Y)
	{
		assert(X.shape.get_shape_all() == Y.shape.get_shape_all());
		return dot(X.shape.get_shape_all(), X.ptr(), Y.ptr());
	}

	// d = Vx * Vy
	template<typename T>
	inline T dotu(const Tensor<T> &X, const Tensor<T> &Y)
	{
		assert(X.shape.get_shape_all() == Y.shape.get_shape_all());
		return dotu(X.shape.get_shape_all(), X.ptr(), Y.ptr());
	}

	// d = Vx * Vy
	template<typename T>
	inline T dotc(const Tensor<T> &X, const Tensor<T> &Y)
	{
		assert(X.shape.get_shape_all() == Y.shape.get_shape_all());
		return dotc(X.shape.get_shape_all(), X.ptr(), Y.ptr());
	}

	// Vy = alpha * Ma.? * Vx + beta * Vy
	template<typename T>
	inline void gemv(const char transA,
		const T alpha, const Tensor<T> &A, const Tensor<T> &X,
		const T beta, Tensor<T> &Y)
	{
		assert(A.shape.size()==2);
		assert(X.shape.size()==1);
		assert(Y.shape.size()==1);
		if(transA=='N')
		{
			assert(A.shape[0]==Y.shape[0]);
			assert(A.shape[1]==X.shape[0]);
		}
		else if(transA=='T' || transA=='C')
		{
			assert(A.shape[1]==Y.shape[0]);
			assert(A.shape[0]==X.shape[0]);
		}
		else
		{
			throw std::invalid_argument("trans must be 'N', 'T' or 'C'.\nIn "+std::string(__FILE__)+" line "+std::to_string(__LINE__)+". ");
		}

		gemv(transA, A.shape[0], A.shape[1],
			alpha, A.ptr(), X.ptr(),
			beta, Y.ptr());
	}
	// Vy = alpha * Ma.? * Vx
	template<typename T>
	inline Tensor<T> gemv(const char transA,
		const T alpha, const Tensor<T> &A, const Tensor<T> &X)
	{
		constexpr T beta = 0.0;
		Tensor<T> Y({(transA=='N') ? A.shape[0] : A.shape[1]});
		gemv(transA, alpha, A, X, beta, Y);
		return Y;
	}

	// Mc = alpha * Ma.? * Mb.? + beta * Mc
	template<typename T>
	inline void gemm(const char transA, const char transB,
		const T alpha, const Tensor<T> &A, const Tensor<T> &B,
		const T beta, Tensor<T> &C)
	{
		assert(A.shape.size()==2);
		assert(B.shape.size()==2);
		assert(C.shape.size()==2);
		if(transA=='N')
			assert(A.shape[0]==C.shape[0]);
		else if(transA=='T' || transA=='C')
			assert(A.shape[1]==C.shape[0]);
		else
			throw std::invalid_argument("trans must be 'N', 'T' or 'C'.\nIn "+std::string(__FILE__)+" line "+std::to_string(__LINE__)+". ");
		if(transB=='N')
			assert(B.shape[1]==C.shape[1]);
		else if(transB=='T' || transB=='C')
			assert(B.shape[0]==C.shape[1]);
		else
			throw std::invalid_argument("trans must be 'N', 'T' or 'C'.\nIn "+std::string(__FILE__)+" line "+std::to_string(__LINE__)+". ");

		const std::size_t m = C.shape[0];
		const std::size_t n = C.shape[1];
		const std::size_t k = (transA=='N') ? A.shape[1] : A.shape[0];

		if(transB=='N')
			assert(k==B.shape[0]);
		else
			assert(k==B.shape[1]);

		gemm(transA, transB, m, n, k,
			alpha, A.ptr(), B.ptr(),
			beta, C.ptr());
	}
	// Mc = alpha * Ma.? * Mb.?
	template<typename T>
	inline Tensor<T> gemm(const char transA, const char transB,
		const T alpha, const Tensor<T> &A, const Tensor<T> &B)
	{
		constexpr T beta = 0.0;
		Tensor<T> C({
			(transA=='N') ? A.shape[0] : A.shape[1],
			(transB=='N') ? B.shape[1] : B.shape[0] });
		gemm(transA, transB,
			alpha, A, B,
			beta, C);
		return C;
	}

	// Mc = alpha * Ma   * Ma.T + beta * C
	// Mc = alpha * Ma.T * Ma   + beta * C
	template<typename T>
	inline void syrk(const char uploC, const char transA,
		const double alpha, const Tensor<T> &A,
		const double beta, Tensor<T> &C)
	{
		assert(A.shape.size()==2);
		assert(C.shape.size()==2);
		assert(C.shape[0]==C.shape[1]);
		if(transA=='N')
			assert(A.shape[0]==C.shape[0]);
		else if(transA=='T' || transA=='C')
			assert(A.shape[1]==C.shape[0]);
		else
			throw std::invalid_argument("trans must be 'N', 'T' or 'C'.\nIn "+std::string(__FILE__)+" line "+std::to_string(__LINE__)+". ");

		const std::size_t n = C.shape[0];
		const std::size_t k = (transA=='N') ? A.shape[1] : A.shape[0];
		syrk(uploC, transA, n, k,
			alpha, A.ptr(),
			beta, C.ptr());
	}
	// Mc = alpha * Ma   * Ma.T
	// Mc = alpha * Ma.T * Ma
	template<typename T>
	inline Tensor<T> syrk(const char uploC, const char transA,
		const double alpha, const Tensor<T> &A)
	{
		constexpr double beta = 0.0;
		const std::size_t n = (transA=='N') ? A.shape[0] : A.shape[1];
		Tensor<T> C({n,n});
		syrk(uploC, transA,	alpha, A, beta, C);
		return C;
	}
}


#ifdef __MKL_RI

namespace Blas_Interface
{
	template<typename T>
	inline void imatcopy (const char trans, const T alpha, Tensor<T> &AB)
	{
		assert(AB.shape.size()==2);
		imatcopy ('R', trans, AB.shape[0], AB.shape[1], alpha, AB.ptr());
		switch(std::toupper(trans))
		{
			case 'N':	case 'R':	break;
			case 'T':	case 'C':	AB=AB.reshape({AB.shape[1], AB.shape[0]});	break;
			default:	throw std::invalid_argument("trans cannot be "+std::to_string(trans)+". "+std::string(__FILE__)+" line "+std::to_string(__LINE__));
		}
	}
	template<typename T>
	inline Tensor<T> omatcopy (char trans, const T alpha, const Tensor<T> &A)
	{
		assert(A.shape.size()==2);
		Tensor<T> B;
		switch(std::toupper(trans))
		{
			case 'N':	case 'R':	B = Tensor<T>({A.shape[0], A.shape[1]});	break;
			case 'T':	case 'C':	B = Tensor<T>({A.shape[1], A.shape[0]});	break;
			default:	throw std::invalid_argument("trans cannot be "+std::to_string(trans)+". "+std::string(__FILE__)+" line "+std::to_string(__LINE__));
		}
		omatcopy ('R', trans, A.shape[0], A.shape[1], alpha, A.ptr(), B.ptr());
		return B;
	}
}

#endif

}

#include "Tensor.hpp"
#include "Tensor-multiply.hpp"