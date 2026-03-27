// ===================
//  Author: Peize Lin
//  date: 2022.12.25
// ===================

#pragma once

#include "Lapack_Interface-Contiguous.h"
#include "Tensor.h"
#include <cassert>
#include <string>

namespace RI
{

namespace Lapack_Interface
{
	// potrf computes the Cholesky factorization of a real symmetric positive definite matrix
	template<typename T>
	inline int potrf( const char &uplo, Tensor<T> &A )
	{
		assert(A.shape.size()==2);
		assert(A.shape[0]==A.shape[1]);
		return potrf(uplo, A.shape[0], A.ptr());
	}

	// potri takes potrf's output to perform matrix inversion
	template<typename T>
	inline int potri( const char &uplo, Tensor<T> &A )
	{
		assert(A.shape.size()==2);
		assert(A.shape[0]==A.shape[1]);
		return potri(uplo, A.shape[0], A.ptr());
	}

	// solve the eigenproblem Ax=ex, where A is Hermitian
	template<typename T>
	inline int heev(const char &jobz, const char &uplo,
		Tensor<T> &A, std::vector<Global_Func::To_Real_t<T>> &W)
	{
		assert(A.shape.size()==2);
		assert(A.shape[0]==A.shape[1]);
		assert(A.shape[0]==W.size());
		return heev(jobz, uplo, A.shape[0], A.ptr(), W.data());
	}
}

}