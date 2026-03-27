// ===================
//  Author: Peize Lin
//  date: 2022.12.25
// ===================

#pragma once

#include "Lapack_Interface.h"
#include "Global_Func-2.h"

namespace RI
{

namespace Lapack_Interface
{
	// potrf computes the Cholesky factorization of a real symmetric positive definite matrix
	template<typename T>
	inline int potrf( const char &uplo, const int &n, T*const A )
	{
		return potrf(uplo, n, A, n);
	}

	// potri takes potrf's output to perform matrix inversion
	template<typename T>
	inline int potri( const char &uplo, const int &n, T*const A )
	{
		return potri(uplo, n, A, n);
	}

	// solve the eigenproblem Ax=ex, where A is Hermitian
	template<typename T>
	inline int heev(const char &jobz, const char &uplo,
		const int &n, T*const A, Global_Func::To_Real_t<T>*const W)
	{
		return heev(jobz, uplo, n, A, n, W);
	}
}

}