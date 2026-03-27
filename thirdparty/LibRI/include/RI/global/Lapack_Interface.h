// ===================
//  Author: Peize Lin
//  date: 2022.12.25
// ===================

#pragma once

#include "Lapack-Fortran.h"

#include <string>
#include <stdexcept>


#ifdef __MKL_RI
#include <mkl_trans.h>
#endif

#define LAPACK_INFO_CHECK(x) if(const int info=(x))	throw std::runtime_error("info="+std::to_string(info)+".\n"+std::string(__FILE__)+" line "+std::to_string(__LINE__));

namespace RI
{

namespace Lapack_Interface
{
	// potrf computes the Cholesky factorization of a real symmetric positive definite matrix
	inline int potrf( const char &uplo, const int &n, float*const A, const int &lda )
	{
		int info;
		const char uplo_changed = Blas_Interface::change_uplo(uplo);
		spotrf_( &uplo_changed, &n, A, &lda, &info );
		return info;
	}
	inline int potrf( const char &uplo, const int &n, double*const A, const int &lda )
	{
		int info;
		const char uplo_changed = Blas_Interface::change_uplo(uplo);
		dpotrf_( &uplo_changed, &n, A, &lda, &info );
		return info;
	}
	inline int potrf( const char &uplo, const int &n, std::complex<float>*const A, const int &lda )
	{
		int info;
		const char uplo_changed = Blas_Interface::change_uplo(uplo);
		cpotrf_( &uplo_changed, &n, A, &lda, &info );
		return info;
	}
	inline int potrf( const char &uplo, const int &n, std::complex<double>*const A, const int &lda )
	{
		int info;
		const char uplo_changed = Blas_Interface::change_uplo(uplo);
		zpotrf_( &uplo_changed, &n, A, &lda, &info );
		return info;
	}

	// potri takes potrf's output to perform matrix inversion
	inline int potri( const char &uplo, const int &n, float*const A, const int &lda )
	{
		int info;
		const char uplo_changed = Blas_Interface::change_uplo(uplo);
		spotri_( &uplo_changed, &n, A, &lda, &info);
		return info;
	}
	inline int potri( const char &uplo, const int &n, double*const A, const int &lda )
	{
		int info;
		const char uplo_changed = Blas_Interface::change_uplo(uplo);
		dpotri_( &uplo_changed, &n, A, &lda, &info);
		return info;
	}
	inline int potri( const char &uplo, const int &n, std::complex<float>*const A, const int &lda )
	{
		int info;
		const char uplo_changed = Blas_Interface::change_uplo(uplo);
		cpotri_( &uplo_changed, &n, A, &lda, &info);
		return info;
	}
	inline int potri( const char &uplo, const int &n, std::complex<double>*const A, const int &lda )
	{
		int info;
		const char uplo_changed = Blas_Interface::change_uplo(uplo);
		zpotri_( &uplo_changed, &n, A, &lda, &info);
		return info;
	}

	// solve the eigenproblem Ax=ex, where A is Symmetric
	inline int syev(const char &jobz, const char &uplo,
		const int &n, float*const A, const int &lda, float*const W,
		float*const WORK, const int &lwork)
	{
		int info;
		const char uplo_changed = Blas_Interface::change_uplo(uplo);
		ssyev_(&jobz, &uplo_changed, &n, A, &lda, W, WORK, &lwork, &info);
		return info;
	}
	inline int syev(const char &jobz, const char &uplo,
		const int &n, double*const A, const int &lda, double*const W,
		double*const WORK, const int &lwork)
	{
		int info;
		const char uplo_changed = Blas_Interface::change_uplo(uplo);
		dsyev_(&jobz, &uplo_changed, &n, A, &lda, W, WORK, &lwork, &info);
		return info;
	}
	// solve the eigenproblem Ax=ex, where A is Hermitian
	inline int heev(const char &jobz, const char &uplo,
		const int &n, std::complex<float>*const A, const int &lda, float*const W,
		std::complex<float>*const WORK, const int &lwork, float*const RWORK)
	{
		int info;
		const char uplo_changed = Blas_Interface::change_uplo(uplo);
		cheev_(&jobz, &uplo_changed, &n, A, &lda, W, WORK, &lwork, RWORK, &info);
		return info;
	}
	inline int heev(const char &jobz, const char &uplo,
		const int &n, std::complex<double>*const A, const int &lda, double*const W,
		std::complex<double>*const WORK, const int &lwork, double*const RWORK)
	{
		int info;
		const char uplo_changed = Blas_Interface::change_uplo(uplo);
		zheev_(&jobz, &uplo_changed, &n, A, &lda, W, WORK, &lwork, RWORK, &info);
		return info;
	}

	// solve the eigenproblem Ax=ex, where A is Hermitian
	template<typename T,
		typename std::enable_if< std::is_arithmetic<T>::value,int>::type =0>
	inline int heev(const char &jobz, const char &uplo,
		const int &n, T*const A, const int &lda, T*const W)
	{
		T work_tmp=100;
		constexpr int minus_one = -1;
		LAPACK_INFO_CHECK(syev(jobz, uplo, n, A, lda, W, &work_tmp, minus_one));		// get best lwork

		const int lwork = work_tmp;
		std::vector<T> WORK(std::max(1,lwork));
		return syev(jobz, uplo, n, A, lda, W, WORK.data(), lwork);
	}
	template<typename T,
		typename std::enable_if< std::is_arithmetic<T>::value,int>::type =0>
	inline int heev(const char &jobz, const char &uplo,
		const int &n, std::complex<T>*const A, const int &lda, T*const W)
	{
		std::vector<T> RWORK(std::max(1,3*n-2));

		std::complex<T> work_tmp;
		constexpr int minus_one = -1;
		LAPACK_INFO_CHECK(heev(jobz, uplo, n, A, lda, W, &work_tmp, minus_one, RWORK.data()));		// get best lwork

		const int lwork = std::real(work_tmp);
		std::vector<std::complex<T>> WORK(std::max(1,lwork));
		return heev(jobz, uplo, n, A, lda, W, WORK.data(), lwork, RWORK.data());
	}
}

}

#undef LAPACK_INFO_CHECK