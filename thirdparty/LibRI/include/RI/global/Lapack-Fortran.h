// ===================
//  Author: Peize Lin
//  date: 2022.12.25
// ===================

#pragma once

#include <complex>

namespace RI
{

extern "C"
{
	// potrf computes the Cholesky factorization of a real symmetric positive definite matrix
	void spotrf_(const char*const uplo, const int*const n, float*const A, const int*const lda, int*const info);
	void dpotrf_(const char*const uplo, const int*const n, double*const A, const int*const lda, int*const info);
	void cpotrf_(const char*const uplo, const int*const n, std::complex<float>*const A, const int*const lda, int*const info);
	void zpotrf_(const char*const uplo, const int*const n, std::complex<double>*const A, const int*const lda, int*const info);

	// potri takes potrf's output to perform matrix inversion
	void spotri_(const char*const uplo, const int*const n, float*const A, const int*const lda, int*const info);
	void dpotri_(const char*const uplo, const int*const n, double*const A, const int*const lda, int*const info);
	void cpotri_(const char*const uplo, const int*const n, std::complex<float>*const A, const int*const lda, int*const info);
	void zpotri_(const char*const uplo, const int*const n, std::complex<double>*const A, const int*const lda, int*const info);

	// solve the eigenproblem Ax=ex, where A is Symmetric
	void ssyev_(const char*const jobz, const char*const uplo,
		const int*const n, float*const A, const int*const lda, float*const W,
		float*const WORK, const int*const lwork, int*const info);
	void dsyev_(const char*const jobz, const char*const uplo,
		const int*const n, double*const A, const int*const lda, double*const W,
		double*const WORK, const int*const lwork, int*const info);
	// solve the eigenproblem Ax=ex, where A is Hermitian
	void cheev_(const char*const jobz, const char*const uplo,
		const int*const n, std::complex<float>*const A, const int*const lda, float*const W,
		std::complex<float>*const WORK, const int*const lwork, float*const RWORK, int*const info);
	void zheev_(const char*const jobz, const char*const uplo,
		const int*const n, std::complex<double>*const A, const int*const lda, double*const W,
		std::complex<double>*const WORK, const int*const lwork, double*const RWORK, int*const info);
}

}