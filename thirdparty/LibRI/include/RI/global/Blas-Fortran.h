// ===================
//  Author: Peize Lin
//  date: 2021.10.31
// ===================

#pragma once

#include <complex>

namespace RI
{

extern "C"
{
	// nrm2 = ||x||_2
	float snrm2_(const int*const n, const float*const X, const int*const incX);
	double dnrm2_(const int*const n, const double*const X, const int*const incX);
	float scnrm2_(const int*const n, const std::complex<float>*const X, const int*const incX);
	double dznrm2_(const int*const n, const std::complex<double>*const X, const int*const incX);

	// Vy = alpha * Vx + Vy
	void saxpy_(const int*const n, const float*const alpha, const float*const X, const int*const incX, float*const Y, const int*const incY);
	void daxpy_(const int*const n, const double*const alpha, const double*const X, const int*const incX, double*const Y, const int*const incY);
	void caxpy_(const int*const n, const std::complex<float>*const alpha, const std::complex<float>*const X, const int*const incX, std::complex<float>*const Y, const int*const incY);
	void zaxpy_(const int*const n, const std::complex<double>*const alpha, const std::complex<double>*const X, const int*const incX, std::complex<double>*const Y, const int*const incY);

	// d = Vx * Vy
	float sdot_(const int*const n, const float*const X, const int*const incX, const float*const Y, const int*const incY);
	double ddot_(const int*const n, const double*const X, const int*const incX, const double*const Y, const int*const incY);

	// d = Vx * Vy
	//	reason for passing results as argument instead of returning it:
	//	https://www.numbercrunch.de/blog/2014/07/lost-in-translation/
	// void cdotu_(std::complex<float>*const result, const int*const n, const std::complex<float>*const X, const int*const incX, const std::complex<float>*const Y, const int*const incY);
	// void zdotu_(std::complex<double>*const result, const int*const n, const std::complex<double>*const X, const int*const incX, const std::complex<double>*const Y, const int*const incY);

	// d = Vx * Vy
	// void cdotc_(std::complex<float>*const result, const int*const n, const std::complex<float>*const X, const int*const incX, const std::complex<float>*const Y, const int*const incY);
	// void zdotc_(std::complex<double>*const result, const int*const n, const std::complex<double>*const X, const int*const incX, const std::complex<double>*const Y, const int*const incY);

	// Vy = alpha * Ma.? * Vx + beta * Vy
	void sgemv_(const char*const transA, const int*const m, const int*const n,
		const float*const alpha, const float*const A, const int*const ldA, const float*const X, const int*const incX,
		const float*const beta, float*const Y, const int*const incY);
	void dgemv_(const char*const transA, const int*const m, const int*const n,
		const double*const alpha, const double*const A, const int*const ldA, const double*const X, const int*const incX,
		const double*const beta, double*const Y, const int*const incY);
	void cgemv_(const char*const transA, const int*const m, const int*const n,
		const std::complex<float>*const alpha, const std::complex<float>*const A, const int*const ldA, const std::complex<float>*const X, const int*const incX,
		const std::complex<float>*const beta, std::complex<float>*const Y, const int*const incY);
	void zgemv_(const char*const transA, const int*const m, const int*const n,
		const std::complex<double>*const alpha, const std::complex<double>*const A, const int*const ldA, const std::complex<double>*const X, const int*const incX,
		const std::complex<double>*const beta, std::complex<double>*const Y, const int*const incY);

	// Mc = alpha * Ma.? * Mb.? + beta * Mc
	void sgemm_(const char*const transA, const char*const transB, const int*const m, const int*const n, const int*const k,
		const float*const alpha, const float*const A, const int*const ldA, const float*const B, const int*const ldB,
		const float*const beta, float*const C, const int*const ldC);
	void dgemm_(const char*const transA, const char*const transB, const int*const m, const int*const n, const int*const k,
		const double*const alpha, const double*const A, const int*const ldA, const double*const B, const int*const ldB,
		const double*const beta, double*const C, const int*const ldC);
	void cgemm_(const char*const transA, const char*const transB, const int*const m, const int*const n, const int*const k,
		const std::complex<float>*const alpha, const std::complex<float>*const A, const int*const ldA, const std::complex<float>*const B, const int*const ldB,
		const std::complex<float>*const beta, std::complex<float>*const C, const int*const ldC);
	void zgemm_(const char*const transA, const char*const transB, const int*const m, const int*const n, const int*const k,
		const std::complex<double>*const alpha, const std::complex<double>*const A, const int*const ldA, const std::complex<double>*const B, const int*const ldB,
		const std::complex<double>*const beta, std::complex<double>*const C, const int*const ldC);

	// Mc = alpha * Ma   * Ma.T + beta * C
	// Mc = alpha * Ma.T * Ma   + beta * C
	void dsyrk_(const char*const uploC, const char*const transA, const int*const n, const int*const k,
		const double*const alpha, const double*const A, const int*const ldA,
		const double*const beta, double*const C, const int*const ldC);
}

}