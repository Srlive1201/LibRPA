// ===================
//  Author: Peize Lin
//  date: 2021.10.31
// ===================

#pragma once

#include "Blas_Interface.h"
#include "Global_Func-2.h"

namespace RI
{

namespace Blas_Interface
{
	// nrm2 = ||x||_2
	template<typename T>
	inline Global_Func::To_Real_t<T> nrm2(const int n, const T*const X)
	{
		constexpr int incX = 1;
		return nrm2(n, X, incX);
	}

	// Vy = alpha * Vx + Vy
	template<typename T>
	inline void axpy(const int &n, const T &alpha, const T*const X, T*const Y)
	{
		constexpr int incX = 1;
		constexpr int incY = 1;
		axpy(n, alpha, X, incX, Y, incY);
	}

	// d = Vx * Vy
	template<typename T>
	inline T dot(const int n, const T*const X, const T*const Y)
	{
		constexpr int incX = 1;
		constexpr int incY = 1;
		return dot(n, X, incX, Y, incY);
	}

	// d = Vx * Vy
	template<typename T>
	inline T dotu(const int n, const T*const X, const T*const Y)
	{
		constexpr int incX = 1;
		constexpr int incY = 1;
		return dotu(n, X, incX, Y, incY);
	}

	// d = Vx * Vy
	template<typename T>
	inline T dotc(const int n, const T*const X, const T*const Y)
	{
		constexpr int incX = 1;
		constexpr int incY = 1;
		return dotc(n, X, incX, Y, incY);
	}

	// Vy = alpha * Ma.? * Vx + beta * Vy
	template<typename T>
	inline void gemv(const char transA, const int m, const int n,
		const T alpha, const T*const A, const T*const X,
		const T beta, T*const Y)
	{
		const int ldA = n;
		constexpr int incX = 1;
		constexpr int incY = 1;
		gemv(transA, m, n,
			alpha, A, ldA, X, incX,
			beta, Y, incY);
	}

	// Mc = alpha * Ma.? * Mb.? + beta * Mc
	template<typename T>
	inline void gemm(const char transA, const char transB, const int m, const int n, const int k,
		const T alpha, const T*const A, const T*const B,
		const T beta, T*const C)
	{
		const int ldA = (transA=='N') ? k : m;
		const int ldB = (transB=='N') ? n : k;
		const int ldC = n;
		gemm(transA, transB, m, n, k,
			alpha, A, ldA, B, ldB,
			beta, C, ldC);
	}

	// Mc = alpha * Ma   * Ma.T + beta * C
	// Mc = alpha * Ma.T * Ma   + beta * C
	template<typename T>
	inline void syrk(const char uploC, const char transA, const int n, const int k,
		const double alpha, const T*const A,
		const double beta, T*const C)
	{
		const int ldA = (transA=='N') ? k : n;
		const int ldC = n;
		syrk(uploC, transA, n, k,
			alpha, A, ldA,
			beta, C, ldC);
	}
}



#ifdef __MKL_RI

namespace Blas_Interface
{
	inline size_t get_lda_matcopy(const char ordering, size_t rows, size_t cols)
	{
		switch(std::toupper(ordering))
		{
			case 'R':	return cols;
			case 'C':	return rows;
			default:	throw std::invalid_argument("ordering cannot be "+std::to_string(ordering)+". "+std::string(__FILE__)+" line "+std::to_string(__LINE__));
		}
	}
	inline size_t get_ldb_matcopy(const char ordering, const char trans, size_t rows, size_t cols)
	{
		switch(std::toupper(ordering))
		{
			case 'R':
				switch(std::toupper(trans))
				{
					case 'N':	case 'R':	return cols;
					case 'T':	case 'C':	return rows;
					default:	throw std::invalid_argument("trans cannot be "+std::to_string(trans)+". "+std::string(__FILE__)+" line "+std::to_string(__LINE__));
				}
			case 'C':
				switch(std::toupper(trans))
				{
					case 'N':	case 'R':	return rows;
					case 'T':	case 'C':	return cols;
					default:	throw std::invalid_argument("trans cannot be "+std::to_string(trans)+". "+std::string(__FILE__)+" line "+std::to_string(__LINE__));
				}
			default:	throw std::invalid_argument("ordering cannot be "+std::to_string(ordering)+". "+std::string(__FILE__)+" line "+std::to_string(__LINE__));
		}
	}
	template<typename T>
	inline void imatcopy (const char ordering, const char trans, size_t rows, size_t cols, const T alpha, T * AB)
	{
		const size_t lda = get_lda_matcopy(ordering, rows, cols);
		const size_t ldb = get_ldb_matcopy(ordering, trans, rows, cols);
		imatcopy (ordering, trans, rows, cols, alpha, AB, lda, ldb);
	}
	template<typename T>
	inline void omatcopy (char ordering, char trans, size_t rows, size_t cols, const T alpha, const T * A, T * B)
	{
		const size_t lda = get_lda_matcopy(ordering, rows, cols);
		const size_t ldb = get_ldb_matcopy(ordering, trans, rows, cols);
		omatcopy (ordering, trans, rows, cols, alpha, A, lda, B, ldb);
	}
}

#endif

}