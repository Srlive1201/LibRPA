// ===================
//  Author: Peize Lin
//  date: 2021.10.31
// ===================

#pragma once

#include "Blas-Fortran.h"

#include <string>
#include <stdexcept>

#ifdef __MKL_RI
#include <mkl_trans.h>
#endif

namespace RI
{

namespace Blas_Interface
{
	inline char change_uplo(const char &uplo)
	{
		switch(uplo)
		{
			case 'U': return 'L';
			case 'L': return 'U';
			default: throw std::invalid_argument("trans must be 'U' or 'L'.\nIn "+std::string(__FILE__)+" line "+std::to_string(__LINE__)+". ");
		}
	}

	inline char change_trans_NT(const char &trans)
	{
		switch(trans)
		{
			case 'N': return 'T';
			case 'T': return 'N';
			default: throw std::invalid_argument("trans must be 'N' or 'T'.\nIn "+std::string(__FILE__)+" line "+std::to_string(__LINE__)+". ");
		}
	}

	// nrm2 = ||x||_2
	inline float nrm2(const int n, const float*const X, const int incX)
	{
		return snrm2_(&n, X, &incX);
	}
	inline double nrm2(const int n, const double*const X, const int incX)
	{
		return dnrm2_(&n, X, &incX);
	}
	inline float nrm2(const int n, const std::complex<float>*const X, const int incX)
	{
		return scnrm2_(&n, X, &incX);
	}
	inline double nrm2(const int n, const std::complex<double>*const X, const int incX)
	{
		return dznrm2_(&n, X, &incX);
	}

	// Vy = alpha * Vx + Vy
	inline void axpy(const int &n, const float &alpha, const float*const X, const int &incX, float*const Y, const int &incY)
	{
		return saxpy_(&n, &alpha, X, &incX, Y, &incY);
	}
	inline void axpy(const int &n, const double &alpha, const double*const X, const int &incX, double*const Y, const int &incY)
	{
		return daxpy_(&n, &alpha, X, &incX, Y, &incY);
	}
	inline void axpy(const int &n, const std::complex<float> &alpha, const std::complex<float>*const X, const int &incX, std::complex<float>*const Y, const int &incY)
	{
		return caxpy_(&n, &alpha, X, &incX, Y, &incY);
	}
	inline void axpy(const int &n, const std::complex<double> &alpha, const std::complex<double>*const X, const int &incX, std::complex<double>*const Y, const int &incY)
	{
		return zaxpy_(&n, &alpha, X, &incX, Y, &incY);
	}

	// d = Vx * Vy
	inline float dot(const int n, const float*const X, const int incX, const float*const Y, const int incY)
	{
		return sdot_(&n, X, &incX, Y, &incY);
	}
	inline double dot(const int n, const double*const X, const int incX, const double*const Y, const int incY)
	{
		return ddot_(&n, X, &incX, Y, &incY);
	}

	// d = Vx * Vy
	template<typename Tdata,
		typename std::enable_if<std::is_floating_point<Tdata>::value,bool>::type=0>
	inline Tdata dotu(const int n, const Tdata*const X, const int incX, const Tdata*const Y, const int incY)
	{
		return dot(n, X, incX, Y, incY);
	}
	template<typename Tdata,
		typename std::enable_if<std::is_floating_point<Tdata>::value,bool>::type=0>
	inline std::complex<Tdata> dotu(const int n, const std::complex<Tdata>*const X, const int incX, const std::complex<Tdata>*const Y, const int incY)
	{
		//cdotu_(&result, &n, X, &incX, Y, &incY);
		//zdotu_(&result, &n, X, &incX, Y, &incY);

		const int incX2 = 2 * incX;
		const int incY2 = 2 * incY;
		const Tdata*const x = reinterpret_cast<const Tdata*const>(X);
		const Tdata*const y = reinterpret_cast<const Tdata*const>(Y);
		//Re(result)=Re(x)*Re(y)-Im(x)*Im(y)
		//Im(result)=Re(x)*Im(y)+Im(x)*Re(y)
		return std::complex<Tdata>(
			dot(n, x, incX2, y,   incY2) - dot(n, x+1, incX2, y+1, incY2),
			dot(n, x, incX2, y+1, incY2) + dot(n, x+1, incX2, y,   incY2));
	}

	// d = Vx.conj() * Vy
	template<typename Tdata,
		typename std::enable_if<std::is_floating_point<Tdata>::value,bool>::type=0>
	inline Tdata dotc(const int n, const Tdata*const X, const int incX, const Tdata*const Y, const int incY)
	{
		return dot(n, X, incX, Y, incY);
	}
	template<typename Tdata,
		typename std::enable_if<std::is_floating_point<Tdata>::value,bool>::type=0>
	inline std::complex<Tdata> dotc(const int n, const std::complex<Tdata>*const X, const int incX, const std::complex<Tdata>*const Y, const int incY)
	{
		//cdotc_(&result, &n, X, &incX, Y, &incY);
		//zdotc_(&result, &n, X, &incX, Y, &incY);

		const int incX2 = 2 * incX;
		const int incY2 = 2 * incY;
		const Tdata*const x = reinterpret_cast<const Tdata*const>(X);
		const Tdata*const y = reinterpret_cast<const Tdata*const>(Y);
		// Re(result)=Re(X)*Re(Y)+Im(X)*Im(Y)
		// Im(result)=Re(X)*Im(Y)-Im(X)*Re(Y)
		return  std::complex<Tdata>(
			dot(n, x, incX2, y,   incY2) + dot(n, x+1, incX2, y+1, incY2),
			dot(n, x, incX2, y+1, incY2) - dot(n, x+1, incX2, y,   incY2));
	}

	// Vy = alpha * Ma.? * Vx + beta * Vy
	inline void gemv(const char transA, const int m, const int n,
		const float alpha, const float*const A, const int ldA, const float*const X, const int incX,
		const float beta, float*const Y, const int incY)
	{
		const char transA_changed = change_trans_NT(transA);
		sgemv_(&transA_changed, &n, &m,
			&alpha, A, &ldA, X, &incX,
			&beta, Y, &incY);
	}
	inline void gemv(const char transA, const int m, const int n,
		const double alpha, const double*const A, const int ldA, const double*const X, const int incX,
		const double beta, double*const Y, const int incY)
	{
		const char transA_changed = change_trans_NT(transA);
		dgemv_(&transA_changed, &n, &m,
			&alpha, A, &ldA, X, &incX,
			&beta, Y, &incY);
	}
	inline void gemv(const char transA, const int m, const int n,
		const std::complex<float> alpha, const std::complex<float>*const A, const int ldA, const std::complex<float>*const X, const int incX,
		const std::complex<float> beta, std::complex<float>*const Y, const int incY)
	{
		const char transA_changed = change_trans_NT(transA);
		cgemv_(&transA_changed, &n, &m,
			&alpha, A, &ldA, X, &incX,
			&beta, Y, &incY);
	}
	inline void gemv(const char transA, const int m, const int n,
		const std::complex<double> alpha, const std::complex<double>*const A, const int ldA, const std::complex<double>*const X, const int incX,
		const std::complex<double> beta, std::complex<double>*const Y, const int incY)
	{
		const char transA_changed = change_trans_NT(transA);
		zgemv_(&transA_changed, &n, &m,
			&alpha, A, &ldA, X, &incX,
			&beta, Y, &incY);
	}

	// Mc = alpha * Ma.? * Mb.? + beta * Mc
	inline void gemm(const char transA, const char transB, const int m, const int n, const int k,
		const float alpha, const float*const A, const int ldA, const float*const B, const int ldB,
		const float beta, float*const C, const int ldC)
	{
		sgemm_(&transB, &transA, &n, &m, &k,
			&alpha, B, &ldB, A, &ldA,
			&beta, C, &ldC);
	}
	inline void gemm(const char transA, const char transB, const int m, const int n, const int k,
		const double alpha, const double*const A, const int ldA, const double*const B, const int ldB,
		const double beta, double*const C, const int ldC)
	{
		dgemm_(&transB, &transA, &n, &m, &k,
			&alpha, B, &ldB, A, &ldA,
			&beta, C, &ldC);
	}
	inline void gemm(const char transA, const char transB, const int m, const int n, const int k,
		const std::complex<float> alpha, const std::complex<float>*const A, const int ldA, const std::complex<float>*const B, const int ldB,
		const std::complex<float> beta, std::complex<float>*const C, const int ldC)
	{
		cgemm_(&transB, &transA, &n, &m, &k,
			&alpha, B, &ldB, A, &ldA,
			&beta, C, &ldC);
	}
	inline void gemm(const char transA, const char transB, const int m, const int n, const int k,
		const std::complex<double> alpha, const std::complex<double>*const A, const int ldA, const std::complex<double>*const B, const int ldB,
		const std::complex<double> beta, std::complex<double>*const C, const int ldC)
	{
		zgemm_(&transB, &transA, &n, &m, &k,
			&alpha, B, &ldB, A, &ldA,
			&beta, C, &ldC);
	}

	// Mc = alpha * Ma   * Ma.T + beta * C
	// Mc = alpha * Ma.T * Ma   + beta * C
	inline void syrk(const char uploC, const char transA, const int n, const int k,
		const double alpha, const double*const A, const int ldA,
		const double beta, double*const C, const int ldC)
	{
		const char uploC_changed = change_uplo(uploC);
		const char transA_changed = change_trans_NT(transA);
		dsyrk_(&uploC_changed, &transA_changed, &n, &k,
			&alpha, A, &ldA,
			&beta, C, &ldC);
	}
}



#ifdef __MKL_RI

namespace Blas_Interface
{
	inline void imatcopy (const char ordering, const char trans, size_t rows, size_t cols, const float alpha, float * AB, size_t lda, size_t ldb)
	{
		mkl_simatcopy (ordering, trans, rows, cols, alpha, AB, lda, ldb);
	}
	inline void imatcopy (const char ordering, const char trans, size_t rows, size_t cols, const double alpha, double * AB, size_t lda, size_t ldb)
	{
		mkl_dimatcopy (ordering, trans, rows, cols, alpha, AB, lda, ldb);
	}
	inline void imatcopy (const char ordering, const char trans, size_t rows, size_t cols, const std::complex<float> alpha, std::complex<float> * AB, size_t lda, size_t ldb)
	{
		mkl_cimatcopy (ordering, trans, rows, cols, alpha, AB, lda, ldb);
	}
	inline void imatcopy (const char ordering, const char trans, size_t rows, size_t cols, const std::complex<double> alpha, std::complex<double> * AB, size_t lda, size_t ldb)
	{
		mkl_zimatcopy (ordering, trans, rows, cols, alpha, AB, lda, ldb);
	}

	inline void omatcopy (char ordering, char trans, size_t rows, size_t cols, const float alpha, const float * A, size_t lda, float * B, size_t ldb)
	{
		mkl_somatcopy (ordering, trans, rows, cols, alpha, A, lda, B, ldb);
	}
	inline void omatcopy (char ordering, char trans, size_t rows, size_t cols, const double alpha, const double * A, size_t lda, double * B, size_t ldb)
	{
		mkl_domatcopy (ordering, trans, rows, cols, alpha, A, lda, B, ldb);
	}
	inline void omatcopy (char ordering, char trans, size_t rows, size_t cols, const std::complex<float> alpha, const std::complex<float> * A, size_t lda, std::complex<float> * B, size_t ldb)
	{
		mkl_comatcopy (ordering, trans, rows, cols, alpha, A, lda, B, ldb);
	}
	inline void omatcopy (char ordering, char trans, size_t rows, size_t cols, const std::complex<double> alpha, const std::complex<double> * A, size_t lda, std::complex<double> * B, size_t ldb)
	{
		mkl_zomatcopy (ordering, trans, rows, cols, alpha, A, lda, B, ldb);
	}
}

#endif

}