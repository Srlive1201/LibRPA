#pragma once

#include "../../global/Tensor.h"

namespace RI
{
	namespace Sym
	{
		template <typename T>
		inline void T1_HR(T* TA, const T* A, const Tensor<T>& T1, const int& n2)
		{
			// C' = T1^\dagger * C
			const int& n1 = T1.shape[0];
			Blas_Interface::gemm('C', 'N', n1, n2, n1,
				T(1), T1.ptr(), n1, A, n2, T(0), TA, n2);
			// zgemm_(&notrans, &dagger, &n12, &nabf, &nabf,
			// 	&alpha, A, &n12, T1.ptr(), &nabf, &beta, TA, &n12);
		}

		template <typename T>
		inline void T1_HR_T2(T* TAT, const T* A, const Tensor<T>& T1, const Tensor<T>& T2)
		{
			// H' = T1^\dagger * H * T2
			const int& n2 = T2.shape[0], & n1 = T1.shape[0];
			const RI::Shape_Vector& shape = { static_cast<size_t>(n1),static_cast<size_t>(n2) };
			RI::Tensor<T> AT2(shape);
			Blas_Interface::gemm('N', 'N', n1, n2, n2,
				T(1), A, n2, T2.ptr(), n2, T(0), AT2.ptr(), n2);
			Blas_Interface::gemm('C', 'N', n1, n2, n1,
				T(1), T1.ptr(), n1, AT2.ptr(), n2, T(0), TAT, n2);
			// col-major version
			// zgemm_(&notrans, &notrans, &n2, &n1, &n2,
			// 	&alpha, T2.ptr(), &n2, A, &n2, &beta, AT2.ptr(), &n2);
			// zgemm_(&notrans, &dagger, &n2, &n1, &n1,
			// 	&alpha, AT2.ptr(), &n2, T1.ptr(), &n1, &beta, TAT, &n2);
		}

		template<typename T>
		inline void T1_DR_T2(T* TAT, const T* A, const Tensor<T>& T1, const Tensor<T>& T2)
		{
			// D' = T1^T * D * T2^*  = T1^T * [T2^\dagger * D^T]^T
			const int& n2 = T2.shape[0], & n1 = T1.shape[0];
			const RI::Shape_Vector& shape = { static_cast<size_t>(n1),static_cast<size_t>(n2) };
			RI::Tensor<T> AT2(shape);
			BlasConnector::gemm('C', 'T', n2, n1, n2,
				T(1), T2.ptr(), n2, A, n2, T(0), AT2.ptr(), n1);
			BlasConnector::gemm('T', 'T', n1, n2, n1,
				T(1), T1.ptr(), n1, AT2.ptr(), n1, T(0), TAT, n2);
			// col-major version
			// zgemm_(&transpose, &dagger, &nw1, &nw2, &nw2,
			// 	&alpha, A, &nw2, T2.ptr(), &nw2, &beta, AT2.ptr(), &nw1);
			// zgemm_(&transpose, &transpose, &nw2, &nw1, &nw1,
			// 	&alpha, AT2.ptr(), &nw1, T1.ptr(), &nw1, &beta, TAT, &nw2);
		}

	}

}