// ===================
//  Author: Peize Lin
//  date: 2022.12.27
// ===================

#pragma once

#include "Tensor_Algorithm.h"
#include "Lapack_Interface-Tensor.h"

#include <string>
#include <stdexcept>
#define LAPACK_INFO_CHECK(x) if(const int info=(x))	throw std::runtime_error("info="+std::to_string(info)+".\n"+std::string(__FILE__)+" line "+std::to_string(__LINE__));

namespace RI
{

namespace Tensor_Algorithm
{
	template<typename T>
	void copy_matrix_triangle(const char &uplo, Tensor<T> &m)
	{
		assert(m.shape.size()==2);
		assert(m.shape[0]==m.shape[1]);
		if(uplo=='U')
		{
			for( std::size_t i0=0; i0!=m.shape[0]; ++i0 )
				for( std::size_t i1=0; i1<i0; ++i1 )
					m(i0,i1) = m(i1,i0);
		}
		else if(uplo=='L')
		{
			for( std::size_t i0=0; i0!=m.shape[0]; ++i0 )
				for( std::size_t i1=i0+1; i1<m.shape[1]; ++i1 )
					m(i0,i1) = m(i1,i0);
		}
		else
		{
			throw std::invalid_argument("uplo must be 'U' or 'L'.\nIn "+std::string(__FILE__)+" line "+std::to_string(__LINE__));
		}
	}


	// m must be a Hermitian positive definite matrix
	// paras = {"uplo":"U"/"L"}
	template<typename T>
	Tensor<T> inverse_matrix_potri(
		const Tensor<T> &m_in,
		const std::map<std::string,std::string> &paras)
	{
		Tensor<T> m = m_in.copy();
		assert(m.shape.size()==2);
		assert(m.shape[0]==m.shape[1]);
		const char uplo = paras.at("uplo").c_str()[0];

		LAPACK_INFO_CHECK(Lapack_Interface::potrf(uplo, m));
		LAPACK_INFO_CHECK(Lapack_Interface::potri(uplo, m));
		copy_matrix_triangle(uplo, m);
		return m;
	}


	// m must be a Hermitian matrix
	// paras = {"uplo":"U"/"L", "absolute_eigen_value_threshold":"", "relative_eigen_value_threshold":""}
	template<typename T>
	Tensor<T> inverse_matrix_heev(
		const Tensor<T> &m_in,
		const std::map<std::string,std::string> &paras)
	{
		Tensor<T> m = m_in.copy();
		assert(m.shape.size()==2);
		assert(m.shape[0]==m.shape[1]);
		using T_real = Global_Func::To_Real_t<T>;
		const char uplo = paras.at("uplo").c_str()[0];

		std::vector<T_real> eigen_values(m.shape[0]);
		LAPACK_INFO_CHECK(Lapack_Interface::heev('V', uplo, m, eigen_values));

		T_real threshold = 0;
		if(paras.find("absolute_eigen_value_threshold")!=paras.end())
		{
			threshold = std::stod(paras.at("absolute_eigen_value_threshold"));
		}
		else if(paras.find("relative_eigen_value_threshold")!=paras.end())
		{
			T_real eigen_value_max = 0;
			for( const T_real &eigen_value : eigen_values )
				eigen_value_max = std::max( std::abs(eigen_value), eigen_value_max );
			threshold = eigen_value_max * std::stod(paras.at("relative_eigen_value_threshold"));
		}

		Tensor<T> em({m.shape[0], m.shape[1]});
		std::size_t ie=0;
		for( std::size_t i=0; i<m.shape[0]; ++i )
			if( std::abs(eigen_values[i]) >= threshold )
			{
				Blas_Interface::axpy(m.shape[1], T(1.0)/eigen_values[i], m.ptr()+i*m.shape[1], em.ptr()+ie*em.shape[1]);
				++ie;
			}
		Tensor<T> mI(m.shape);
		Blas_Interface::gemm( 'C','N', em.shape[1], em.shape[1], ie, T(1.0), em.ptr(), m.ptr(), T(0.0), mI.ptr() );
		return mI;
	}


/*	template<typename T>
	void inverse_matrix_inplace(
		Tensor<T> &m,
		const std::map<std::string,std::string> &paras)
	{
		if(paras.at("method")=="potri")
			inverse_matrix_potri(m, paras);
		else if(paras.at("method")=="heev")
			inverse_matrix_heev(m, paras);
		else
			throw std::invalid_argument(std::string(__FILE__)+" line "+std::to_string(__LINE__));
	}
	template<typename T>
	Tensor<T> inverse_matrix(
		const Tensor<T> &m,
		const std::map<std::string,std::string> &paras)
	{
		// m must be a Hermitian positive definite matrix
		Tensor<T> m_inv = m.copy();
		inverse_matrix_inplace(m_inv, paras);
		return m_inv;
	}
*/
}

}

#undef LAPACK_INFO_CHECK