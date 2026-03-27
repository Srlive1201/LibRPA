// ===================
//  Author: Peize Lin
//  date: 2022.12.27
// ===================

#pragma once

#include "Tensor.h"

#include <map>

namespace RI
{

namespace Tensor_Algorithm
{
	// m must be a Hermitian positive definite matrix
	// paras = {"uplo":"U"/"L"}
	template<typename T>
	extern Tensor<T> inverse_matrix_potri(
		const Tensor<T> &m_in,
		const std::map<std::string,std::string> &paras = {{"uplo","U"}});

	// m must be a Hermitian matrix
	// paras = {"uplo":"U"/"L", "absolute_eigen_value_threshold":"", "relative_eigen_value_threshold":""}
	template<typename T>
	extern Tensor<T> inverse_matrix_heev(
		const Tensor<T> &m_in,
		const std::map<std::string,std::string> &paras = {{"uplo","U"}});
}

}

#include "Tensor_Algorithm.hpp"
