// ===================
//  Author: Peize Lin
//  date: 2022.02.18
// ===================

#pragma once

#include <array>

namespace RI
{

namespace Array_Operator
{
	template<typename T, std::size_t N>
	extern std::array<T,N> operator%(const std::array<T,N> &v1, const std::array<T,N> &v2);

	template<typename T, std::size_t N>
	extern std::array<T,N> operator+(const std::array<T,N> &v1, const std::array<T,N> &v2);

	template<typename T, std::size_t N>
	extern std::array<T,N> operator-(const std::array<T,N> &v1, const std::array<T,N> &v2);

	template<typename T, std::size_t N>
	extern std::array<T,N> operator-(const std::array<T,N> &v_in);

	template<typename T, std::size_t N>
	extern std::array<T,N> operator*(const T &s, const std::array<T,N> &v_in);
}

}

#include "Array_Operator.hpp"