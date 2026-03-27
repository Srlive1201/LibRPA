// ===================
//  Author: Peize Lin
//  date: 2022.06.02
// ===================

#pragma once

#include <complex>
#include <type_traits>

namespace RI
{

namespace Global_Func
{
	template<typename T> struct To_Real	{ using type=T; };
	template<typename T> struct To_Real<std::complex<T>> { using type=T; };
	template<typename T> using To_Real_t = typename Global_Func::To_Real<T>::type;

	template<typename T> struct To_Complex	{ using type=std::complex<T>; };
	template<typename T> struct To_Complex<std::complex<T>> { using type=std::complex<T>; };
	template<typename T> using To_Complex_t = typename To_Complex<T>::type;



	template<typename> struct is_complex_helper : std::false_type {};
	template<typename T> struct is_complex_helper<std::complex<T>> : std::true_type {};
	template<typename T> struct is_complex : is_complex_helper<typename std::remove_const<typename std::remove_reference<T>::type>::type> {};

	// t = convert(t)
	template<
		typename Tout, typename Tin,
		typename std::enable_if<std::is_same<Tin,Tout>::value,int>::type =0>
	Tout convert(const Tin &t)
	{ return t; }

	// complex = convert(real)
	template<
		typename Tout, typename Tin,
		typename std::enable_if<!Global_Func::is_complex<Tin >::value,int>::type =0,
		typename std::enable_if< Global_Func::is_complex<Tout>::value,int>::type =0>
	Tout convert(const Tin &t)
	{ return Tout(t,0); }

	// real = convert(complex)
	template<
		typename Tout, typename Tin,
		typename std::enable_if< Global_Func::is_complex<Tin >::value,int>::type =0,
		typename std::enable_if<!Global_Func::is_complex<Tout>::value,int>::type =0>
	Tout convert(const Tin &t)
	{ return t.real(); }
}

}