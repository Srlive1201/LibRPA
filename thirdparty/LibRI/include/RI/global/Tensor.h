// ===================
//  Author: Peize Lin
//  date: 2021.10.31
// ===================

#pragma once

#include "Shape_Vector.h"
#include "Global_Func-2.h"
#include <memory>
#include <vector>
#include <valarray>
#include <initializer_list>
#include <cassert>

// Attention: Unless using t2=t1.copy() explicitly, t2 and t1 share the same memory for all t2=t1, t2(t1), ...

namespace RI
{

template<typename T>
class Tensor
{
public:
	Shape_Vector shape;
	std::shared_ptr<std::valarray<T>> data=nullptr;

	explicit inline Tensor (const Shape_Vector &shape_in);
	explicit inline Tensor (const Shape_Vector &shape_in, std::shared_ptr<std::valarray<T>> data_in);

	Tensor()=default;
	Tensor(const Tensor<T> &t_in)=default;
	Tensor(Tensor<T> &&t_in)=default;
	Tensor<T> &operator=(const Tensor<T> &t_in)=default;
	Tensor<T> &operator=(Tensor<T> &&t_in)=default;

	inline Tensor reshape (const Shape_Vector &shape_in) const;

	Tensor copy() const;

	inline T& operator() (const std::size_t i0) const;
	inline T& operator() (const std::size_t i0, const std::size_t i1) const;
	inline T& operator() (const std::size_t i0, const std::size_t i1, const std::size_t i2) const;
	inline T& operator() (const std::size_t i0, const std::size_t i1, const std::size_t i2, const std::size_t i3) const;

	Tensor transpose() const;
	Tensor dagger() const;

	// ||d||_p = (|d_1|^p+|d_2|^p+...)^{1/p}
	// if(p==std::numeric_limits<double>::max())    ||d||_max = max_i |d_i|
	Global_Func::To_Real_t<T> norm(const double p) const;

	T* ptr()const{ return &(*this->data)[0]; }

	bool empty() const { return shape.empty(); }

	Tensor & operator += (const Tensor &);
	Tensor & operator -= (const Tensor &);
	Tensor operator-() const;

	template <class Archive> void serialize( Archive & ar ){ ar(shape, data); }		// for cereal

	// Deprecated. Using shape.get_shape_all() instead
	std::size_t get_shape_all() const { return this->shape.get_shape_all(); }
};


template<typename T>
extern Tensor<T> operator+ (const Tensor<T> &t1, const Tensor<T> &t2);
template<typename T>
extern Tensor<T> operator- (const Tensor<T> &t1, const Tensor<T> &t2);

template<typename T>
extern Tensor<T> operator* (const Tensor<T> &t1, const Tensor<T> &t2);
template<typename T>
extern Tensor<T> operator* (const T &t1, const Tensor<T> &t2);
template<typename T>
extern Tensor<T> operator* (const Tensor<T> &t1, const T &t2);


namespace Global_Func
{
	template< typename Tout, typename Tin >
	Tensor<Tout> convert(const Tensor<Tin> &t);
}


template<typename T, std::size_t N0>
extern Tensor<T> to_Tensor(const std::array<T,N0> &a);
template<typename T, std::size_t N0, std::size_t N1>
extern Tensor<T> to_Tensor(const std::array<std::array<T,N1>,N0> &a);
template<typename T, std::size_t N0, std::size_t N1, std::size_t N2>
extern Tensor<T> to_Tensor(const std::array<std::array<std::array<T,N2>,N1>,N0> &a);
template<typename T, std::size_t N0, std::size_t N1, std::size_t N2, std::size_t N3>
extern Tensor<T> to_Tensor(const std::array<std::array<std::array<std::array<T,N3>,N2>,N1>,N0> &a);

template<typename T, std::size_t N0>
extern std::array<T,N0> to_array(const Tensor<T> &t);
template<typename T, std::size_t N0, std::size_t N1>
extern std::array<std::array<T,N1>,N0> to_array(const Tensor<T> &t);
template<typename T, std::size_t N0, std::size_t N1, std::size_t N2>
extern std::array<std::array<std::array<T,N2>,N1>,N0> to_array(const Tensor<T> &t);
template<typename T, std::size_t N0, std::size_t N1, std::size_t N2, std::size_t N3>
extern std::array<std::array<std::array<std::array<T,N3>,N2>,N1>,N0> to_array(const Tensor<T> &t);
}

#include "Blas_Interface-Tensor.h"