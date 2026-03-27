// ===================
//  Author: Peize Lin
//  date: 2021.10.31
// ===================

#pragma once

#include "Tensor.h"
#include "Blas_Interface-Tensor.h"

#include <complex>
#include <memory>
#include <vector>
#include <functional>
#include <cassert>
#include <limits>

namespace RI
{

template<typename T>
Tensor<T>::Tensor (const Shape_Vector &shape_in)
{
	this->shape = shape_in;
	if(!this->shape.empty())
	{
		const std::size_t shape_all = this->shape.get_shape_all();
		if(shape_all)
			this->data = std::make_shared<std::valarray<T>>(0,shape_all);
	}
}

template<typename T>
Tensor<T>::Tensor (const Shape_Vector &shape_in, std::shared_ptr<std::valarray<T>> data_in)
{
	assert( std::accumulate(shape_in.begin(), shape_in.end(), static_cast<std::size_t>(1), std::multiplies<std::size_t>() ) == data_in->size() );
	this->shape = shape_in;
	this->data = data_in;
}

template<typename T>
Tensor<T> Tensor<T>::reshape (const Shape_Vector &shape_in) const
{
	assert(
		std::accumulate(shape_in.begin(), shape_in.end(), static_cast<std::size_t>(1), std::multiplies<std::size_t>())
		== this->shape.get_shape_all() );
	return Tensor<T>(shape_in, this->data);
}

template<typename T>
Tensor<T> Tensor<T>::copy() const
{
	Tensor<T> t(this->shape);
	*t.data = *this->data;
	return t;
}


template<typename T>
T& Tensor<T>::operator() (const std::size_t i0) const
{
	assert(this->shape.size()==1);
	assert(i0>=0);	assert(i0<this->shape[0]);
	return (*this->data)[i0];
}
template<typename T>
T& Tensor<T>::operator() (const std::size_t i0, const std::size_t i1) const
{
	assert(this->shape.size()==2);
	assert(i0>=0);	assert(i0<this->shape[0]);
	assert(i1>=0);	assert(i1<this->shape[1]);
	return (*this->data)[i0*this->shape[1]+i1];
}
template<typename T>
T& Tensor<T>::operator() (const std::size_t i0, const std::size_t i1, const std::size_t i2) const
{
	assert(this->shape.size()==3);
	assert(i0>=0);	assert(i0<this->shape[0]);
	assert(i1>=0);	assert(i1<this->shape[1]);
	assert(i2>=0);	assert(i2<this->shape[2]);
	return (*this->data)[(i0*this->shape[1]+i1)*this->shape[2]+i2];
}
template<typename T>
T& Tensor<T>::operator() (const std::size_t i0, const std::size_t i1, const std::size_t i2, const std::size_t i3) const
{
	assert(this->shape.size()==4);
	assert(i0>=0);	assert(i0<this->shape[0]);
	assert(i1>=0);	assert(i1<this->shape[1]);
	assert(i2>=0);	assert(i2<this->shape[2]);
	assert(i3>=0);	assert(i3<this->shape[3]);
	return (*this->data)[((i0*this->shape[1]+i1)*this->shape[2]+i2)*this->shape[3]+i3];
}

template<typename T1, typename T2>
bool same_shape (const Tensor<T1> &t1, const Tensor<T2> &t2)
{
	if(t1.shape.size() != t2.shape.size())
		return false;
	for(std::size_t ishape=0; ishape<t1.shape.size(); ++ishape)
		if(t1.shape[ishape] != t2.shape[ishape])
			return false;
	return true;
}

template<typename T>
Tensor<T> Tensor<T>::operator-() const
{
	Tensor<T> t(this->shape);
	*t.data = -*this->data;
	return t;
}


template<typename T>
Tensor<T> operator+ (const Tensor<T> &t1, const Tensor<T> &t2)
{
	assert(same_shape(t1,t2));
	Tensor<T> t(t1.shape);
	*t.data = *t1.data + *t2.data;
	return t;
}
template<typename T>
Tensor<T> operator- (const Tensor<T> &t1, const Tensor<T> &t2)
{
	assert(same_shape(t1,t2));
	Tensor<T> t(t1.shape);
	*t.data = *t1.data - *t2.data;
	return t;
}

template<typename T>
Tensor<T> operator* (const T &t1, const Tensor<T> &t2)
{
	Tensor<T> t(t2.shape);
	*t.data = t1 * *t2.data;
	return t;
}

template<typename T>
Tensor<T> operator* (const Tensor<T> &t1, const T &t2)
{
	Tensor<T> t(t1.shape);
	*t.data = *t1.data * t2;
	return t;
}

template<typename T>
Tensor<T> &Tensor<T>::operator+= (const Tensor &t)
{
	assert(same_shape(*this,t));
	*this->data += *t.data;
	return *this;
}

template<typename T>
Tensor<T> &Tensor<T>::operator-= (const Tensor &t)
{
	assert(same_shape(*this,t));
	*this->data -= *t.data;
	return *this;
}

template<typename T>
Tensor<T> Tensor<T>::transpose() const
{
	assert(this->shape.size()==2);
	Tensor<T> t({this->shape[1], this->shape[0]});
	for(std::size_t i0=0; i0<this->shape[0]; ++i0)
		for(std::size_t i1=0; i1<this->shape[1]; ++i1)
			t(i1,i0) = (*this)(i0,i1);
	return t;
}

template<typename T>
Tensor<T> Tensor<T>::dagger() const
{
	assert(this->shape.size() == 2);
	Tensor<T> t({ this->shape[1], this->shape[0] });
	for (std::size_t i0 = 0; i0 < this->shape[0]; ++i0)
		for (std::size_t i1 = 0; i1 < this->shape[1]; ++i1)
			t(i1, i0) = std::conj((*this)(i0, i1));
	return t;
}

template<typename T>
Global_Func::To_Real_t<T> Tensor<T>::norm(const double p) const
{
	using T_res = Global_Func::To_Real_t<T>;
	if(p==2)
		return Blas_Interface::nrm2(*this);
	else if(p==1)
	{
		T_res s = 0;
		for(std::size_t i=0; i<this->data->size(); ++i)
			s += std::abs((*this->data)[i]);
		return s;
	}
	else if(p==std::numeric_limits<double>::max())
	{
		T_res s = 0;
		for(std::size_t i=0; i<this->data->size(); ++i)
			s = std::max(std::real(s), std::abs((*this->data)[i]));
		return s;
	}
	else
	{
		T_res s = 0;
		for(std::size_t i=0; i<this->data->size(); ++i)
			s += std::pow(std::abs((*this->data)[i]), p);
		return std::pow(s,1.0/p);
	}
}


namespace Global_Func
{
	template<typename Tout, typename Tin>
	Tensor<Tout> convert(const Tensor<Tin> &t_in)
	{
		Tensor<Tout> t_out(t_in.shape);
		for(std::size_t i=0; i<t_out.data->size(); ++i)
			(*t_out.data)[i] = Global_Func::convert<Tout>((*t_in.data)[i]);
		return t_out;
	}
}


template<typename T, std::size_t N0>
Tensor<T> to_Tensor(const std::array<T,N0> &a)
{
	Tensor<T> t({N0});
	for(std::size_t i0=0; i0<N0; ++i0)
		t(i0) = a[i0];
	return t;
}
template<typename T, std::size_t N0, std::size_t N1>
Tensor<T> to_Tensor(const std::array<std::array<T,N1>,N0> &a)
{
	Tensor<T> t({N0,N1});
	for(std::size_t i0=0; i0<N0; ++i0)
		for(std::size_t i1=0; i1<N1; ++i1)
			t(i0,i1) = a[i0][i1];
	return t;
}
template<typename T, std::size_t N0, std::size_t N1, std::size_t N2>
Tensor<T> to_Tensor(const std::array<std::array<std::array<T,N2>,N1>,N0> &a)
{
	Tensor<T> t({N0,N1,N2});
	for(std::size_t i0=0; i0<N0; ++i0)
		for(std::size_t i1=0; i1<N1; ++i1)
			for(std::size_t i2=0; i2<N2; ++i2)
				t(i0,i1,i2) = a[i0][i1][i2];
	return t;
}
template<typename T, std::size_t N0, std::size_t N1, std::size_t N2, std::size_t N3>
Tensor<T> to_Tensor(const std::array<std::array<std::array<std::array<T,N3>,N2>,N1>,N0> &a)
{
	Tensor<T> t({N0,N1,N2,N3});
	for(std::size_t i0=0; i0<N0; ++i0)
		for(std::size_t i1=0; i1<N1; ++i1)
			for(std::size_t i2=0; i2<N2; ++i2)
				for(std::size_t i3=0; i3<N3; ++i3)
					t(i0,i1,i2,i3) = a[i0][i1][i2][i3];
	return t;
}

template<typename T, std::size_t N0>
std::array<T,N0> to_array(const Tensor<T> &t)
{
	assert(t.shape.size()==1);
	assert(t.shape[0]==N0);
	std::array<T,N0> a;
	for(std::size_t i0=0; i0<N0; ++i0)
		a[i0] = t(i0);
	return a;
}
template<typename T, std::size_t N0, std::size_t N1>
std::array<std::array<T,N1>,N0> to_array(const Tensor<T> &t)
{
	assert(t.shape.size()==2);
	assert(t.shape[0]==N0);
	assert(t.shape[1]==N1);
	std::array<std::array<T,N1>,N0> a;
	for(std::size_t i0=0; i0<N0; ++i0)
		for(std::size_t i1=0; i1<N1; ++i1)
			a[i0][i1] = t(i0,i1);
	return a;
}
template<typename T, std::size_t N0, std::size_t N1, std::size_t N2>
std::array<std::array<std::array<T,N2>,N1>,N0> to_array(const Tensor<T> &t)
{
	assert(t.shape.size()==3);
	assert(t.shape[0]==N0);
	assert(t.shape[1]==N1);
	assert(t.shape[2]==N2);
	std::array<std::array<T,N1>,N0> a;
	for(std::size_t i0=0; i0<N0; ++i0)
		for(std::size_t i1=0; i1<N1; ++i1)
			for(std::size_t i2=0; i2<N2; ++i2)
				a[i0][i1][i2] = t(i0,i1,i2);
	return a;
}
template<typename T, std::size_t N0, std::size_t N1, std::size_t N2, std::size_t N3>
std::array<std::array<std::array<std::array<T,N3>,N2>,N1>,N0> to_array(const Tensor<T> &t)
{
	assert(t.shape.size()==4);
	assert(t.shape[0]==N0);
	assert(t.shape[1]==N1);
	assert(t.shape[2]==N2);
	assert(t.shape[3]==N3);
	std::array<std::array<T,N1>,N0> a;
	for(std::size_t i0=0; i0<N0; ++i0)
		for(std::size_t i1=0; i1<N1; ++i1)
			for(std::size_t i2=0; i2<N2; ++i2)
				for(std::size_t i3=0; i3<N3; ++i3)
					a[i0][i1][i2][i3] = t(i0,i1,i2,i3);
	return a;
}

}