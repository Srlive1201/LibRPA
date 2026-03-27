// ===================
//  Author: Peize Lin
//  date: 2023.02.24
// ===================

#pragma once

#include <cereal/cereal.hpp>
#include <initializer_list>
#include <cassert>
#include <numeric>

namespace RI
{

class Shape_Vector
{
public:
	Shape_Vector()=default;
	Shape_Vector(const Shape_Vector &v_in)=default;
	Shape_Vector(Shape_Vector &&v_in)=default;
	Shape_Vector &operator=(const Shape_Vector &v_in)=default;
	Shape_Vector &operator=(Shape_Vector &&v_in)=default;
	Shape_Vector(const std::initializer_list<std::size_t> &v_in)
		:size_(v_in.size())
	{
		assert(v_in.size()<=sizeof(v)/sizeof(*v));
		std::size_t* ptr_this = this->v;
		for(auto ptr_in=v_in.begin(); ptr_in<v_in.end(); )
			*(ptr_this++) = *(ptr_in++);
	}

	const std::size_t* begin() const noexcept { return this->v; }
	const std::size_t* end() const noexcept { return this->v+size_; }
	std::size_t size() const noexcept { return size_; }
	bool empty() const noexcept{ return !size_; }

	std::size_t& operator[] (const std::size_t i)
	{
		assert(i<size_);
		return this->v[i];
	}
	const std::size_t& operator[] (const std::size_t i) const
	{
		assert(i<size_);
		return this->v[i];
	}

	std::size_t get_shape_all() const
	{
		return std::accumulate(this->begin(), this->end(), static_cast<std::size_t>(1), std::multiplies<std::size_t>() );
	}

	template <class Archive> void serialize( Archive & ar ){ ar(cereal::binary_data(this->v,sizeof(v)), size_); }		// for cereal

public:		//private:
	std::size_t v[4];
	std::size_t size_=0;
};

}