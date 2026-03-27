// ===================
//  Author: Peize Lin
//  date: 2021.10.31
// ===================

#pragma once

#include "RI/global/Tensor.h"

namespace RI
{

static std::ostream &operator<<(std::ostream &os, const Shape_Vector &v)
{
	for(std::size_t i=0; i<v.size(); ++i)
		os<<v[i]<<"\t";
	return os;
}

template<typename T>
std::ostream &operator<<(std::ostream &os, const Tensor<T> &t)
{
	switch(t.shape.size())
	{
		case 0:
		{
			os<<std::endl;
			return os;
		}
		case 1:
		{
			for(std::size_t i0=0; i0<t.shape[0]; ++i0)
				os<<t(i0)<<"\t";
			os<<std::endl;
			return os;
		}
		case 2:
		{
			for(std::size_t i0=0; i0<t.shape[0]; ++i0)
			{
				for(std::size_t i1=0; i1<t.shape[1]; ++i1)
	//				os<<t(i0,i1)<<"\t";						// test
					os<<( std::abs(t(i0,i1))>1E-10 ? t(i0,i1) : 0 )<<"\t";
				os<<std::endl;
			}
			return os;
		}
		case 3:
		{
			os<<"["<<std::endl;
			for(std::size_t i0=0; i0<t.shape[0]; ++i0)
			{
				for(std::size_t i1=0; i1<t.shape[1]; ++i1)
				{
					for(std::size_t i2=0; i2<t.shape[2]; ++i2)
						os<<t(i0,i1,i2)<<"\t";
					os<<std::endl;
				}
				os<<std::endl;
			}
			os<<"]"<<std::endl;
			return os;
		}
		case 4:
		{
			os<<"{"<<std::endl;
			for(std::size_t i0=0; i0<t.shape[0]; ++i0)
			{
				os<<"["<<std::endl;
				for(std::size_t i1=0; i1<t.shape[1]; ++i1)
				{
					for(std::size_t i2=0; i2<t.shape[2]; ++i2)
					{
						for(std::size_t i3=0; i3<t.shape[3]; ++i3)
							os<<t(i0,i1,i2,i3)<<"\t";
						os<<std::endl;
					}
					os<<std::endl;
				}
				os<<"]"<<std::endl;
			}
			os<<"}"<<std::endl;
			return os;
		}
		default:
			throw std::invalid_argument(std::string(__FILE__)+" line "+std::to_string(__LINE__));
	}
}

}
//int main1()
//{
//	RI::Tensor<double> t1({2,3});
//	for(int i=0; i<t1.shape[0]; ++i)
//		for(int j=0; j<t1.shape[1]; ++j)
//			t1(i,j) = 10*i+j;
//	std::cout<<t1.shape<<std::endl<<t1<<std::endl;
//	RI::Tensor<double> t2 = t1;
//	t2(1,1)=200;
//	std::cout<<t1.shape<<std::endl<<t1<<std::endl;
//	RI::Tensor<double> t3 = std::move(t2);
//	t3(0,0)=100;
//	std::cout<<t1.shape<<std::endl<<t1<<std::endl;
//	std::cout<<t2.shape<<std::endl<<t2<<std::endl;
//	std::cout<<t3.shape<<std::endl<<t3<<std::endl;
//}