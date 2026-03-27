// ===================
//  Author: Peize Lin
//  date: 2022.07.27
// ===================

#pragma once

#include "RI/global/Map_Operator-2.h"
#include "RI/global/Map_Operator-3.h"
#include "unittests/print_stl.h"

#include <string>

namespace Map_Operator_Test
{
	static void test_union_map1()
	{
		std::map<int,double> m1, m2;
		m1[2]=2;	m1[1]=1;
		m2[0]=0;	m2[1]=10;
		const std::function<double(const double&,const double&)> plus = std::plus<double>();
		std::cout<<RI::Map_Operator::zip_union(m1,m2,plus)<<std::endl;
		/*{
			0: 0,
			1: 11,
			2: 2
		}*/
	}

	static void test_union_map2()
	{
		std::map<int,std::map<std::string,double>> m1, m2;
		m1[2]["a"]=2;	m1[1]["b"]=1;
		m2[0]["c"]=0;	m2[1]["b"]=10;	m1[1]["c"]=20;
		const std::function<double(const double&,const double&)> plus = std::plus<double>();
		std::cout<<RI::Map_Operator::zip_union(m1,m2,plus)<<std::endl;
		/*{
			0: {"c":0},
			1: {"b":11, "c":20},
			2: {"a":2}
		}*/
	}

	static void test_intersection_map1()
	{
		std::map<int,double> m1, m2;
		m1[2]=2;	m1[1]=1;
		m2[0]=0;	m2[1]=10;
		const std::function<double(const double&,const double&)> plus = std::plus<double>();
		std::cout<<RI::Map_Operator::zip_intersection(m1,m2,plus)<<std::endl;
		/*{
			1: 11
		}*/
	}

	static void test_intersection_map2()
	{
		std::map<int,std::map<std::string,double>> m1, m2;
		m1[2]["a"]=2;	m1[1]["b"]=1;
		m2[0]["c"]=0;	m2[1]["b"]=10;	m1[1]["c"]=20;
		const std::function<double(const double&,const double&)> plus = std::plus<double>();
		std::cout<<RI::Map_Operator::zip_intersection(m1,m2,plus)<<std::endl;
		/*{
			1: {"b":11},
		}*/
	}

	static void test_intersection_map3()
	{
		std::map<int,std::map<std::string,double>> m1;
		std::map<int,std::map<std::string,std::string>> m2;
		m1[2]["a"]=2;	m1[1]["b"]=1;
		m2[0]["c"]="0";	m2[1]["b"]="10";	m2[1]["c"]="20";
		const std::function<std::string(const double&, const std::string&)> plus = [](const double &d1, const std::string &d2) -> std::string	{ return std::to_string(d1+std::stod(d2))+"a"; };
		std::cout<<RI::Map_Operator::zip_intersection(m1,m2,plus)<<std::endl;
		/*{
			1: {"b":"11a"}
		}*/
	}

	static void test_transform_map()
	{
		std::map<int,std::map<std::string,double>> m;
		m[0]["c"]=0;	m[1]["b"]=10;	m[1]["c"]=20;
		constexpr double frac = 2;
		const std::function<double(const double&)> multiply_frac = [](const double &t){ return t*frac; };
		std::cout<<RI::Map_Operator::transform(m, multiply_frac)<<std::endl;
		/*{
			0: {"c":0},
			1: {"b":20, "c":40}
		}*/
	}

	static void test_for_each_map()
	{
		std::map<int,std::map<std::string,double>> m;
		m[0]["c"]=0;	m[1]["b"]=10;	m[1]["c"]=20;
		constexpr double frac = 2;
		const std::function<void(double&)> multiply_frac = [](double &t){ t*=frac; };
		RI::Map_Operator::for_each(m, multiply_frac);
		std::cout<<m<<std::endl;
		/*{
			0: {"c":0},
			1: {"b":20, "c":40}
		}*/
	}
}