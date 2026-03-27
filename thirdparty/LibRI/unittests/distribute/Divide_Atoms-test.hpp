// ===================
//  Author: Peize Lin
//  date: 2022.07.13
// ===================

#pragma once

#include "RI/distribute/Divide_Atoms.h"
#include "unittests/print_stl.h"

namespace Divide_Atoms_Test
{
	static void test_divide_atoms()
	{
		const std::size_t group_size = 6;
		std::vector<std::size_t> atoms(31);
		for(std::size_t i=0; i<atoms.size(); ++i)
			atoms[i]=i;
		for(std::size_t i=0; i<group_size; ++i)
			std::cout<<RI::Divide_Atoms::divide_atoms(i, group_size, atoms)<<std::endl;
	}
	/*
		0|	1|	2|	3|	4|	5|
		6|	7|	8|	9|	10|
		11|	12|	13|	14|	15|
		16|	17|	18|	19|	20|
		21|	22|	23|	24|	25|
		26|	27|	28|	29|	30|
	*/


	static void test_divide_atoms_with_period()
	{
		const std::size_t group_size = 6;
		std::vector<std::size_t> atoms(31);
		for(std::size_t i=0; i<atoms.size(); ++i)
			atoms[i]=i;
		const std::array<std::size_t,1> period = {2};
		for(std::size_t i=0; i<group_size; ++i)
			std::cout<<RI::Divide_Atoms::divide_atoms(i, group_size, atoms, period)<<std::endl;
	}
	/*
		{ 0, 0	 }|	{ 0, 1	 }|	{ 1, 0	 }|	{ 1, 1	 }|	{ 2, 0	 }|	{ 2, 1	 }|	{ 3, 0	 }|	{ 3, 1	 }|	{ 4, 0	 }|	{ 4, 1	 }|	{ 5, 0	 }|	{ 5, 1	 }|
		{ 6, 0	 }|	{ 6, 1	 }|	{ 7, 0	 }|	{ 7, 1	 }|	{ 8, 0	 }|	{ 8, 1	 }|	{ 9, 0	 }|	{ 9, 1	 }|	{ 10, 0	 }|	{ 10, 1	 }|
		{ 11, 0	 }|	{ 11, 1	 }|	{ 12, 0	 }|	{ 12, 1	 }|	{ 13, 0	 }|	{ 13, 1	 }|	{ 14, 0	 }|	{ 14, 1	 }|	{ 15, 0	 }|	{ 15, 1	 }|
		{ 16, 0	 }|	{ 16, 1	 }|	{ 17, 0	 }|	{ 17, 1	 }|	{ 18, 0	 }|	{ 18, 1	 }|	{ 19, 0	 }|	{ 19, 1	 }|	{ 20, 0	 }|	{ 20, 1	 }|
		{ 21, 0	 }|	{ 21, 1	 }|	{ 22, 0	 }|	{ 22, 1	 }|	{ 23, 0	 }|	{ 23, 1	 }|	{ 24, 0	 }|	{ 24, 1	 }|	{ 25, 0	 }|	{ 25, 1	 }|
		{ 26, 0	 }|	{ 26, 1	 }|	{ 27, 0	 }|	{ 27, 1	 }|	{ 28, 0	 }|	{ 28, 1	 }|	{ 29, 0	 }|	{ 29, 1	 }|	{ 30, 0	 }|	{ 30, 1	 }|
	*/


	static void test_divide_atoms_periods()
	{
		const std::size_t group_size = 6;
		std::vector<std::size_t> atoms(31);
		for(std::size_t i=0; i<atoms.size(); ++i)
			atoms[i]=i;
		const std::array<std::size_t,1> period = {2};
		for(std::size_t i=0; i<group_size; ++i)
			std::cout<<RI::Divide_Atoms::divide_atoms_periods(i, group_size, atoms, period)<<std::endl;
	}
	/*
		{ 0, 0	 }|	{ 0, 1	 }|	{ 1, 0	 }|	{ 1, 1	 }|	{ 2, 0	 }|	{ 2, 1	 }|	{ 3, 0	 }|	{ 3, 1	 }|	{ 4, 0	 }|	{ 4, 1	 }|	{ 5, 0	 }|
		{ 5, 1	 }|	{ 6, 0	 }|	{ 6, 1	 }|	{ 7, 0	 }|	{ 7, 1	 }|	{ 8, 0	 }|	{ 8, 1	 }|	{ 9, 0	 }|	{ 9, 1	 }|	{ 10, 0	 }|	{ 10, 1	 }|
		{ 11, 0	 }|	{ 11, 1	 }|	{ 12, 0	 }|	{ 12, 1	 }|	{ 13, 0	 }|	{ 13, 1	 }|	{ 14, 0	 }|	{ 14, 1	 }|	{ 15, 0	 }|	{ 15, 1	 }|
		{ 16, 0	 }|	{ 16, 1	 }|	{ 17, 0	 }|	{ 17, 1	 }|	{ 18, 0	 }|	{ 18, 1	 }|	{ 19, 0	 }|	{ 19, 1	 }|	{ 20, 0	 }|	{ 20, 1	 }|
		{ 21, 0	 }|	{ 21, 1	 }|	{ 22, 0	 }|	{ 22, 1	 }|	{ 23, 0	 }|	{ 23, 1	 }|	{ 24, 0	 }|	{ 24, 1	 }|	{ 25, 0	 }|	{ 25, 1	 }|
		{ 26, 0	 }|	{ 26, 1	 }|	{ 27, 0	 }|	{ 27, 1	 }|	{ 28, 0	 }|	{ 28, 1	 }|	{ 29, 0	 }|	{ 29, 1	 }|	{ 30, 0	 }|	{ 30, 1	 }|
	*/
}