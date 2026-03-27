// ===================
//  Author: Peize Lin
//  date: 2022.06.06
// ===================

#pragma once

#include "RI/ri/Label.h"
#include "RI/ri/Label_Tools.h"
#include <string>
#include <array>
#include <unordered_map>
#include <fstream>

namespace CS_Matrix_Test
{
	template<typename T>
	void print_uplimits3(const std::string &file_name, const std::unordered_map<Label::ab,std::array<T,3>> &datas)
	{
		for(const auto & data : datas)
		{
			for(int i=0; i<3; ++i)
				if(!data.second[i].empty())
				{
					std::ofstream ofs(file_name+"_"+RI::Label_Tools::get_name(data.first)+"_"+std::to_string(i));
					ofs<<data.second[i]<<std::endl;
				}
		}
	}

	template<typename T>
	void print_uplimits2(const std::string &file_name, const std::unordered_map<Label::ab,T> &datas)
	{
		for(const auto & data : datas)
		{
			if(!data.second.empty())
			{
				std::ofstream ofs(file_name+"_"+RI::Label_Tools::get_name(data.first));
				ofs<<data.second<<std::endl;
			}
		}
	}
}