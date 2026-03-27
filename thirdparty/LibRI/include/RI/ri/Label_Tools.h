// ===================
//  Author: Peize Lin
//  date: 2022.10.21
// ===================

#pragma once

#include "Label.h"
#include <string>
#include <vector>
#include <set>
#include <stdexcept>

namespace RI
{

namespace Label_Tools
{
	inline std::string get_name(const Label::ab &label)
	{
		switch(label)
		{
			case Label::ab::a:		return "a";
			case Label::ab::b:		return "b";
			case Label::ab::a0b0:	return "a0b0";
			case Label::ab::a0b1:	return "a0b1";
			case Label::ab::a0b2:	return "a0b2";
			case Label::ab::a1b0:	return "a1b0";
			case Label::ab::a1b1:	return "a1b1";
			case Label::ab::a1b2:	return "a1b2";
			case Label::ab::a2b0:	return "a2b0";
			case Label::ab::a2b1:	return "a2b1";
			case Label::ab::a2b2:	return "a2b2";
			default:	throw std::invalid_argument(std::string(__FILE__)+" line "+std::to_string(__LINE__));
		}
	}

	inline std::string get_name(const Label::ab_ab &label)
	{
		switch(label)
		{
			case Label::ab_ab::a1b1_a2b2:		return "a1b1_a2b2";
			case Label::ab_ab::a1b0_a2b2:		return "a1b0_a2b2";
			case Label::ab_ab::a1b0_a2b1:		return "a1b0_a2b1";
			case Label::ab_ab::a0b1_a2b2:		return "a0b1_a2b2";
			case Label::ab_ab::a0b0_a2b2:		return "a0b0_a2b2";
			case Label::ab_ab::a0b0_a2b1:		return "a0b0_a2b1";
			case Label::ab_ab::a0b1_a1b2:		return "a0b1_a1b2";
			case Label::ab_ab::a0b0_a1b2:		return "a0b0_a1b2";
			case Label::ab_ab::a0b0_a1b1:		return "a0b0_a1b1";
			case Label::ab_ab::a1b2_a2b1:		return "a1b2_a2b1";
			case Label::ab_ab::a1b2_a2b0:		return "a1b2_a2b0";
			case Label::ab_ab::a1b1_a2b0:		return "a1b1_a2b0";
			case Label::ab_ab::a0b2_a2b1:		return "a0b2_a2b1";
			case Label::ab_ab::a0b2_a2b0:		return "a0b2_a2b0";
			case Label::ab_ab::a0b1_a2b0:		return "a0b1_a2b0";
			case Label::ab_ab::a0b2_a1b1:		return "a0b2_a1b1";
			case Label::ab_ab::a0b2_a1b0:		return "a0b2_a1b0";
			case Label::ab_ab::a0b1_a1b0:		return "a0b1_a1b0";
			default:	throw std::invalid_argument(std::string(__FILE__)+" line "+std::to_string(__LINE__));
		}
	}

	inline std::string get_name(const Label::Aab &label)
	{
		switch(label)
		{
			case Label::Aab::a:			return "a";
			case Label::Aab::b:			return "b";
			case Label::Aab::a01b01:	return "a01b01";
			case Label::Aab::a01b2:		return "a01b2";
			case Label::Aab::a2b01:		return "a2b01";
			case Label::Aab::a2b2:		return "a2b2";
			default:	throw std::invalid_argument(std::string(__FILE__)+" line "+std::to_string(__LINE__));
		}
	}

	inline std::string get_name(const Label::Aab_Aab &label)
	{
		switch(label)
		{
			case Label::Aab_Aab::a01b01_a01b01:	return "a01b01_a01b01";
			case Label::Aab_Aab::a01b01_a01b2:	return "a01b01_a01b2";
			case Label::Aab_Aab::a01b01_a2b01:	return "a01b01_a2b01";
			case Label::Aab_Aab::a01b01_a2b2:	return "a01b01_a2b2";
			case Label::Aab_Aab::a01b2_a2b01:	return "a01b2_a2b01";
			default:	throw std::invalid_argument(std::string(__FILE__)+" line "+std::to_string(__LINE__));
		}
	}

	template<typename Tlabel>
	std::string get_name(const std::vector<Tlabel> &label_list)
	{
		std::string name = "";
		for(const auto &label : label_list)
			name += Label_Tools::get_name(label) + "_";
		return name.substr(0, name.size()-1);
	}

	inline int get_a(const Label::ab &label)
	{
		switch(label)
		{
			case Label::ab::a0b0:	case Label::ab::a0b1:	case Label::ab::a0b2:	return 0;
			case Label::ab::a1b0:	case Label::ab::a1b1:	case Label::ab::a1b2:	return 1;
			case Label::ab::a2b0:	case Label::ab::a2b1:	case Label::ab::a2b2:	return 2;
			default:	throw std::invalid_argument(std::string(__FILE__)+" line "+std::to_string(__LINE__));
		}
	}

	inline int get_b(const Label::ab &label)
	{
		switch(label)
		{
			case Label::ab::a0b0:	case Label::ab::a1b0:	case Label::ab::a2b0:	return 0;
			case Label::ab::a0b1:	case Label::ab::a1b1:	case Label::ab::a2b1:	return 1;
			case Label::ab::a0b2:	case Label::ab::a1b2:	case Label::ab::a2b2:	return 2;
			default:	throw std::invalid_argument(std::string(__FILE__)+" line "+std::to_string(__LINE__));
		}
	}

	inline int get_unused_a(const Label::ab_ab &label)
	{
		switch(label)
		{
			case Label::ab_ab::a1b1_a2b2:	case Label::ab_ab::a1b2_a2b1:
			case Label::ab_ab::a1b0_a2b2:	case Label::ab_ab::a1b2_a2b0:
			case Label::ab_ab::a1b0_a2b1:	case Label::ab_ab::a1b1_a2b0:
				return 0;
			case Label::ab_ab::a0b1_a2b2:	case Label::ab_ab::a0b2_a2b1:
			case Label::ab_ab::a0b0_a2b2:	case Label::ab_ab::a0b2_a2b0:
			case Label::ab_ab::a0b0_a2b1:	case Label::ab_ab::a0b1_a2b0:
				return 1;
			case Label::ab_ab::a0b1_a1b2:	case Label::ab_ab::a0b2_a1b1:
			case Label::ab_ab::a0b0_a1b2:	case Label::ab_ab::a0b2_a1b0:
			case Label::ab_ab::a0b0_a1b1:	case Label::ab_ab::a0b1_a1b0:
				return 2;
			default:	throw std::invalid_argument(std::string(__FILE__)+" line "+std::to_string(__LINE__));
		}
	}

	inline int get_unused_b(const Label::ab_ab &label)
	{
		switch(label)
		{
			case Label::ab_ab::a1b1_a2b2:	case Label::ab_ab::a1b2_a2b1:
			case Label::ab_ab::a0b1_a2b2:	case Label::ab_ab::a0b2_a2b1:
			case Label::ab_ab::a0b1_a1b2:	case Label::ab_ab::a0b2_a1b1:
				return 0;
			case Label::ab_ab::a1b0_a2b2:	case Label::ab_ab::a1b2_a2b0:
			case Label::ab_ab::a0b0_a2b2:	case Label::ab_ab::a0b2_a2b0:
			case Label::ab_ab::a0b0_a1b2:	case Label::ab_ab::a0b2_a1b0:
				return 1;
			case Label::ab_ab::a1b0_a2b1:	case Label::ab_ab::a1b1_a2b0:
			case Label::ab_ab::a0b0_a2b1:	case Label::ab_ab::a0b1_a2b0:
			case Label::ab_ab::a0b0_a1b1:	case Label::ab_ab::a0b1_a1b0:
				return 2;
			default:	throw std::invalid_argument(std::string(__FILE__)+" line "+std::to_string(__LINE__));
		}
	}

	inline Label::Aab to_Aab(const Label::ab &label)
	{
		switch(label)
		{
			case Label::ab::a:
				return Label::Aab::a;
			case Label::ab::b:
				return Label::Aab::b;
			case Label::ab::a0b0:	case Label::ab::a0b1:	case Label::ab::a1b0:	case Label::ab::a1b1:
				return Label::Aab::a01b01;
			case Label::ab::a0b2:	case Label::ab::a1b2:
				return Label::Aab::a01b2;
			case Label::ab::a2b0:	case Label::ab::a2b1:
				return Label::Aab::a2b01;
			case Label::ab::a2b2:
				return Label::Aab::a2b2;
			default:	throw std::invalid_argument(std::string(__FILE__)+" line "+std::to_string(__LINE__));
		}
	}

	inline Label::Aab_Aab to_Aab_Aab(const Label::ab_ab &label)
	{
		switch(label)
		{
			case Label::ab_ab::a0b0_a1b1:	case Label::ab_ab::a0b1_a1b0:
				return Label::Aab_Aab::a01b01_a01b01;
			case Label::ab_ab::a0b1_a1b2:	case Label::ab_ab::a0b2_a1b1:	case Label::ab_ab::a0b0_a1b2:	case Label::ab_ab::a0b2_a1b0:
				return Label::Aab_Aab::a01b01_a01b2;
			case Label::ab_ab::a1b0_a2b1:	case Label::ab_ab::a1b1_a2b0:	case Label::ab_ab::a0b0_a2b1:	case Label::ab_ab::a0b1_a2b0:
				return Label::Aab_Aab::a01b01_a2b01;
			case Label::ab_ab::a1b1_a2b2:	case Label::ab_ab::a0b1_a2b2:	case Label::ab_ab::a1b0_a2b2:	case Label::ab_ab::a0b0_a2b2:
				return Label::Aab_Aab::a01b01_a2b2;
			case Label::ab_ab::a1b2_a2b1:	case Label::ab_ab::a0b2_a2b1:	case Label::ab_ab::a1b2_a2b0:	case Label::ab_ab::a0b2_a2b0:
				return Label::Aab_Aab::a01b2_a2b01;
			default:	throw std::invalid_argument(std::string(__FILE__)+" line "+std::to_string(__LINE__));
		}
	}

	inline std::set<Label::Aab_Aab> to_Aab_Aab_set(const std::vector<Label::ab_ab> &labels)
	{
		std::set<Label::Aab_Aab> Aab_Aab_set;
		for(const Label::ab_ab &label : labels)
			Aab_Aab_set.insert(to_Aab_Aab(label));
		return Aab_Aab_set;
	}
}

}