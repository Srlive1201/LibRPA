#pragma once

#include "../../ri/Filter_Atom.h"
#include "Symmetry_Filter.h"
#include "../../ri/Label_Tools.h"

#include <stdexcept>

namespace RI
{

template<typename TA, typename TC, typename Tdata>
class Filter_Atom_Symmetry: public Filter_Atom<TA, std::pair<TA, TC>>
{
  public:
	using TAC = std::pair<TA, TC>;

	Filter_Atom_Symmetry(
		const TC &period,
		const std::map<std::pair<TA,TA>, std::set<TC>> &irsec)
		:symmetry(period, irsec){}

	virtual bool filter_for1(const Label::ab_ab &label, const TA &A1) const override
	{
		switch(label)
		{
			// a01b01_a2b01:                                      Aa01
				case Label::ab_ab::a0b0_a2b1:
				case Label::ab_ab::a0b1_a2b0:
				case Label::ab_ab::a1b0_a2b1:
				case Label::ab_ab::a1b1_a2b0:
			// a01b01_a2b2:                                       Aa01
				case Label::ab_ab::a0b0_a2b2:
				case Label::ab_ab::a0b1_a2b2:
				case Label::ab_ab::a1b0_a2b2:
				case Label::ab_ab::a1b1_a2b2:
			// a01b2_a2b01:                                       Aa01
				case Label::ab_ab::a0b2_a2b0:
				case Label::ab_ab::a0b2_a2b1:
				case Label::ab_ab::a1b2_a2b0:
				case Label::ab_ab::a1b2_a2b1:
					return !this->symmetry.is_Aa_in_irreducible_sector(A1);
			default:
				throw std::invalid_argument("label "+Label_Tools::get_name(label)+" error in "+std::string(__FILE__)+" line "+std::to_string(__LINE__));
		}
	}
	virtual bool filter_for1(const Label::ab_ab &label, const TAC &A1) const override
	{
		switch(label)
		{
			// a01b01_a01b01:                                      Aa2
				case Label::ab_ab::a0b0_a1b1:
				case Label::ab_ab::a0b1_a1b0:
					return !this->symmetry.is_Aa_in_irreducible_sector(A1.first);
			// a01b01_a01b2:                                       Ab01
				case Label::ab_ab::a0b0_a1b2:
				case Label::ab_ab::a0b1_a1b2:
				case Label::ab_ab::a0b2_a1b0:
				case Label::ab_ab::a0b2_a1b1:
					return !this->symmetry.is_Ab_in_irreducible_sector(A1.first);
			default:
				throw std::invalid_argument("label "+Label_Tools::get_name(label)+" error in "+std::string(__FILE__)+" line "+std::to_string(__LINE__));
		}
	}

	virtual bool filter_for2(const Label::ab_ab &label, const TA &A1, const TAC &A2) const override
	{
		switch(label)
		{
			// a01b2_a2b01:                                       Aa01,          Ab01
				case Label::ab_ab::a0b2_a2b0:
				case Label::ab_ab::a0b2_a2b1:
				case Label::ab_ab::a1b2_a2b0:
				case Label::ab_ab::a1b2_a2b1:
					return !this->symmetry.in_irreducible_sector(A1, A2);
			// a01b01_a2b01:                                      Aa01,          Ab01
			// a01b01_a2b2:                                       Aa01,          Ab2
			default:
					return false;
		}
	}

	virtual bool filter_for32(const Label::ab_ab &label, const TA &A1, const TAC &A2, const TAC &A3) const override
	{
		switch(label)
		{
			// a01b01_a2b01:                                       Aa01,          Ab01,          Ab2
				case Label::ab_ab::a0b0_a2b1:
				case Label::ab_ab::a0b1_a2b0:
				case Label::ab_ab::a1b0_a2b1:
				case Label::ab_ab::a1b1_a2b0:
			// a01b01_a2b2:                                        Aa01,          Ab2,           Ab01
				case Label::ab_ab::a0b0_a2b2:
				case Label::ab_ab::a0b1_a2b2:
				case Label::ab_ab::a1b0_a2b2:
				case Label::ab_ab::a1b1_a2b2:
					return !this->symmetry.in_irreducible_sector(A1, A3);
			// a01b2_a2b01:                                        Aa01,          Ab01,          Aa2
			default:
					return false;
		}
	}
	virtual bool filter_for32(const Label::ab_ab &label, const TAC &A1, const TA &A2, const TAC &A3) const override
	{
		switch(label)
		{
			// a01b01_a01b2:                                        Ab01,         Aa01,          Aa2
				case Label::ab_ab::a0b0_a1b2:
				case Label::ab_ab::a0b1_a1b2:
				case Label::ab_ab::a0b2_a1b0:
				case Label::ab_ab::a0b2_a1b1:
					return !this->symmetry.in_irreducible_sector(A3, A1);
			default:
				throw std::invalid_argument("label "+Label_Tools::get_name(label)+" error in "+std::string(__FILE__)+" line "+std::to_string(__LINE__));
		}
	}
	virtual bool filter_for32(const Label::ab_ab &label, const TAC &A1, const TAC &A2, const TAC &A3) const override
	{
		switch(label)
		{
			// a01b01_a01b01:                                       Aa2,           Ab01,          Ab2
				case Label::ab_ab::a0b0_a1b1:
				case Label::ab_ab::a0b1_a1b0:
					return !this->symmetry.in_irreducible_sector(A1, A3);
			default:
					throw std::invalid_argument("label "+Label_Tools::get_name(label)+" error in "+std::string(__FILE__)+" line "+std::to_string(__LINE__));
		}
	}

  public:	// private
	Symmetry_Filter<TA,TC,Tdata> symmetry;
};

}