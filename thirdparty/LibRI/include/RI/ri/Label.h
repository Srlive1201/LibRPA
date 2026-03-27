// ===================
//  Author: Peize Lin
//  date: 2022.06.02
// ===================

#pragma once

#include <array>
#include <stdexcept>

namespace RI
{

namespace Label
{
	enum class ab
	{
		a, b,
		a0b0, a0b1, a0b2,
		a1b0, a1b1, a1b2,
		a2b0, a2b1, a2b2
	};
	constexpr std::array<ab,11> array_ab =
	{
		ab::a, ab::b,
		ab::a0b0, ab::a0b1, ab::a0b2,
		ab::a1b0, ab::a1b1, ab::a1b2,
		ab::a2b0, ab::a2b1, ab::a2b2
	};

	enum class Aab
	{
		a, b,
		a01b01, a01b2,
		a2b01,  a2b2
	};
	constexpr std::array<Aab,6> array_Aab =
	{
		Aab::a, Aab::b,
		Aab::a01b01, Aab::a01b2,
		Aab::a2b01,  Aab::a2b2
	};

	enum class ab_ab
	{
		a1b1_a2b2, a1b2_a2b1,	// a0,b0
		a1b0_a2b2, a1b2_a2b0,	// a0,b1
		a1b0_a2b1, a1b1_a2b0,	// a0,b2
		a0b1_a2b2, a0b2_a2b1,	// a1,b0
		a0b0_a2b2, a0b2_a2b0,	// a1,b1
		a0b0_a2b1, a0b1_a2b0,	// a1,b2
		a0b1_a1b2, a0b2_a1b1,	// a2,b0
		a0b0_a1b2, a0b2_a1b0,	// a2,b1
		a0b0_a1b1, a0b1_a1b0	// a2,b2
	};
	constexpr std::array<ab_ab,18> array_ab_ab =
	{
		ab_ab::a1b1_a2b2, ab_ab::a1b2_a2b1,
		ab_ab::a1b0_a2b2, ab_ab::a1b2_a2b0,
		ab_ab::a1b0_a2b1, ab_ab::a1b1_a2b0,
		ab_ab::a0b1_a2b2, ab_ab::a0b2_a2b1,
		ab_ab::a0b0_a2b2, ab_ab::a0b2_a2b0,
		ab_ab::a0b0_a2b1, ab_ab::a0b1_a2b0,
		ab_ab::a0b1_a1b2, ab_ab::a0b2_a1b1,
		ab_ab::a0b0_a1b2, ab_ab::a0b2_a1b0,
		ab_ab::a0b0_a1b1, ab_ab::a0b1_a1b0
	};

	enum class Aab_Aab
	{
		a01b01_a01b01,
		a01b01_a01b2,
		a01b01_a2b01,
		a01b01_a2b2,
		a01b2_a2b01,
	};
	constexpr std::array<Aab_Aab,5> array_Aab_Aab =
	{
		Aab_Aab::a01b01_a01b01,
		Aab_Aab::a01b01_a01b2,
		Aab_Aab::a01b01_a2b01,
		Aab_Aab::a01b01_a2b2,
		Aab_Aab::a01b2_a2b01,
	};
}

}