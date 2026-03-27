// ===================
//  Author: Peize Lin
//  date: 2022.06.06
// ===================

#pragma once

#include "CS_Matrix.h"
#include "CS_Matrix_Tools.h"
#include "../global/Global_Func-1.h"
#include "../global/Array_Operator.h"

#include <stdexcept>

namespace RI
{

template<typename TA, typename TC, typename Tdata> template<typename T_Data_Wrapper>
typename CS_Matrix<TA,TC,Tdata>::Step CS_Matrix<TA,TC,Tdata>::set_label_A(
	const Label::ab_ab &label_in,
	const TA &Aa01, const TAC &Aa2, const TAC &Ab01, const TAC &Ab2,
	const TC &period,
	const T_Data_Wrapper &data_wrapper) const
{
	using namespace Array_Operator;

	// 每种实现如 LRI_01_01 要写一份
	auto get_Ab = [&Ab01, &Ab2](const int xb) -> TAC
	{
		switch(xb)
		{
			case 0:	case 1:	return Ab01;
			case 2:			return Ab2;
			default:		throw std::invalid_argument("get_Ab");
		}
	};
	auto get_Aa_Ab = [&Aa01, &Aa2, &period, &get_Ab](const int xa, const int xb) -> std::pair<TA,TAC>
	{
		const TAC Ab = get_Ab(xb);
		switch(xa)
		{
			case 0:	case 1:	return {Aa01, Ab};
			case 2:			return {Aa2.first, {Ab.first, (Ab.second-Aa2.second)%period}};
			default:		throw std::invalid_argument("get_Aa_Ab");
		}
	};
	auto get_uplimit_tensor2 = [&get_Aa_Ab, &data_wrapper](const Label::ab &label_ab) -> const Tdata&
	{
		const int xa = Label_Tools::get_a(label_ab);
		const int xb = Label_Tools::get_b(label_ab);
		const std::pair<TA,TAC> &Aa_Ab = get_Aa_Ab(xa,xb);
		const TA &Aa=Aa_Ab.first;
		const TAC &Ab=Aa_Ab.second;
		return Global_Func::find( data_wrapper(label_ab).csm_uplimits.square_tensor2, Aa, Ab);
	};

	Step step;

	step.label = label_in;

	const int index_a = Label_Tools::get_unused_a(step.label);
	step.a_square = Global_Func::find( data_wrapper(Label::ab::a).csm_uplimits.square_tensor3[index_a], Aa01, Aa2 );
	step.a_norm   = Global_Func::find( data_wrapper(Label::ab::a).csm_uplimits.norm_tensor3  [index_a], Aa01, Aa2 );

	const int index_b = Label_Tools::get_unused_b(step.label);
	step.b_square = Global_Func::find( data_wrapper(Label::ab::b).csm_uplimits.square_tensor3[index_b], Ab01.first, TAC{Ab2.first, (Ab2.second-Ab01.second)%period} );
	step.b_norm   = Global_Func::find( data_wrapper(Label::ab::b).csm_uplimits.norm_tensor3  [index_b], Ab01.first, TAC{Ab2.first, (Ab2.second-Ab01.second)%period} );

	const std::pair<Label::ab, Label::ab> label_split_tmp = CS_Matrix_Tools::split_label(step.label);
	step.first_square  = get_uplimit_tensor2( label_split_tmp.first );
	step.second_square = get_uplimit_tensor2( label_split_tmp.second );

	return step;
}

}