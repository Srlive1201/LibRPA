// ===================
//  Author: Peize Lin
//  date: 2022.06.02
// ===================

#pragma once

#include "LRI.h"
#include "../global/Array_Operator.h"
#include "../global/Global_Func-1.h"
#include "CS_Matrix.h"
#include <memory.h>

template<typename TA, typename Tcell, size_t Ndim, typename Tdata>
auto LRI<TA,Tcell,Ndim,Tdata>::cal(const std::vector<Label::ab_ab> &labels)
-> std::map<TA, std::map<TAC, Tensor<Tdata>>>
{
	using namespace Array_Operator;

	std::map<TA, std::map<TAC, Tensor<Tdata>>> Ds_result;

	auto tensor3_merge = [](const Tensor<Tdata> &D, const bool flag_01_2) -> Tensor<Tdata>
	{
		assert(D.shape.size()==3);
		if(flag_01_2)
			return Tensor<Tdata>( {D.shape[0]*D.shape[1], D.shape[2]}, D.data );
		else
			return Tensor<Tdata>( {D.shape[0], D.shape[1]*D.shape[2]}, D.data );
	};

	auto tensor3_transpose = [](const Tensor<Tdata> &D) -> Tensor<Tdata>
	{
		assert(D.shape.size()==3);
		Tensor<Tdata> D_new({D.shape[1], D.shape[0], D.shape[2]});
		for(size_t i0=0; i0<D.shape[0]; ++i0)
			for(size_t i1=0; i1<D.shape[1]; ++i1)
			{
				memcpy(
					D_new.ptr()+(i1*D.shape[0]+i0)*D.shape[2],
					D.ptr()+(i0*D.shape[1]+i1)*D.shape[2],
					D.shape[2]*sizeof(Tdata));
			}
		return D_new;
	};

	auto get_D_result = [&](
		const Label::ab_ab &label,
		const TA &Aa01, const TAC &Aa2, const TAC &Ab01, const TAC &Ab2)
		-> Tensor<Tdata>&
	{
		const bool flag_a_01 = (Label::get_unused_a(label)!=2);
		const bool flag_b_01 = (Label::get_unused_b(label)!=2);
		if(flag_a_01)
			if(flag_b_01)
				return Ds_result[Aa01][Ab01];
			else
				return Ds_result[Aa01][Ab2];
		else
			if(flag_b_01)
				return Ds_result[Aa2.first][{Ab01.first, (Ab01.second-Aa2.second)%period}];
			else
				return Ds_result[Aa2.first][{Ab2.first, (Ab2.second-Aa2.second)%period}];
	};

	auto add_D = [](const Tensor<Tdata> &D_add, Tensor<Tdata> &D_result)
	{
		if(D_result.empty())
			D_result = D_add;
		else
			//D_result += D_add;
			D_result = D_result + D_add;
	};

	auto judge_x = [](const Label::ab_ab &label) -> int
	{
		switch(label)
		{
			case Label::ab_ab::a1b0_a2b2:	case Label::ab_ab::a1b2_a2b0:
			case Label::ab_ab::a0b0_a2b2:	case Label::ab_ab::a0b2_a2b0:
			case Label::ab_ab::a0b0_a1b2:	case Label::ab_ab::a0b2_a1b0:
				return 0;
			case Label::ab_ab::a1b1_a2b2:	case Label::ab_ab::a1b2_a2b1:
			case Label::ab_ab::a0b1_a2b2:	case Label::ab_ab::a0b2_a2b1:
			case Label::ab_ab::a0b1_a1b2:	case Label::ab_ab::a0b2_a1b1:
				return 1;
			default:
				return -1;
		}
	};

	auto get_abx = [](const Label::ab_ab &label) -> const Label::ab
	{
		switch(label)
		{
			case Label::ab_ab::a0b0_a1b2:	case Label::ab_ab::a0b0_a2b2:	return Label::ab::a0b0;
			case Label::ab_ab::a0b1_a1b2:	case Label::ab_ab::a0b1_a2b2:	return Label::ab::a0b1;
			case Label::ab_ab::a0b2_a1b0:	case Label::ab_ab::a1b0_a2b2:	return Label::ab::a1b0;
			case Label::ab_ab::a0b2_a1b1:	case Label::ab_ab::a1b1_a2b2:	return Label::ab::a1b1;
			case Label::ab_ab::a0b2_a2b0:	case Label::ab_ab::a1b2_a2b0:	return Label::ab::a2b0;
			case Label::ab_ab::a0b2_a2b1:	case Label::ab_ab::a1b2_a2b1:	return Label::ab::a2b1;
			default:	throw std::invalid_argument("get_abx");
		}
	};

	const bool flag_D_b_transpose = [&labels, &judge_x]() -> bool
	{
		for(const Label::ab_ab &label : labels)
			if(judge_x(label)==0)
				return true;
		return false;
	}();

	for(const TA &Aa01 : this->parallel->get_list_Aa01())
	{
		for(const TAC &Aa2 : this->parallel->get_list_Aa2(Aa01))
		{
			const Tensor<Tdata> D_a = Global_Func::find( Ds_ab[Label::ab::a],
				Aa01, Aa2);
			if(D_a.empty())	continue;
			for(const TAC &Ab01 : this->parallel->get_list_Ab01(Aa01, Aa2))
			{
				std::unordered_map<Label::ab_ab, Tensor<Tdata>> Ds_b01;
				std::unordered_map<Label::ab_ab, Tdata_real> Ds_b01_csm;
				for(const TAC &Ab2 : this->parallel->get_list_Ab2(Aa01, Aa2, Ab01))
				{
					const Tensor<Tdata> D_b = Global_Func::find( Ds_ab[Label::ab::b],
						Ab01.first, TAC{Ab2.first, (Ab2.second-Ab01.second)%period});

					if(D_b.empty())	continue;
					const Tensor<Tdata> D_b_transpose = flag_D_b_transpose ? tensor3_transpose(D_b) : Tensor<Tdata>{};

					int x;
					Label::ab label_abx;
					Tensor<Tdata> D_bx;				// D_bx(iby,ibx,ib2)

					//auto F_b01 = [&labels, &Aa01, &Aa2, &Ab01, &Ab2, &Ds_b01, &get_D_result, &add_D, &csm](
					auto F_b01 = [&](
						const Label::ab_ab &label_one,
						const std::function<Tensor<Tdata>()> &cal_D_mul2,
						const std::function<Tensor<Tdata>(const Tensor<Tdata>&)> &cal_D_mul3,
						const std::function<Tensor<Tdata>(const Tensor<Tdata>&)> &cal_D_mul4)
					{
						for(const Label::ab_ab &label : labels)
						{
							if(label==label_one)
							{
								if(this->csm.threshold.at(label))
									this->csm.set_label_A(label, Aa01, Aa2, Ab01, Ab2, period);

								if(this->csm.threshold.at(label))
									if(this->csm.filter_4())
										continue;

								if(Ds_b01[label].empty())
								{
									const Tensor<Tdata> D_mul2 = cal_D_mul2();
									if(this->csm.threshold.at(label))
									{
										const Tdata_real D_mul2_csm = D_mul2.norm(2);
										if(this->csm.filter_3(D_mul2_csm))
											continue;
									}

									Ds_b01[label] = cal_D_mul3(D_mul2);
									if(this->csm.threshold.at(label))
										Ds_b01_csm[label] = Ds_b01[label].norm(2);
								}
								if(this->csm.threshold.at(label))
									if(this->csm.filter_2(Ds_b01_csm[label]))
										continue;
								const Tensor<Tdata> &D_mul3 = Ds_b01[label];

								const Tensor<Tdata> D_mul4 = cal_D_mul4(D_mul3);
								if(this->csm.threshold.at(label))
								{
									const Tdata_real D_mul4_csm = D_mul4.norm(2);
									if(this->csm.filter_1(D_mul4_csm))
										continue;
								}
								Tensor<Tdata> &D_result = get_D_result(label, Aa01, Aa2, Ab01, Ab2);
								add_D(D_mul4, D_result);
							}
						}
					};

					//auto F_bx2 = [&labels, &Aa01, &Aa2, &Ab01, &Ab2, &Ds_b01, &D_b, &D_b_transpose, &x, &label_abx, &D_bx, &judge_x, &get_abx, &get_D_result, &add_D](
					auto F_bx2 = [&](
						const std::set<Label::ab_ab> &label_two,
						const std::function<Tensor<Tdata>()> &cal_D_mul2,
						const std::function<Tensor<Tdata>(const Tensor<Tdata>&)> &cal_D_mul3,
						const std::function<Tensor<Tdata>(const Tensor<Tdata>&)> &cal_D_mul4)
					{
						for(const Label::ab_ab &label : labels)
						{
							if(Global_Func::in_set(label,label_two))
							{
								if(this->csm.threshold.at(label))
									this->csm.set_label_A(label, Aa01, Aa2, Ab01, Ab2, period);

								if(this->csm.threshold.at(label))
									if(this->csm.filter_4())
										continue;

								x = judge_x(label);
								label_abx = get_abx(label);
								// D_bx(iby,ibx,ib2)
								D_bx = (x==1) ? D_b : D_b_transpose;

								if(Ds_b01[label].empty())
								{
									Ds_b01[label] = cal_D_mul2();
									if(this->csm.threshold.at(label))
										Ds_b01_csm[label] = Ds_b01[label].norm(2);
								}
								if(this->csm.threshold.at(label))
									if(this->csm.filter_3(Ds_b01_csm[label]))
										continue;
								const Tensor<Tdata> &D_mul2 = Ds_b01[label];

								const Tensor<Tdata> D_mul3 = cal_D_mul3(D_mul2);
								if(this->csm.threshold.at(label))
								{
									const Tdata_real D_mul3_csm = D_mul3.norm(2);
									if(this->csm.filter_2(D_mul3_csm))	continue;
								}

								const Tensor<Tdata> D_mul4 = cal_D_mul4(D_mul3);
								if(this->csm.threshold.at(label))
								{
									const Tdata_real D_mul4_csm = D_mul4.norm(2);
									if(this->csm.filter_1(D_mul4_csm)) continue;
								}
								Tensor<Tdata> &D_result = get_D_result(label, Aa01, Aa2, Ab01, Ab2);
								add_D(D_mul4, D_result);
							}
						}
					};

					F_b01( Label::ab_ab::a0b0_a1b1,
						// D_mul2(ia1,ia2,ib0) = D_a(ia0,ia1,ia2) * D_ab(ia0,ib0)
						[&](){ return
						Blas_Interface::gemm(
							'T', 'N', 1.0,
							tensor3_merge(D_a,false),
							Ds_ab[Label::ab::a0b0][Aa01][Ab01]
							).reshape({D_a.shape[1], D_a.shape[2], D_b.shape[0]});},
						// D_mul3(ia2,ib0,ib1) = D_mul2(ia1,ia2,ib0) * D_ab(ia1,ib1)
						[&](const Tensor<Tdata> &D_mul2){ return
						Blas_Interface::gemm(
							'T', 'N', 1.0,
							tensor3_merge(D_mul2,false),
							Ds_ab[Label::ab::a1b1][Aa01][Ab01]
							).reshape({D_a.shape[2], D_b.shape[0], D_b.shape[1]});},
						// D_mul4(ia2,ib2) = D_mul3(ia2,ib0,ib1) * D_b(ib0,ib1,ib2);
						[&](const Tensor<Tdata> &D_mul3){ return
						Blas_Interface::gemm(
							'N', 'N', 1.0,
							tensor3_merge(D_mul3,false),
							tensor3_merge(D_b,true));}
					);

					F_b01( Label::ab_ab::a0b1_a1b0,
						// D_mul2(ia0,ia2,ib0) = D_a(ia0,ia1,ia2) * D_ab(ia1,ib0)
						[&](){ return
						Blas_Interface::gemm(
							'T', 'N', 1.0,
							tensor3_merge(tensor3_transpose(D_a),false),
							Ds_ab[Label::ab::a1b0][Aa01][Ab01]
							).reshape({D_a.shape[0], D_a.shape[2], D_b.shape[0]});},
						// D_mul3(ia2,ib0,ib1) = D_mul2(ia0,ia2,ib0) * D_ab(ia0,ib1)
						[&](const Tensor<Tdata> &D_mul2){ return
						Blas_Interface::gemm(
							'T', 'N', 1.0,
							tensor3_merge(D_mul2,false),
							Ds_ab[Label::ab::a0b1][Aa01][Ab01]
							).reshape({D_a.shape[2], D_b.shape[0], D_b.shape[1]});},
						// D_mul4(ia2,ib2) = D_mul3(ia2,ib0,ib1) * D_b(ib0,ib1,ib2);
						[&](const Tensor<Tdata> &D_mul3){ return
						Blas_Interface::gemm(
							'N', 'N', 1.0,
							tensor3_merge(D_mul3,false),
							tensor3_merge(D_b,true));}
					);

					F_b01( Label::ab_ab::a0b0_a2b1,
						// D_mul2(ib0,ia1,ia2) = D_ab(ia0,ib0) * D_a(ia0,ia1,ia2)
						[&](){ return
						Blas_Interface::gemm(
							'T', 'N', 1.0,
							Ds_ab[Label::ab::a0b0][Aa01][Ab01],
							tensor3_merge(D_a,false)
							).reshape({D_b.shape[0], D_a.shape[1], D_a.shape[2]});},
						// D_mul3(ia1,ib0,ib1) = D_mul2(ib0,ia1,ia2) * D_ab(ia2,ib1)
						[&](const Tensor<Tdata> &D_mul2){ return
						Blas_Interface::gemm(
							'N', 'N', 1.0,
							tensor3_merge(tensor3_transpose(D_mul2),true),
							Ds_ab[Label::ab::a2b1][Aa2.first][{Ab01.first, (Ab01.second-Aa2.second)%period}]
							).reshape({D_a.shape[1], D_b.shape[0], D_b.shape[1]});},
						// D_mul4(ia1,ib2) = D_mul3(ia1,ib0,ib1) * D_b(ib0,ib1,ib2);
						[&](const Tensor<Tdata> &D_mul3){ return
						Blas_Interface::gemm(
							'N', 'N', 1.0,
							tensor3_merge(D_mul3,false),
							tensor3_merge(D_b,true));}
					);

					F_b01( Label::ab_ab::a0b1_a2b0,
						// D_mul2(ia0,ia1,ib0) = D_a(ia0,ia1,ia2) * D_ab(ia2,ib0)
						[&](){ return
						Blas_Interface::gemm(
							'N', 'N', 1.0,
							tensor3_merge(D_a,true),
							Ds_ab[Label::ab::a2b0][Aa2.first][{Ab01.first, (Ab01.second-Aa2.second)%period}]
							).reshape({D_a.shape[0], D_a.shape[1], D_b.shape[0]});},
						// D_mul3(ia1,ib0,ib1) = D_mul2(ia0,ia1,ib0) * D_ab(ia0,ib1)
						[&](const Tensor<Tdata> &D_mul2){ return
						Blas_Interface::gemm(
							'T', 'N', 1.0,
							tensor3_merge(D_mul2,false),
							Ds_ab[Label::ab::a0b1][Aa01][Ab01]
							).reshape({D_a.shape[1], D_b.shape[0], D_b.shape[1]});},
						// D_mul4(ia1,ib2) = D_mul3(ia1,ib0,ib1) * D_b(ib0,ib1,ib2);
						[&](const Tensor<Tdata> &D_mul3){ return
						Blas_Interface::gemm(
							'N', 'N', 1.0,
							tensor3_merge(D_mul3,false),
							tensor3_merge(D_b,true));}
					);

					F_b01( Label::ab_ab::a1b0_a2b1,
						// D_mul2(ib1,ia0,ia1) = D_ab(ia2,ib1) * D_a(ia0,ia1,ia2)
						[&](){ return
						Blas_Interface::gemm(
							'T', 'T', 1.0,
							Ds_ab[Label::ab::a2b1][Aa2.first][{Ab01.first,(Ab01.second-Aa2.second)%period}],
							tensor3_merge(D_a,true)
							).reshape({D_b.shape[1], D_a.shape[0], D_a.shape[1]});},
						// D_mul3(ib0,ib1,ia0) = D_ab(ia1,ib0) * D_mul2(ib1,ia0,ia1)
						[&](const Tensor<Tdata> &D_mul2){ return
						Blas_Interface::gemm(
							'T', 'T', 1.0,
							Ds_ab[Label::ab::a1b0][Aa01][Ab01],
							tensor3_merge(D_mul2,true)
							).reshape({D_b.shape[0], D_b.shape[1], D_a.shape[0]});},
						// D_mul4(ia0,ib2) = D_mul3(ib0,ib1,ia0) * D_b(ib0,ib1,ib2);
						[&](const Tensor<Tdata> &D_mul3){ return
						Blas_Interface::gemm(
							'T', 'N', 1.0,
							tensor3_merge(D_mul3,true),
							tensor3_merge(D_b,true));}
					);

					F_b01( Label::ab_ab::a1b1_a2b0,
						// D_mul2(ib0,ia0,ia1) = D_ab(ia2,ib0) * D_a(ia0,ia1,ia2)
						[&](){ return
						Blas_Interface::gemm(
							'T', 'T', 1.0,
							Ds_ab[Label::ab::a2b0][Aa2.first][{Ab01.first,(Ab01.second-Aa2.second)%period}],
							tensor3_merge(D_a,true)
							).reshape({D_b.shape[0], D_a.shape[0], D_a.shape[1]});},
						// D_mul3(ia0,ib0,ib1) = D_mul2(ib0,ia0,ia1) * D_ab(ia1,ib1)
						[&](const Tensor<Tdata> &D_mul2){ return
						Blas_Interface::gemm(
							'N', 'N', 1.0,
							tensor3_merge(tensor3_transpose(D_mul2),true),
							Ds_ab[Label::ab::a1b1][Aa01][Ab01]
							).reshape({D_a.shape[0], D_b.shape[0], D_b.shape[1]});},
						// D_mul4(ia0,ib2) = D_mul3(ia0,ib0,ib1) * D_b(ib0,ib1,ib2);
						[&](const Tensor<Tdata> &D_mul3){ return
						Blas_Interface::gemm(
							'N', 'N', 1.0,
							tensor3_merge(D_mul3,false),
							tensor3_merge(D_b,true));}
					);

					F_bx2( {Label::ab_ab::a0b0_a1b2, Label::ab_ab::a0b1_a1b2},
						// D_mul2(ia1,ia2,ibx) = D_a(ia0,ia1,ia2) * D_ab(ia0,ibx)
						[&](){ return
						Blas_Interface::gemm(
							'T', 'N', 1.0,
							tensor3_merge(D_a,false),
							Ds_ab[label_abx][Aa01][Ab01]
							).reshape({D_a.shape[1], D_a.shape[2], D_b.shape[x]});},
						// D_mul3(ia2,ibx,ib2) = D_mul2(ia1,ia2,ibx) * D_ab(ia1,ib2)
						[&](const Tensor<Tdata> &D_mul2){ return
						Blas_Interface::gemm(
							'T', 'N', 1.0,
							tensor3_merge(D_mul2,false),
							Ds_ab[Label::ab::a1b2][Aa01][Ab2]
							).reshape({D_a.shape[2], D_b.shape[x], D_b.shape[2]});},
						// D_mul4(ia2,iby) = D_mul3(ia2,ibx,ib2) * D_bx(iby,ibx,ib2);
						[&](const Tensor<Tdata> &D_mul3){ return
						Blas_Interface::gemm(
							'N', 'T', 1.0,
							tensor3_merge(D_mul3,false),
							tensor3_merge(D_bx,false));}
					);

					F_bx2( {Label::ab_ab::a0b2_a1b0, Label::ab_ab::a0b2_a1b1},
						// D_mul2(ia0,ia2,ibx) = D_a(ia0,ia1,ia2) * D_ab(ia1,ibx)
						[&](){ return
						Blas_Interface::gemm(
							'T', 'N', 1.0,
							tensor3_merge(tensor3_transpose(D_a),false),
							Ds_ab[label_abx][Aa01][Ab01]
							).reshape({D_a.shape[0], D_a.shape[2], D_b.shape[x]});},
						// D_mul3(ia2,ibx,ib2) = D_mul2(ia0,ia2,ibx) * D_ab(ia0,ib2)
						[&](const Tensor<Tdata> &D_mul2){ return
						Blas_Interface::gemm(
							'T', 'N', 1.0,
							tensor3_merge(D_mul2,false),
							Ds_ab[Label::ab::a0b2][Aa01][Ab2]
							).reshape({D_a.shape[2], D_b.shape[x], D_b.shape[2]});},
						// D_mul4(ia2,iby) = D_mul3(ia2,ibx,ib2) * D_bx(iby,ibx,ib2);
						[&](const Tensor<Tdata> &D_mul3){ return
						Blas_Interface::gemm(
							'N', 'T', 1.0,
							tensor3_merge(D_mul3,false),
							tensor3_merge(D_bx,false));}
					);

					F_bx2( {Label::ab_ab::a0b0_a2b2, Label::ab_ab::a0b1_a2b2},
						// D_mul2(ia1,ibx,ia2) = D_ab(ia0,ibx) * D_a(ia0,ia1,ia2)
						[&](){ return
						tensor3_transpose(Blas_Interface::gemm(
							'T', 'N', 1.0,
							Ds_ab[label_abx][Aa01][Ab01],
							tensor3_merge(D_a,false)
							).reshape({D_b.shape[x], D_a.shape[1], D_a.shape[2]}));},
						// D_mul3(ia1,ibx,ib2) = D_mul2(ia1,ibx,ia2) * D_ab(ia2,ib2)
						[&](const Tensor<Tdata> &D_mul2){ return
						Blas_Interface::gemm(
							'N', 'N', 1.0,
							tensor3_merge(D_mul2,true),
							Ds_ab[Label::ab::a2b2][Aa2.first][{Ab2.first, (Ab2.second-Aa2.second)%period}]
							).reshape({D_a.shape[1], D_b.shape[x], D_b.shape[2]});},
						// D_mul4(ia1,iby) = D_mul3(ia1,ibx,ib2) * D_bx(iby,ibx,ib2);
						[&](const Tensor<Tdata> &D_mul3){ return
						Blas_Interface::gemm(
							'N', 'T', 1.0,
							tensor3_merge(D_mul3,false),
							tensor3_merge(D_bx,false));}
					);

					F_bx2( {Label::ab_ab::a0b2_a2b0, Label::ab_ab::a0b2_a2b1},
						// D_mul2(ia0,ia1,ibx) = D_a(ia0,ia1,ia2) * D_ab(ia2,ibx)
						[&](){ return
						Blas_Interface::gemm(
							'N', 'N', 1.0,
							tensor3_merge(D_a,true),
							Ds_ab[label_abx][Aa2.first][{Ab01.first, (Ab01.second-Aa2.second)%period}]
							).reshape({D_a.shape[0], D_a.shape[1], D_b.shape[x]});},
						// D_mul3(ia1,ibx,ib2) = D_mul2(ia0,ia1,ibx) * D_ab(ia0,ib2)
						[&](const Tensor<Tdata> &D_mul2){ return
						Blas_Interface::gemm(
							'T', 'N', 1.0,
							tensor3_merge(D_mul2,false),
							Ds_ab[Label::ab::a0b2][Aa01][Ab2]
							).reshape({D_a.shape[1], D_b.shape[x], D_b.shape[2]});},
						// D_mul4(ia1,iby) = D_mul3(ia1,ibx,ib2) * D_bx(iby,ibx,ib2);
						[&](const Tensor<Tdata> &D_mul3){ return
						Blas_Interface::gemm(
							'N', 'T', 1.0,
							tensor3_merge(D_mul3,false),
							tensor3_merge(D_bx,false));}
					);

					F_bx2( {Label::ab_ab::a1b0_a2b2, Label::ab_ab::a1b1_a2b2},
						// D_mul2(ia0,ibx,ia2) = D_ab(ia1,ibx) * D_a(ia0,ia1,ia2)
						[&](){ return
						tensor3_transpose(Blas_Interface::gemm(
							'T', 'N', 1.0,
							Ds_ab[label_abx][Aa01][Ab01],
							tensor3_merge(tensor3_transpose(D_a),false)
							).reshape({D_b.shape[x], D_a.shape[0], D_a.shape[2]}));},
						// D_mul3(ia0,ibx,ib2) = D_mul2(ia0,ibx,ia2) * D_ab(ia2,ib2)
						[&](const Tensor<Tdata> &D_mul2){ return
						Blas_Interface::gemm(
							'N', 'N', 1.0,
							tensor3_merge(D_mul2,true),
							Ds_ab[Label::ab::a2b2][Aa2.first][{Ab2.first, (Ab2.second-Aa2.second)%period}]
							).reshape({D_a.shape[0], D_b.shape[x], D_b.shape[2]});},
						// D_mul4(ia0,iby) = D_mul3(ia0,ibx,ib2) * D_bx(iby,ibx,ib2);
						[&](const Tensor<Tdata> &D_mul3){ return
						Blas_Interface::gemm(
							'N', 'T', 1.0,
							tensor3_merge(D_mul3,false),
							tensor3_merge(D_bx,false));}
					);

					F_bx2( {Label::ab_ab::a1b2_a2b0, Label::ab_ab::a1b2_a2b1},
						// D_mul2(ia0,ibx,ia1) = D_ab(ia2,ibx) * D_a(ia0,ia1,ia2)
						[&](){ return
						tensor3_transpose(Blas_Interface::gemm(
							'T', 'T', 1.0,
							Ds_ab[label_abx][Aa2.first][{Ab01.first, (Ab01.second-Aa2.second)%period}],
							tensor3_merge(D_a,true)
							).reshape({D_b.shape[x], D_a.shape[0], D_a.shape[1]}));},
						// D_mul3(ia0,ibx,ib2) = D_mul2(ia0,ibx,ia1) * D_ab(ia1,ib2)
						[&](const Tensor<Tdata> &D_mul2){ return
						Blas_Interface::gemm(
							'N', 'N', 1.0,
							tensor3_merge(D_mul2,true),
							Ds_ab[Label::ab::a1b2][Aa01][Ab2]
							).reshape({D_a.shape[0], D_b.shape[x], D_b.shape[2]});},
						// D_mul4(ia0,iby) = D_mul3(ia0,ibx,ib2) * D_bx(iby,ibx,ib2);
						[&](const Tensor<Tdata> &D_mul3){ return
						Blas_Interface::gemm(
							'N', 'T', 1.0,
							tensor3_merge(D_mul3,false),
							tensor3_merge(D_bx,false));}
					);
				} // for Ab2
			} // for Ab01
		} // for Aa2
	} // for Aa01

	return Ds_result;
}