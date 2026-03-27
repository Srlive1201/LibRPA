// ===================
//  Author: Peize Lin
//  date: 2023.08.02
// ===================

#pragma once

#include "LRI.h"
#include "LRI_Cal_Aux.h"
#include "../global/Array_Operator.h"
#include "../global/Tensor_Multiply.h"

#include <omp.h>
#ifdef __MKL_RI
#include <mkl_service.h>
#endif

namespace RI
{

template<typename TA, typename Tcell, std::size_t Ndim, typename Tdata>
void LRI<TA,Tcell,Ndim,Tdata>::cal_loop3(
	const std::vector<Label::ab_ab> &labels,
	std::map<TA, std::map<TAC, Tensor<Tdata>>> &Ds_result,
	const double fac_add_Ds)
{
	using namespace Array_Operator;

	const Data_Pack_Wrapper<TA,TC,Tdata> data_wrapper(this->data_pool, this->data_ab_name);
	const LRI_Cal_Tools<TA,TC,Tdata> tools(this->period, this->data_pool, this->data_ab_name);

	std::map<TA, std::map<TAC, Tensor<Tdata>>> Ds_a_transpose, Ds_b_transpose;
	std::tie(Ds_a_transpose, Ds_b_transpose) = tools.cal_Ds_transpose(labels);

  #ifdef __MKL_RI
	const std::size_t mkl_threads = mkl_get_max_threads();
	mkl_set_num_threads(1);
  #endif

	std::map<TA, omp_lock_t> lock_Ds_result_add_map = LRI_Cal_Aux::init_lock_result(labels, this->parallel->list_A, Ds_result);

	#pragma omp parallel
	{
		std::map<TA, std::map<TAC, Tensor<Tdata>>> Ds_result_thread;

		for(const Label::ab_ab &label : labels)
		{
			const std::vector<TA>  list_Aa01_Da = LRI_Cal_Aux::filter_list_map( this->parallel->list_A.at(Label_Tools::to_Aab_Aab(label)).a01, data_wrapper(Label::ab::a).Ds_ab );
			const std::vector<TAC> list_Ab01_Db = LRI_Cal_Aux::filter_list_map( this->parallel->list_A.at(Label_Tools::to_Aab_Aab(label)).b01, data_wrapper(Label::ab::b).Ds_ab );
			const std::vector<TAC> list_Aa2_Da  = LRI_Cal_Aux::filter_list_set( this->parallel->list_A.at(Label_Tools::to_Aab_Aab(label)).a2,  data_wrapper(Label::ab::a).index_Ds_ab[0] );
			const std::vector<TAC> list_Ab2_Db  = LRI_Cal_Aux::filter_list_set( this->parallel->list_A.at(Label_Tools::to_Aab_Aab(label)).b2,  data_wrapper(Label::ab::b).index_Ds_ab[0] );
			switch(label)
			{

			  // Aab_Aab::a01b01_a01b01

				case Label::ab_ab::a0b0_a1b1:
				{
					const std::vector<TA >  list_Aa01 = LRI_Cal_Aux::filter_list_map( LRI_Cal_Aux::filter_list_map(
						list_Aa01_Da,
						data_wrapper(Label::ab::a0b0).Ds_ab ),
						data_wrapper(Label::ab::a1b1).Ds_ab );
					const std::vector<TAC> &list_Aa2 =
						list_Aa2_Da;
					const std::vector<TAC>  list_Ab01 = LRI_Cal_Aux::filter_list_set( LRI_Cal_Aux::filter_list_set(
						list_Ab01_Db,
						data_wrapper(Label::ab::a0b0).index_Ds_ab[0]),
						data_wrapper(Label::ab::a1b1).index_Ds_ab[0]);
					const std::vector<TAC> &list_Ab2 =
						list_Ab2_Db;

					for(const TAC &Aa2 : list_Aa2)
					{
						if(this->filter_atom->filter_for1(label,Aa2))	continue;
						std::map<TAC,Tensor<Tdata>> Ds_result_fixed;

						#pragma omp for schedule(dynamic) nowait
						for(std::size_t ib01=0; ib01<list_Ab01.size(); ++ib01)
						{
							const TAC &Ab01 = list_Ab01[ib01];
							if(this->filter_atom->filter_for2(label,Aa2,Ab01))	continue;
							// D_mul = D_a * D_a0b0 * D_a1b1
							Tensor<Tdata> D_mul;
							for(const TA &Aa01 : list_Aa01)
							{
								if(this->filter_atom->filter_for31(label,Aa2,Ab01,Aa01))	continue;
								const Tensor<Tdata> &D_a = tools.get_Ds_ab(Label::ab::a, Aa01, Aa2);
								if(D_a.empty())	continue;
								const Tensor<Tdata> &D_a0b0 = tools.get_Ds_ab(Label::ab::a0b0, Aa01, Ab01);
								if(D_a0b0.empty())	continue;
								const Tensor<Tdata> &D_a1b1 = tools.get_Ds_ab(Label::ab::a1b1, Aa01, Ab01);
								if(D_a1b1.empty())	continue;

								// a1a2b0 = a0a1a2 * a0b0
								const Tensor<Tdata> D_tmp1 = Tensor_Multiply::x1x2y1_ax1x2_ay1(D_a, D_a0b0);
								// a2b0b1 = a1a2b0 * a1b1
								Tensor<Tdata> D_tmp2 = Tensor_Multiply::x1x2y1_ax1x2_ay1(D_tmp1, D_a1b1);
								LRI_Cal_Aux::add_Ds(std::move(D_tmp2), D_mul);
							}
							if(D_mul.empty())	continue;

							// D_result = D_mul * D_b
							for(const TAC &Ab2 : list_Ab2)
							{
								if(this->filter_atom->filter_for32(label,Aa2,Ab01,Ab2))	continue;
								const Tensor<Tdata> &D_b = tools.get_Ds_ab(Label::ab::b, Ab01, Ab2);
								if(D_b.empty())	continue;
								// a2b2 = a2b0b1 * b0b1b2
								Tensor<Tdata> D_tmp3 = Tensor_Multiply::x0y2_x0ab_aby2(D_mul, D_b);
								LRI_Cal_Aux::add_Ds(std::move(D_tmp3), Ds_result_fixed[Ab2]);
							}
						} // end for Ab01

						if(!Ds_result_fixed.empty())
							LRI_Cal_Aux::add_Ds(LRI_Cal_Aux::Ds_translate(std::move(Ds_result_fixed), Aa2.second, this->period),
												Ds_result_thread[Aa2.first]);
						LRI_Cal_Aux::add_Ds_omp_try_map(Ds_result_thread, Ds_result, lock_Ds_result_add_map, fac_add_Ds);
					} // end for Ab2
				} break; // end case a0b0_a1b1

				case Label::ab_ab::a0b1_a1b0:
				{
					const std::vector<TA >  list_Aa01 = LRI_Cal_Aux::filter_list_map( LRI_Cal_Aux::filter_list_map(
						list_Aa01_Da,
						data_wrapper(Label::ab::a0b1).Ds_ab ),
						data_wrapper(Label::ab::a1b0).Ds_ab );
					const std::vector<TAC> &list_Aa2 =
						list_Aa2_Da;
					const std::vector<TAC>  list_Ab01 = LRI_Cal_Aux::filter_list_set( LRI_Cal_Aux::filter_list_set(
						list_Ab01_Db,
						data_wrapper(Label::ab::a0b1).index_Ds_ab[0]),
						data_wrapper(Label::ab::a1b0).index_Ds_ab[0]);
					const std::vector<TAC> &list_Ab2 =
						list_Ab2_Db;

					for(const TAC &Aa2 : list_Aa2)
					{
						if(this->filter_atom->filter_for1(label,Aa2))	continue;
						std::map<TAC,Tensor<Tdata>> Ds_result_fixed;

						#pragma omp for schedule(dynamic) nowait
						for(std::size_t ib01=0; ib01<list_Ab01.size(); ++ib01)
						{
							const TAC &Ab01 = list_Ab01[ib01];
							if(this->filter_atom->filter_for2(label,Aa2,Ab01))	continue;
							// D_mul = D_a * D_a0b1 * D_a1b0
							Tensor<Tdata> D_mul;
							for(const TA &Aa01 : list_Aa01)
							{
								if(this->filter_atom->filter_for31(label,Aa2,Ab01,Aa01))	continue;
								const Tensor<Tdata> &D_a = Global_Func::find(Ds_a_transpose, Aa01, Aa2);
								if(D_a.empty())	continue;
								const Tensor<Tdata> &D_a0b1 = tools.get_Ds_ab(Label::ab::a0b1, Aa01, Ab01);
								if(D_a0b1.empty())	continue;
								const Tensor<Tdata> &D_a1b0 = tools.get_Ds_ab(Label::ab::a1b0, Aa01, Ab01);
								if(D_a1b0.empty())	continue;

								// a0a2b0 = a1a0a2 * a1b0
								const Tensor<Tdata> D_tmp1 = Tensor_Multiply::x1x2y1_ax1x2_ay1(D_a, D_a1b0);
								// a2b0b1 = a0a2b0 * a0b1
								Tensor<Tdata> D_tmp2 = Tensor_Multiply::x1x2y1_ax1x2_ay1(D_tmp1, D_a0b1);
								LRI_Cal_Aux::add_Ds(std::move(D_tmp2), D_mul);
							}
							if(D_mul.empty())	continue;

							// D_result = D_mul * D_b
							for(const TAC &Ab2 : list_Ab2)
							{
								if(this->filter_atom->filter_for32(label,Aa2,Ab01,Ab2))	continue;
								const Tensor<Tdata> &D_b = tools.get_Ds_ab(Label::ab::b, Ab01, Ab2);
								if(D_b.empty())	continue;
								// a2b2 = a2b0b1 * b0b1b2
								Tensor<Tdata> D_tmp3 = Tensor_Multiply::x0y2_x0ab_aby2(D_mul, D_b);
								LRI_Cal_Aux::add_Ds(std::move(D_tmp3), Ds_result_fixed[Ab2]);
							}
						} // end for Ab01

						if(!Ds_result_fixed.empty())
							LRI_Cal_Aux::add_Ds(LRI_Cal_Aux::Ds_translate(std::move(Ds_result_fixed), Aa2.second, this->period),
												Ds_result_thread[Aa2.first]);
						LRI_Cal_Aux::add_Ds_omp_try_map(Ds_result_thread, Ds_result, lock_Ds_result_add_map, fac_add_Ds);
					} // end for Ab2
				} break; // end case a0b1_a1b0

			  // Aab_Aab::a01b01_a01b2

				case Label::ab_ab::a0b0_a1b2:
				{
					const std::vector<TA >  list_Aa01 = LRI_Cal_Aux::filter_list_map( LRI_Cal_Aux::filter_list_map(
						list_Aa01_Da,
						data_wrapper(Label::ab::a0b0).Ds_ab ),
						data_wrapper(Label::ab::a1b2).Ds_ab );
					const std::vector<TAC> &list_Aa2 =
						list_Aa2_Da;
					const std::vector<TAC>  list_Ab01 = LRI_Cal_Aux::filter_list_set(
						list_Ab01_Db,
						data_wrapper(Label::ab::a0b0).index_Ds_ab[0]);
					const std::vector<TAC>  list_Ab2 = LRI_Cal_Aux::filter_list_set(
						list_Ab2_Db,
						data_wrapper(Label::ab::a1b2).index_Ds_ab[0]);

					for(const TAC &Ab01 : list_Ab01)
					{
						if(this->filter_atom->filter_for1(label,Ab01))	continue;
						std::map<TAC,Tensor<Tdata>> Ds_result_fixed;

						#pragma omp for schedule(dynamic) nowait
						for(std::size_t ia01=0; ia01<list_Aa01.size(); ++ia01)
						{
							const TA &Aa01 = list_Aa01[ia01];
							if(this->filter_atom->filter_for2(label,Ab01,Aa01))	continue;
							const Tensor<Tdata> &D_a0b0 = tools.get_Ds_ab(Label::ab::a0b0, Aa01, Ab01);
							if(D_a0b0.empty())	continue;
							// D_mul = D_b * D_a1b2
							Tensor<Tdata> D_mul;
							for(const TAC &Ab2 : list_Ab2)
							{
								if(this->filter_atom->filter_for31(label,Ab01,Aa01,Ab2))	continue;
								const Tensor<Tdata> &D_b = tools.get_Ds_ab(Label::ab::b, Ab01, Ab2);
								if(D_b.empty())	continue;
								const Tensor<Tdata> &D_a1b2 = tools.get_Ds_ab(Label::ab::a1b2, Aa01, Ab2);
								if(D_a1b2.empty())	continue;

								// b0b1a1 = b0b1b2 * a1b2
								Tensor<Tdata> D_tmp1 = Tensor_Multiply::x0x1y0_x0x1a_y0a(D_b, D_a1b2);
								LRI_Cal_Aux::add_Ds(std::move(D_tmp1), D_mul);
							}
							if(D_mul.empty())	continue;

							// D_result = D_mul * D_a * D_a0b0
							for(const TAC &Aa2 : list_Aa2)
							{
								if(this->filter_atom->filter_for32(label,Ab01,Aa01,Aa2))	continue;
								const Tensor<Tdata> &D_a = Global_Func::find(Ds_a_transpose, Aa01, Aa2);
								if(D_a.empty())	continue;
								// b1a1a0 = b0b1a1 * a0b0
								const Tensor<Tdata> D_tmp2 = Tensor_Multiply::x1x2y0_ax1x2_y0a(D_mul, D_a0b0);
								// a2b1 = a1a0a2 * b1a1a0
								Tensor<Tdata> D_tmp3 = Tensor_Multiply::x2y0_abx2_y0ab(D_a, D_tmp2);
								LRI_Cal_Aux::add_Ds(std::move(D_tmp3), Ds_result_fixed[Aa2]);
							}
						} // end for Aa01

						if(!Ds_result_fixed.empty())
							LRI_Cal_Aux::add_Ds(LRI_Cal_Aux::Ds_exchange(std::move(Ds_result_fixed), Ab01, this->period),
												Ds_result_thread);
						LRI_Cal_Aux::add_Ds_omp_try_map(Ds_result_thread, Ds_result, lock_Ds_result_add_map, fac_add_Ds);
					} // end for Ab01
				} break; // end case a0b0_a1b2

				case Label::ab_ab::a0b1_a1b2:
				{
					const std::vector<TA >  list_Aa01 = LRI_Cal_Aux::filter_list_map( LRI_Cal_Aux::filter_list_map(
						list_Aa01_Da,
						data_wrapper(Label::ab::a0b1).Ds_ab ),
						data_wrapper(Label::ab::a1b2).Ds_ab );
					const std::vector<TAC> &list_Aa2 =
						list_Aa2_Da;
					const std::vector<TAC>  list_Ab01 = LRI_Cal_Aux::filter_list_set(
						list_Ab01_Db,
						data_wrapper(Label::ab::a0b1).index_Ds_ab[0]);
					const std::vector<TAC>  list_Ab2 = LRI_Cal_Aux::filter_list_set(
						list_Ab2_Db,
						data_wrapper(Label::ab::a1b2).index_Ds_ab[0]);

					for(const TAC &Ab01 : list_Ab01)
					{
						if(this->filter_atom->filter_for1(label,Ab01))	continue;
						std::map<TAC,Tensor<Tdata>> Ds_result_fixed;

						#pragma omp for schedule(dynamic) nowait
						for(std::size_t ia01=0; ia01<list_Aa01.size(); ++ia01)
						{
							const TA &Aa01 = list_Aa01[ia01];
							if(this->filter_atom->filter_for2(label,Ab01,Aa01))	continue;
							const Tensor<Tdata> &D_a0b1 = tools.get_Ds_ab(Label::ab::a0b1, Aa01, Ab01);
							if(D_a0b1.empty())	continue;
							// D_mul = D_b * D_a1b2
							Tensor<Tdata> D_mul;
							for(const TAC &Ab2 : list_Ab2)
							{
								if(this->filter_atom->filter_for31(label,Ab01,Aa01,Ab2))	continue;
								const Tensor<Tdata> &D_b = tools.get_Ds_ab(Label::ab::b, Ab01, Ab2);
								if(D_b.empty())	continue;
								const Tensor<Tdata> &D_a1b2 = tools.get_Ds_ab(Label::ab::a1b2, Aa01, Ab2);
								if(D_a1b2.empty())	continue;

								// a1b0b1 = a1b2 * b0b1b2
								Tensor<Tdata> D_tmp1 = Tensor_Multiply::x0y0y1_x0a_y0y1a(D_a1b2, D_b);
								LRI_Cal_Aux::add_Ds(std::move(D_tmp1), D_mul);
							}
							if(D_mul.empty())	continue;

							// D_result = D_mul * D_a * D_a0b1
							for(const TAC &Aa2 : list_Aa2)
							{
								if(this->filter_atom->filter_for32(label,Ab01,Aa01,Aa2))	continue;
								const Tensor<Tdata> &D_a = tools.get_Ds_ab(Label::ab::a, Aa01, Aa2);
								if(D_a.empty())	continue;
								// a0a1b0 = a0b1 * a1b0b1
								const Tensor<Tdata> D_tmp2 = Tensor_Multiply::x0y0y1_x0a_y0y1a(D_a0b1, D_mul);
								// a2b0 = a0a1a2 * a0a1b0
								Tensor<Tdata> D_tmp3 = Tensor_Multiply::x2y2_abx2_aby2(D_a, D_tmp2);
								LRI_Cal_Aux::add_Ds(std::move(D_tmp3), Ds_result_fixed[Aa2]);
							}
						} // end for Aa01

						if(!Ds_result_fixed.empty())
							LRI_Cal_Aux::add_Ds(LRI_Cal_Aux::Ds_exchange(std::move(Ds_result_fixed), Ab01, this->period),
												Ds_result_thread);
						LRI_Cal_Aux::add_Ds_omp_try_map(Ds_result_thread, Ds_result, lock_Ds_result_add_map, fac_add_Ds);
					} // end for Ab01
				} break; // end case a0b1_a1b2

				case Label::ab_ab::a0b2_a1b0:
				{
					const std::vector<TA >  list_Aa01 = LRI_Cal_Aux::filter_list_map( LRI_Cal_Aux::filter_list_map(
						list_Aa01_Da,
						data_wrapper(Label::ab::a1b0).Ds_ab ),
						data_wrapper(Label::ab::a0b2).Ds_ab );
					const std::vector<TAC> &list_Aa2 =
						list_Aa2_Da;
					const std::vector<TAC>  list_Ab01 = LRI_Cal_Aux::filter_list_set(
						list_Ab01_Db,
						data_wrapper(Label::ab::a1b0).index_Ds_ab[0]);
					const std::vector<TAC>  list_Ab2 = LRI_Cal_Aux::filter_list_set(
						list_Ab2_Db,
						data_wrapper(Label::ab::a0b2).index_Ds_ab[0]);

					for(const TAC &Ab01 : list_Ab01)
					{
						if(this->filter_atom->filter_for1(label,Ab01))	continue;
						std::map<TAC,Tensor<Tdata>> Ds_result_fixed;

						#pragma omp for schedule(dynamic) nowait
						for(std::size_t ia01=0; ia01<list_Aa01.size(); ++ia01)
						{
							const TA &Aa01 = list_Aa01[ia01];
							if(this->filter_atom->filter_for2(label,Ab01,Aa01))	continue;
							const Tensor<Tdata> &D_a1b0 = tools.get_Ds_ab(Label::ab::a1b0, Aa01, Ab01);
							if(D_a1b0.empty())	continue;
							// D_mul = D_b * D_a0b2
							Tensor<Tdata> D_mul;
							for(const TAC &Ab2 : list_Ab2)
							{
								if(this->filter_atom->filter_for31(label,Ab01,Aa01,Ab2))	continue;
								const Tensor<Tdata> &D_b = tools.get_Ds_ab(Label::ab::b, Ab01, Ab2);
								if(D_b.empty())	continue;
								const Tensor<Tdata> &D_a0b2 = tools.get_Ds_ab(Label::ab::a0b2, Aa01, Ab2);
								if(D_a0b2.empty())	continue;

								// b0b1a0 = b0b1b2 * a0b2
								Tensor<Tdata> D_tmp1 = Tensor_Multiply::x0x1y0_x0x1a_y0a(D_b, D_a0b2);
								LRI_Cal_Aux::add_Ds(std::move(D_tmp1), D_mul);
							}
							if(D_mul.empty())	continue;

							// D_result = D_mul * D_a * D_a1b0
							for(const TAC &Aa2 : list_Aa2)
							{
								if(this->filter_atom->filter_for32(label,Ab01,Aa01,Aa2))	continue;
								const Tensor<Tdata> &D_a = tools.get_Ds_ab(Label::ab::a, Aa01, Aa2);
								if(D_a.empty())	continue;
								// b1a0a1 = b0b1a0 * a1b0
								const Tensor<Tdata> D_tmp2 = Tensor_Multiply::x1x2y0_ax1x2_y0a(D_mul, D_a1b0);
								// a2b1 = a0a1a2 * b1a0a1
								Tensor<Tdata> D_tmp3 = Tensor_Multiply::x2y0_abx2_y0ab(D_a, D_tmp2);
								LRI_Cal_Aux::add_Ds(std::move(D_tmp3), Ds_result_fixed[Aa2]);
							}
						} // end for Aa01

						if(!Ds_result_fixed.empty())
							LRI_Cal_Aux::add_Ds(LRI_Cal_Aux::Ds_exchange(std::move(Ds_result_fixed), Ab01, this->period),
												Ds_result_thread);
						LRI_Cal_Aux::add_Ds_omp_try_map(Ds_result_thread, Ds_result, lock_Ds_result_add_map, fac_add_Ds);
					} // end for Ab01
				} break; // end case a0b2_a1b0

				case Label::ab_ab::a0b2_a1b1:
				{
					const std::vector<TA >  list_Aa01 = LRI_Cal_Aux::filter_list_map( LRI_Cal_Aux::filter_list_map(
						list_Aa01_Da,
						data_wrapper(Label::ab::a1b1).Ds_ab ),
						data_wrapper(Label::ab::a0b2).Ds_ab );
					const std::vector<TAC> &list_Aa2 =
						list_Aa2_Da;
					const std::vector<TAC>  list_Ab01 = LRI_Cal_Aux::filter_list_set(
						list_Ab01_Db,
						data_wrapper(Label::ab::a1b1).index_Ds_ab[0]);
					const std::vector<TAC>  list_Ab2 = LRI_Cal_Aux::filter_list_set(
						list_Ab2_Db,
						data_wrapper(Label::ab::a0b2).index_Ds_ab[0]);

					for(const TAC &Ab01 : list_Ab01)
					{
						if(this->filter_atom->filter_for1(label,Ab01))	continue;
						std::map<TAC,Tensor<Tdata>> Ds_result_fixed;

						#pragma omp for schedule(dynamic) nowait
						for(std::size_t ia01=0; ia01<list_Aa01.size(); ++ia01)
						{
							const TA &Aa01 = list_Aa01[ia01];
							if(this->filter_atom->filter_for2(label,Ab01,Aa01))	continue;
							const Tensor<Tdata> &D_a1b1 = tools.get_Ds_ab(Label::ab::a1b1, Aa01, Ab01);
							if(D_a1b1.empty())	continue;
							// D_mul = D_b * D_a0b2
							Tensor<Tdata> D_mul;
							for(const TAC &Ab2 : list_Ab2)
							{
								if(this->filter_atom->filter_for31(label,Ab01,Aa01,Ab2))	continue;
								const Tensor<Tdata> &D_b = tools.get_Ds_ab(Label::ab::b, Ab01, Ab2);
								if(D_b.empty())	continue;
								const Tensor<Tdata> &D_a0b2 = tools.get_Ds_ab(Label::ab::a0b2, Aa01, Ab2);
								if(D_a0b2.empty())	continue;

								// a0b0b1 = a0b2 * b0b1b2
								Tensor<Tdata> D_tmp1 = Tensor_Multiply::x0y0y1_x0a_y0y1a(D_a0b2, D_b);
								LRI_Cal_Aux::add_Ds(std::move(D_tmp1), D_mul);
							}
							if(D_mul.empty())	continue;

							// D_result = D_mul * D_a * D_a1b1
							for(const TAC &Aa2 : list_Aa2)
							{
								if(this->filter_atom->filter_for32(label,Ab01,Aa01,Aa2))	continue;
								const Tensor<Tdata> &D_a = Global_Func::find(Ds_a_transpose, Aa01, Aa2);
								if(D_a.empty())	continue;
								// a1a0b0 = a1b1 * a0b0b1
								const Tensor<Tdata> D_tmp2 = Tensor_Multiply::x0y0y1_x0a_y0y1a(D_a1b1, D_mul);
								// a2b0 = a1a0a2 * a1a0b0
								Tensor<Tdata> D_tmp3 = Tensor_Multiply::x2y2_abx2_aby2(D_a, D_tmp2);
								LRI_Cal_Aux::add_Ds(std::move(D_tmp3), Ds_result_fixed[Aa2]);
							}
						} // end for Aa01

						if(!Ds_result_fixed.empty())
							LRI_Cal_Aux::add_Ds(LRI_Cal_Aux::Ds_exchange(std::move(Ds_result_fixed), Ab01, this->period),
												Ds_result_thread);
						LRI_Cal_Aux::add_Ds_omp_try_map(Ds_result_thread, Ds_result, lock_Ds_result_add_map, fac_add_Ds);
					} // end for Ab01
				} break; // end case a0b2_a1b1

			  // Aab_Aab::a01b01_a2b01

				case Label::ab_ab::a0b0_a2b1:
				{
					const std::vector<TA >  list_Aa01 = LRI_Cal_Aux::filter_list_map(
						list_Aa01_Da,
						data_wrapper(Label::ab::a0b0).Ds_ab );
					const std::vector<TAC>  list_Aa2 = LRI_Cal_Aux::filter_list_map(
						list_Aa2_Da,
						data_wrapper(Label::ab::a2b1).Ds_ab );
					const std::vector<TAC>  list_Ab01 = LRI_Cal_Aux::filter_list_set( LRI_Cal_Aux::filter_list_set(
						list_Ab01_Db,
						data_wrapper(Label::ab::a0b0).index_Ds_ab[0]),
						data_wrapper(Label::ab::a2b1).index_Ds_ab[0]);
					const std::vector<TAC> &list_Ab2 =
						list_Ab2_Db;

					for(const TA &Aa01 : list_Aa01)
					{
						if(this->filter_atom->filter_for1(label,Aa01))	continue;
						std::map<TAC,Tensor<Tdata>> Ds_result_fixed;

						#pragma omp for schedule(dynamic) nowait
						for(std::size_t ib01=0; ib01<list_Ab01.size(); ++ib01)
						{
							const TAC &Ab01 = list_Ab01[ib01];
							if(this->filter_atom->filter_for2(label,Aa01,Ab01))	continue;
							const Tensor<Tdata> &D_a0b0 = tools.get_Ds_ab(Label::ab::a0b0, Aa01, Ab01);
							if(D_a0b0.empty())	continue;
							// D_mul = D_a * D_a2b1
							Tensor<Tdata> D_mul;
							for(const TAC &Aa2 : list_Aa2)
							{
								if(this->filter_atom->filter_for31(label,Aa01,Ab01,Aa2))	continue;
								const Tensor<Tdata> &D_a = Global_Func::find(Ds_a_transpose, Aa01, Aa2);
								if(D_a.empty())	continue;
								const Tensor<Tdata> &D_a2b1 = tools.get_Ds_ab(Label::ab::a2b1, Aa2, Ab01);
								if(D_a2b1.empty())	continue;

								// b1a1a0 = a2b1 * a1a0a2
								Tensor<Tdata> D_tmp1 = Tensor_Multiply::x1y0y1_ax1_y0y1a(D_a2b1, D_a);
								LRI_Cal_Aux::add_Ds(std::move(D_tmp1), D_mul);
							}
							if(D_mul.empty())	continue;

							// D_result = D_mul * D_a0b0 * D_b
							for(const TAC &Ab2 : list_Ab2)
							{
								if(this->filter_atom->filter_for32(label,Aa01,Ab01,Ab2))	continue;
								const Tensor<Tdata> &D_b = tools.get_Ds_ab(Label::ab::b, Ab01, Ab2);
								if(D_b.empty())	continue;

								// b0b1a1 = a0b0 * b1a1a0
								const Tensor<Tdata> D_tmp2 = Tensor_Multiply::x1y0y1_ax1_y0y1a(D_a0b0, D_mul);
								// a1b2 = b0b1a1 * b0b1b2
								Tensor<Tdata> D_tmp3 = Tensor_Multiply::x2y2_abx2_aby2(D_tmp2, D_b);
								LRI_Cal_Aux::add_Ds(std::move(D_tmp3), Ds_result_fixed[Ab2]);
							}
						} // end for Ab01

						if(!Ds_result_fixed.empty())
							LRI_Cal_Aux::add_Ds(std::move(Ds_result_fixed),
												Ds_result_thread[Aa01]);
						LRI_Cal_Aux::add_Ds_omp_try_map(Ds_result_thread, Ds_result, lock_Ds_result_add_map, fac_add_Ds);
					} // end for Aa01
				} break; // end case a0b0_a2b1

				case Label::ab_ab::a0b1_a2b0:
				{
					const std::vector<TA >  list_Aa01 = LRI_Cal_Aux::filter_list_map(
						list_Aa01_Da,
						data_wrapper(Label::ab::a0b1).Ds_ab );
					const std::vector<TAC>  list_Aa2 = LRI_Cal_Aux::filter_list_map(
						list_Aa2_Da,
						data_wrapper(Label::ab::a2b0).Ds_ab );
					const std::vector<TAC>  list_Ab01 = LRI_Cal_Aux::filter_list_set( LRI_Cal_Aux::filter_list_set(
						list_Ab01_Db,
						data_wrapper(Label::ab::a0b1).index_Ds_ab[0]),
						data_wrapper(Label::ab::a2b0).index_Ds_ab[0]);
					const std::vector<TAC> &list_Ab2 =
						list_Ab2_Db;

					for(const TA &Aa01 : list_Aa01)
					{
						if(this->filter_atom->filter_for1(label,Aa01))	continue;
						std::map<TAC,Tensor<Tdata>> Ds_result_fixed;

						#pragma omp for schedule(dynamic) nowait
						for(std::size_t ib01=0; ib01<list_Ab01.size(); ++ib01)
						{
							const TAC &Ab01 = list_Ab01[ib01];
							if(this->filter_atom->filter_for2(label,Aa01,Ab01))	continue;
							const Tensor<Tdata> &D_a0b1 = tools.get_Ds_ab(Label::ab::a0b1, Aa01, Ab01);
							if(D_a0b1.empty())	continue;
							// D_mul = D_a * D_a2b0
							Tensor<Tdata> D_mul;
							for(const TAC &Aa2 : list_Aa2)
							{
								if(this->filter_atom->filter_for31(label,Aa01,Ab01,Aa2))	continue;
								const Tensor<Tdata> &D_a = tools.get_Ds_ab(Label::ab::a, Aa01, Aa2);
								if(D_a.empty())	continue;
								const Tensor<Tdata> &D_a2b0 = tools.get_Ds_ab(Label::ab::a2b0, Aa2, Ab01);
								if(D_a2b0.empty())	continue;

								// a0a1b0 = a0a1a2 * a2b0
								Tensor<Tdata> D_tmp1 = Tensor_Multiply::x0x1y1_x0x1a_ay1(D_a, D_a2b0);
								LRI_Cal_Aux::add_Ds(std::move(D_tmp1), D_mul);
							}
							if(D_mul.empty())	continue;

							// D_result = D_mul * D_a0b1 * D_b
							for(const TAC &Ab2 : list_Ab2)
							{
								if(this->filter_atom->filter_for32(label,Aa01,Ab01,Ab2))	continue;
								const Tensor<Tdata> &D_b = tools.get_Ds_ab(Label::ab::b, Ab01, Ab2);
								if(D_b.empty())	continue;

								// a1b0b1 = a0a1b0 * a0b1
								const Tensor<Tdata> D_tmp2 = Tensor_Multiply::x1x2y1_ax1x2_ay1(D_mul, D_a0b1);
								// a1b2 = a1b0b1 * b0b1b2
								Tensor<Tdata> D_tmp3 = Tensor_Multiply::x0y2_x0ab_aby2(D_tmp2, D_b);
								LRI_Cal_Aux::add_Ds(std::move(D_tmp3), Ds_result_fixed[Ab2]);
							}
						} // end for Ab01

						if(!Ds_result_fixed.empty())
							LRI_Cal_Aux::add_Ds(std::move(Ds_result_fixed),
												Ds_result_thread[Aa01]);
						LRI_Cal_Aux::add_Ds_omp_try_map(Ds_result_thread, Ds_result, lock_Ds_result_add_map, fac_add_Ds);
					} // end for Aa01
				} break; // end case a0b1_a2b0

				case Label::ab_ab::a1b0_a2b1:
				{
					const std::vector<TA >  list_Aa01 = LRI_Cal_Aux::filter_list_map(
						list_Aa01_Da,
						data_wrapper(Label::ab::a1b0).Ds_ab );
					const std::vector<TAC>  list_Aa2 = LRI_Cal_Aux::filter_list_map(
						list_Aa2_Da,
						data_wrapper(Label::ab::a2b1).Ds_ab );
					const std::vector<TAC>  list_Ab01 = LRI_Cal_Aux::filter_list_set( LRI_Cal_Aux::filter_list_set(
						list_Ab01_Db,
						data_wrapper(Label::ab::a1b0).index_Ds_ab[0]),
						data_wrapper(Label::ab::a2b1).index_Ds_ab[0]);
					const std::vector<TAC> &list_Ab2 =
						list_Ab2_Db;

					for(const TA &Aa01 : list_Aa01)
					{
						if(this->filter_atom->filter_for1(label,Aa01))	continue;
						std::map<TAC,Tensor<Tdata>> Ds_result_fixed;

						#pragma omp for schedule(dynamic) nowait
						for(std::size_t ib01=0; ib01<list_Ab01.size(); ++ib01)
						{
							const TAC &Ab01 = list_Ab01[ib01];
							if(this->filter_atom->filter_for2(label,Aa01,Ab01))	continue;
							const Tensor<Tdata> &D_a1b0 = tools.get_Ds_ab(Label::ab::a1b0, Aa01, Ab01);
							if(D_a1b0.empty())	continue;
							// D_mul = D_a * D_a2b1
							Tensor<Tdata> D_mul;
							for(const TAC &Aa2 : list_Aa2)
							{
								if(this->filter_atom->filter_for31(label,Aa01,Ab01,Aa2))	continue;
								const Tensor<Tdata> &D_a = tools.get_Ds_ab(Label::ab::a, Aa01, Aa2);
								if(D_a.empty())	continue;
								const Tensor<Tdata> &D_a2b1 = tools.get_Ds_ab(Label::ab::a2b1, Aa2, Ab01);
								if(D_a2b1.empty())	continue;

								// b1a0a1 = a2b1 * a0a1a2
								Tensor<Tdata> D_tmp1 = Tensor_Multiply::x1y0y1_ax1_y0y1a(D_a2b1, D_a);
								LRI_Cal_Aux::add_Ds(std::move(D_tmp1), D_mul);
							}
							if(D_mul.empty())	continue;

							// D_result = D_mul * D_a1b0 * D_b
							for(const TAC &Ab2 : list_Ab2)
							{
								if(this->filter_atom->filter_for32(label,Aa01,Ab01,Ab2))	continue;
								const Tensor<Tdata> &D_b = tools.get_Ds_ab(Label::ab::b, Ab01, Ab2);
								if(D_b.empty())	continue;

								// b0b1a0 = a1b0 * b1a0a1
								const Tensor<Tdata> D_tmp2 = Tensor_Multiply::x1y0y1_ax1_y0y1a(D_a1b0, D_mul);
								// a0b2 = b0b1a0 * b0b1b2
								Tensor<Tdata> D_tmp3 = Tensor_Multiply::x2y2_abx2_aby2(D_tmp2, D_b);
								LRI_Cal_Aux::add_Ds(std::move(D_tmp3), Ds_result_fixed[Ab2]);
							}
						} // end for Ab01

						if(!Ds_result_fixed.empty())
							LRI_Cal_Aux::add_Ds(std::move(Ds_result_fixed),
												Ds_result_thread[Aa01]);
						LRI_Cal_Aux::add_Ds_omp_try_map(Ds_result_thread, Ds_result, lock_Ds_result_add_map, fac_add_Ds);
					} // end for Aa01
				} break; // end case a1b0_a2b1

				case Label::ab_ab::a1b1_a2b0:
				{
					const std::vector<TA >  list_Aa01 = LRI_Cal_Aux::filter_list_map(
						list_Aa01_Da,
						data_wrapper(Label::ab::a1b1).Ds_ab );
					const std::vector<TAC>  list_Aa2 = LRI_Cal_Aux::filter_list_map(
						list_Aa2_Da,
						data_wrapper(Label::ab::a2b0).Ds_ab );
					const std::vector<TAC>  list_Ab01 = LRI_Cal_Aux::filter_list_set( LRI_Cal_Aux::filter_list_set(
						list_Ab01_Db,
						data_wrapper(Label::ab::a1b1).index_Ds_ab[0]),
						data_wrapper(Label::ab::a2b0).index_Ds_ab[0]);
					const std::vector<TAC> &list_Ab2 =
						list_Ab2_Db;

					for(const TA &Aa01 : list_Aa01)
					{
						if(this->filter_atom->filter_for1(label,Aa01))	continue;
						std::map<TAC,Tensor<Tdata>> Ds_result_fixed;

						#pragma omp for schedule(dynamic) nowait
						for(std::size_t ib01=0; ib01<list_Ab01.size(); ++ib01)
						{
							const TAC &Ab01 = list_Ab01[ib01];
							if(this->filter_atom->filter_for2(label,Aa01,Ab01))	continue;
							const Tensor<Tdata> &D_a1b1 = tools.get_Ds_ab(Label::ab::a1b1, Aa01, Ab01);
							if(D_a1b1.empty())	continue;
							// D_mul = D_a * D_a2b0
							Tensor<Tdata> D_mul;
							for(const TAC &Aa2 : list_Aa2)
							{
								if(this->filter_atom->filter_for31(label,Aa01,Ab01,Aa2))	continue;
								const Tensor<Tdata> &D_a = Global_Func::find(Ds_a_transpose, Aa01, Aa2);
								if(D_a.empty())	continue;
								const Tensor<Tdata> &D_a2b0 = tools.get_Ds_ab(Label::ab::a2b0, Aa2, Ab01);
								if(D_a2b0.empty())	continue;

								// a1a0b0 = a1a0a2 * a2b0
								Tensor<Tdata> D_tmp1 = Tensor_Multiply::x0x1y1_x0x1a_ay1(D_a, D_a2b0);
								LRI_Cal_Aux::add_Ds(std::move(D_tmp1), D_mul);
							}
							if(D_mul.empty())	continue;

							// D_result = D_mul * D_a1b1 * D_b
							for(const TAC &Ab2 : list_Ab2)
							{
								if(this->filter_atom->filter_for32(label,Aa01,Ab01,Ab2))	continue;
								const Tensor<Tdata> &D_b = tools.get_Ds_ab(Label::ab::b, Ab01, Ab2);
								if(D_b.empty())	continue;

								// a0b0b1 = a1a0b0 * a1b1
								const Tensor<Tdata> D_tmp2 = Tensor_Multiply::x1x2y1_ax1x2_ay1(D_mul, D_a1b1);
								// a0b2 = a0b0b1 * b0b1b2
								Tensor<Tdata> D_tmp3 = Tensor_Multiply::x0y2_x0ab_aby2(D_tmp2, D_b);
								LRI_Cal_Aux::add_Ds(std::move(D_tmp3), Ds_result_fixed[Ab2]);
							}
						} // end for Ab01

						if(!Ds_result_fixed.empty())
							LRI_Cal_Aux::add_Ds(std::move(Ds_result_fixed),
												Ds_result_thread[Aa01]);
						LRI_Cal_Aux::add_Ds_omp_try_map(Ds_result_thread, Ds_result, lock_Ds_result_add_map, fac_add_Ds);
					} // end for Aa01
				} break; // end case a1b1_a2b0

			  // Aab_Aab::a01b01_a2b2

				case Label::ab_ab::a0b0_a2b2:
				{
					const std::vector<TA >  list_Aa01 = LRI_Cal_Aux::filter_list_map(
						list_Aa01_Da,
						data_wrapper(Label::ab::a0b0).Ds_ab );
					const std::vector<TAC>  list_Aa2 = LRI_Cal_Aux::filter_list_map(
						list_Aa2_Da,
						data_wrapper(Label::ab::a2b2).Ds_ab );
					const std::vector<TAC>  list_Ab01 = LRI_Cal_Aux::filter_list_set(
						list_Ab01_Db,
						data_wrapper(Label::ab::a0b0).index_Ds_ab[0]);
					const std::vector<TAC>  list_Ab2 = LRI_Cal_Aux::filter_list_set(
						list_Ab2_Db,
						data_wrapper(Label::ab::a2b2).index_Ds_ab[0]);

					for(const TA &Aa01 : list_Aa01)
					{
						if(this->filter_atom->filter_for1(label,Aa01))	continue;
						std::map<TAC,Tensor<Tdata>> Ds_result_fixed;

						#pragma omp for schedule(dynamic) nowait
						for(std::size_t ib2=0; ib2<list_Ab2.size(); ++ib2)
						{
							const TAC &Ab2 = list_Ab2[ib2];
							if(this->filter_atom->filter_for2(label,Aa01,Ab2))	continue;
							// D_mul = D_a * D_a2b2
							Tensor<Tdata> D_mul;
							for(const TAC &Aa2 : list_Aa2)
							{
								if(this->filter_atom->filter_for31(label,Aa01,Ab2,Aa2))	continue;
								const Tensor<Tdata> &D_a = Global_Func::find(Ds_a_transpose, Aa01, Aa2);
								if(D_a.empty())	continue;
								const Tensor<Tdata> &D_a2b2 = tools.get_Ds_ab(Label::ab::a2b2, Aa2, Ab2);
								if(D_a2b2.empty())	continue;

								// b2a1a0 = a2b2 * a1a0a2
								Tensor<Tdata> D_tmp1 = Tensor_Multiply::x1y0y1_ax1_y0y1a(D_a2b2, D_a);
								LRI_Cal_Aux::add_Ds(std::move(D_tmp1), D_mul);
							}
							if(D_mul.empty())	continue;

							// D_result = D_mul * D_a0b0 * D_b
							for(const TAC &Ab01 : list_Ab01)
							{
								if(this->filter_atom->filter_for32(label,Aa01,Ab2,Ab01))	continue;
								const Tensor<Tdata> &D_b = Global_Func::find(Ds_b_transpose, Ab01.first, TAC{Ab2.first, (Ab2.second-Ab01.second)%this->period});
								if(D_b.empty())	continue;
								const Tensor<Tdata> &D_a0b0 = tools.get_Ds_ab(Label::ab::a0b0, Aa01, Ab01);
								if(D_a0b0.empty())	continue;

								// b0b2a1 = a0b0 * b2a1a0
								const Tensor<Tdata> D_tmp2 = Tensor_Multiply::x1y0y1_ax1_y0y1a(D_a0b0, D_mul);
								// a1b1 = b0b2a1 * b1b0b2
								Tensor<Tdata> D_tmp3 = Tensor_Multiply::x2y0_abx2_y0ab(D_tmp2, D_b);
								LRI_Cal_Aux::add_Ds(std::move(D_tmp3), Ds_result_fixed[Ab01]);
							}
						} // end for Ab01

						if(!Ds_result_fixed.empty())
							LRI_Cal_Aux::add_Ds(std::move(Ds_result_fixed),
												Ds_result_thread[Aa01]);
						LRI_Cal_Aux::add_Ds_omp_try_map(Ds_result_thread, Ds_result, lock_Ds_result_add_map, fac_add_Ds);
					} // end for Aa01
				} break; // end case a0b0_a2b2

				case Label::ab_ab::a0b1_a2b2:
				{
					const std::vector<TA >  list_Aa01 = LRI_Cal_Aux::filter_list_map(
						list_Aa01_Da,
						data_wrapper(Label::ab::a0b1).Ds_ab );
					const std::vector<TAC>  list_Aa2 = LRI_Cal_Aux::filter_list_map(
						list_Aa2_Da,
						data_wrapper(Label::ab::a2b2).Ds_ab );
					const std::vector<TAC>  list_Ab01 = LRI_Cal_Aux::filter_list_set(
						list_Ab01_Db,
						data_wrapper(Label::ab::a0b1).index_Ds_ab[0]);
					const std::vector<TAC>  list_Ab2 = LRI_Cal_Aux::filter_list_set(
						list_Ab2_Db,
						data_wrapper(Label::ab::a2b2).index_Ds_ab[0]);

					for(const TA &Aa01 : list_Aa01)
					{
						if(this->filter_atom->filter_for1(label,Aa01))	continue;
						std::map<TAC,Tensor<Tdata>> Ds_result_fixed;

						#pragma omp for schedule(dynamic) nowait
						for(std::size_t ib2=0; ib2<list_Ab2.size(); ++ib2)
						{
							const TAC &Ab2 = list_Ab2[ib2];
							if(this->filter_atom->filter_for2(label,Aa01,Ab2))	continue;
							// D_mul = D_a * D_a2b2
							Tensor<Tdata> D_mul;
							for(const TAC &Aa2 : list_Aa2)
							{
								if(this->filter_atom->filter_for31(label,Aa01,Ab2,Aa2))	continue;
								const Tensor<Tdata> &D_a = Global_Func::find(Ds_a_transpose, Aa01, Aa2);
								if(D_a.empty())	continue;
								const Tensor<Tdata> &D_a2b2 = tools.get_Ds_ab(Label::ab::a2b2, Aa2, Ab2);
								if(D_a2b2.empty())	continue;

								// b2a1a0 = a2b2 * a1a0a2
								Tensor<Tdata> D_tmp1 = Tensor_Multiply::x1y0y1_ax1_y0y1a(D_a2b2, D_a);
								LRI_Cal_Aux::add_Ds(std::move(D_tmp1), D_mul);
							}
							if(D_mul.empty())	continue;

							// D_result = D_mul * D_a0b1 * D_b
							for(const TAC &Ab01 : list_Ab01)
							{
								if(this->filter_atom->filter_for32(label,Aa01,Ab2,Ab01))	continue;
								const Tensor<Tdata> &D_b = tools.get_Ds_ab(Label::ab::b, Ab01, Ab2);
								if(D_b.empty())	continue;
								const Tensor<Tdata> &D_a0b1 = tools.get_Ds_ab(Label::ab::a0b1, Aa01, Ab01);
								if(D_a0b1.empty())	continue;

								// b1b2a1 = a0b1 * a2a1a0
								const Tensor<Tdata> D_tmp2 = Tensor_Multiply::x1y0y1_ax1_y0y1a(D_a0b1, D_mul);
								// a1b0 = b1b2a1 * b0b1b2
								Tensor<Tdata> D_tmp3 = Tensor_Multiply::x2y0_abx2_y0ab(D_tmp2, D_b);
								LRI_Cal_Aux::add_Ds(std::move(D_tmp3), Ds_result_fixed[Ab01]);
							}
						} // end for Ab01

						if(!Ds_result_fixed.empty())
							LRI_Cal_Aux::add_Ds(std::move(Ds_result_fixed),
												Ds_result_thread[Aa01]);
						LRI_Cal_Aux::add_Ds_omp_try_map(Ds_result_thread, Ds_result, lock_Ds_result_add_map, fac_add_Ds);
					} // end for Aa01
				} break; // end case a0b1_a2b2

				case Label::ab_ab::a1b0_a2b2:
				{
					const std::vector<TA >  list_Aa01 = LRI_Cal_Aux::filter_list_map(
						list_Aa01_Da,
						data_wrapper(Label::ab::a1b0).Ds_ab );
					const std::vector<TAC>  list_Aa2 = LRI_Cal_Aux::filter_list_map(
						list_Aa2_Da,
						data_wrapper(Label::ab::a2b2).Ds_ab );
					const std::vector<TAC>  list_Ab01 = LRI_Cal_Aux::filter_list_set(
						list_Ab01_Db,
						data_wrapper(Label::ab::a1b0).index_Ds_ab[0]);
					const std::vector<TAC>  list_Ab2 = LRI_Cal_Aux::filter_list_set(
						list_Ab2_Db,
						data_wrapper(Label::ab::a2b2).index_Ds_ab[0]);

					for(const TA &Aa01 : list_Aa01)
					{
						if(this->filter_atom->filter_for1(label,Aa01))	continue;
						std::map<TAC,Tensor<Tdata>> Ds_result_fixed;

						#pragma omp for schedule(dynamic) nowait
						for(std::size_t ib2=0; ib2<list_Ab2.size(); ++ib2)
						{
							const TAC &Ab2 = list_Ab2[ib2];
							if(this->filter_atom->filter_for2(label,Aa01,Ab2))	continue;
							// D_mul = D_a * D_a2b2
							Tensor<Tdata> D_mul;
							for(const TAC &Aa2 : list_Aa2)
							{
								if(this->filter_atom->filter_for31(label,Aa01,Ab2,Aa2))	continue;
								const Tensor<Tdata> &D_a = tools.get_Ds_ab(Label::ab::a, Aa01, Aa2);
								if(D_a.empty())	continue;
								const Tensor<Tdata> &D_a2b2 = tools.get_Ds_ab(Label::ab::a2b2, Aa2, Ab2);
								if(D_a2b2.empty())	continue;

								// b2a0a1 = a2b2 * a0a1a2
								Tensor<Tdata> D_tmp1 = Tensor_Multiply::x1y0y1_ax1_y0y1a(D_a2b2, D_a);
								LRI_Cal_Aux::add_Ds(std::move(D_tmp1), D_mul);
							}
							if(D_mul.empty())	continue;

							// D_result = D_mul * D_a1b0 * D_b
							for(const TAC &Ab01 : list_Ab01)
							{
								if(this->filter_atom->filter_for32(label,Aa01,Ab2,Ab01))	continue;
								const Tensor<Tdata> &D_b = Global_Func::find(Ds_b_transpose, Ab01.first, TAC{Ab2.first, (Ab2.second-Ab01.second)%this->period});
								if(D_b.empty())	continue;
								const Tensor<Tdata> &D_a1b0 = tools.get_Ds_ab(Label::ab::a1b0, Aa01, Ab01);
								if(D_a1b0.empty())	continue;

								// b0b2a0 = a1b0 * b2a0a1
								const Tensor<Tdata> D_tmp2 = Tensor_Multiply::x1y0y1_ax1_y0y1a(D_a1b0, D_mul);
								// a0b1 = b0b2a0 * b1b0b2
								Tensor<Tdata> D_tmp3 = Tensor_Multiply::x2y0_abx2_y0ab(D_tmp2, D_b);
								LRI_Cal_Aux::add_Ds(std::move(D_tmp3), Ds_result_fixed[Ab01]);
							}
						} // end for Ab01

						if(!Ds_result_fixed.empty())
							LRI_Cal_Aux::add_Ds(std::move(Ds_result_fixed),
												Ds_result_thread[Aa01]);
						LRI_Cal_Aux::add_Ds_omp_try_map(Ds_result_thread, Ds_result, lock_Ds_result_add_map, fac_add_Ds);
					} // end for Aa01
				} break; // end case a1b0_a2b2

				case Label::ab_ab::a1b1_a2b2:
				{
					const std::vector<TA >  list_Aa01 = LRI_Cal_Aux::filter_list_map(
						list_Aa01_Da,
						data_wrapper(Label::ab::a1b1).Ds_ab );
					const std::vector<TAC>  list_Aa2 = LRI_Cal_Aux::filter_list_map(
						list_Aa2_Da,
						data_wrapper(Label::ab::a2b2).Ds_ab );
					const std::vector<TAC>  list_Ab01 = LRI_Cal_Aux::filter_list_set(
						list_Ab01_Db,
						data_wrapper(Label::ab::a1b1).index_Ds_ab[0]);
					const std::vector<TAC>  list_Ab2 = LRI_Cal_Aux::filter_list_set(
						list_Ab2_Db,
						data_wrapper(Label::ab::a2b2).index_Ds_ab[0]);

					for(const TA &Aa01 : list_Aa01)
					{
						if(this->filter_atom->filter_for1(label,Aa01))	continue;
						std::map<TAC,Tensor<Tdata>> Ds_result_fixed;

						#pragma omp for schedule(dynamic) nowait
						for(std::size_t ib2=0; ib2<list_Ab2.size(); ++ib2)
						{
							const TAC &Ab2 = list_Ab2[ib2];
							if(this->filter_atom->filter_for2(label,Aa01,Ab2))	continue;
							// D_mul = D_a * D_a2b2
							Tensor<Tdata> D_mul;
							for(const TAC &Aa2 : list_Aa2)
							{
								if(this->filter_atom->filter_for31(label,Aa01,Ab2,Aa2))	continue;
								const Tensor<Tdata> &D_a = tools.get_Ds_ab(Label::ab::a, Aa01, Aa2);
								if(D_a.empty())	continue;
								const Tensor<Tdata> &D_a2b2 = tools.get_Ds_ab(Label::ab::a2b2, Aa2, Ab2);
								if(D_a2b2.empty())	continue;

								// b2a0a1 = a2b2 * a0a1a2
								Tensor<Tdata> D_tmp1 = Tensor_Multiply::x1y0y1_ax1_y0y1a(D_a2b2, D_a);
								LRI_Cal_Aux::add_Ds(std::move(D_tmp1), D_mul);
							}
							if(D_mul.empty())	continue;

							// D_result = D_mul * D_a1b1 * D_b
							for(const TAC &Ab01 : list_Ab01)
							{
								if(this->filter_atom->filter_for32(label,Aa01,Ab2,Ab01))	continue;
								const Tensor<Tdata> &D_b = tools.get_Ds_ab(Label::ab::b, Ab01, Ab2);
								if(D_b.empty())	continue;
								const Tensor<Tdata> &D_a1b1 = tools.get_Ds_ab(Label::ab::a1b1, Aa01, Ab01);
								if(D_a1b1.empty())	continue;

								// b1b2a0 = a1b1 * b2a0a1
								const Tensor<Tdata> D_tmp2 = Tensor_Multiply::x1y0y1_ax1_y0y1a(D_a1b1, D_mul);
								// a0b0 = b1b2a0 * b0b1b2
								Tensor<Tdata> D_tmp3 = Tensor_Multiply::x2y0_abx2_y0ab(D_tmp2, D_b);
								LRI_Cal_Aux::add_Ds(std::move(D_tmp3), Ds_result_fixed[Ab01]);
							}
						} // end for Ab01

						if(!Ds_result_fixed.empty())
							LRI_Cal_Aux::add_Ds(std::move(Ds_result_fixed),
												Ds_result_thread[Aa01]);
						LRI_Cal_Aux::add_Ds_omp_try_map(Ds_result_thread, Ds_result, lock_Ds_result_add_map, fac_add_Ds);
					} // end for Aa01
				} break; // end case a1b1_a2b2

			  // Aab_Aab::a01b2_a2b01

				case Label::ab_ab::a1b2_a2b1:
				{
					const std::vector<TA >  list_Aa01 = LRI_Cal_Aux::filter_list_map(
						list_Aa01_Da,
						data_wrapper(Label::ab::a1b2).Ds_ab );
					const std::vector<TAC> &list_Aa2 = LRI_Cal_Aux::filter_list_map(
						list_Aa2_Da,
						data_wrapper(Label::ab::a2b1).Ds_ab );
					const std::vector<TAC>  list_Ab01 = LRI_Cal_Aux::filter_list_set(
						list_Ab01_Db,
						data_wrapper(Label::ab::a2b1).index_Ds_ab[0]);
					const std::vector<TAC>  list_Ab2 = LRI_Cal_Aux::filter_list_set(
						list_Ab2_Db,
						data_wrapper(Label::ab::a1b2).index_Ds_ab[0]);

					for(const TA &Aa01 : list_Aa01)
					{
						if(this->filter_atom->filter_for1(label,Aa01))	continue;
						#pragma omp for schedule(dynamic) nowait
						for(std::size_t ib01=0; ib01<list_Ab01.size(); ++ib01)
						{
							const TAC &Ab01 = list_Ab01[ib01];
							if(this->filter_atom->filter_for2(label,Aa01,Ab01))	continue;
							// D_mul1 = D_b * D_a1b2
							Tensor<Tdata> D_mul1;
							for(const TAC &Ab2 : list_Ab2)
							{
								if(this->filter_atom->filter_for31(label,Aa01,Ab01,Ab2))	continue;
								const Tensor<Tdata> &D_b = tools.get_Ds_ab(Label::ab::b, Ab01, Ab2);
								if(D_b.empty())	continue;
								const Tensor<Tdata> &D_a1b2 = tools.get_Ds_ab(Label::ab::a1b2, Aa01, Ab2);
								if(D_a1b2.empty())	continue;

								// b0b1a1 = b0b1b2 * a1b2
								Tensor<Tdata> D_tmp1 = Tensor_Multiply::x0x1y0_x0x1a_y0a(D_b, D_a1b2);
								LRI_Cal_Aux::add_Ds(std::move(D_tmp1), D_mul1);
							}
							if(D_mul1.empty())	continue;

							// D_mul2 = D_a2b1 * D_a
							Tensor<Tdata> D_mul2;
							for(const TAC &Aa2 : list_Aa2)
							{
								if(this->filter_atom->filter_for32(label,Aa01,Ab01,Aa2))	continue;
								const Tensor<Tdata> &D_a = Global_Func::find(Ds_a_transpose, Aa01, Aa2);
								if(D_a.empty())	continue;
								const Tensor<Tdata> &D_a2b1 = tools.get_Ds_ab(Label::ab::a2b1, Aa2, Ab01);
								if(D_a2b1.empty())	continue;
								// b1a1a0 = a2b1 * a1a0a2
								Tensor<Tdata> D_tmp2 = Tensor_Multiply::x1y0y1_ax1_y0y1a(D_a2b1, D_a);
								LRI_Cal_Aux::add_Ds(std::move(D_tmp2), D_mul2);
							}
							if(D_mul2.empty())	continue;

							// D_result = D_mul2 * D_mul1
							// a0b0 = b1a1a0 * b0b1a1
							Tensor<Tdata> D_mul3 = Tensor_Multiply::x2y0_abx2_y0ab(D_mul2, D_mul1);
							LRI_Cal_Aux::add_Ds(std::move(D_mul3),
												Ds_result_thread[Aa01][Ab01]);
						} // end for Aa01

						LRI_Cal_Aux::add_Ds_omp_try_map(Ds_result_thread, Ds_result, lock_Ds_result_add_map, fac_add_Ds);
					} // end for Ab01
				} break; // end case a1b2_a2b1

				case Label::ab_ab::a0b2_a2b0:
				{
					const std::vector<TA >  list_Aa01 = LRI_Cal_Aux::filter_list_map(
						list_Aa01_Da,
						data_wrapper(Label::ab::a0b2).Ds_ab );
					const std::vector<TAC> &list_Aa2 = LRI_Cal_Aux::filter_list_map(
						list_Aa2_Da,
						data_wrapper(Label::ab::a2b0).Ds_ab );
					const std::vector<TAC>  list_Ab01 = LRI_Cal_Aux::filter_list_set(
						list_Ab01_Db,
						data_wrapper(Label::ab::a2b0).index_Ds_ab[0]);
					const std::vector<TAC>  list_Ab2 = LRI_Cal_Aux::filter_list_set(
						list_Ab2_Db,
						data_wrapper(Label::ab::a0b2).index_Ds_ab[0]);

					for(const TA &Aa01 : list_Aa01)
					{
						if(this->filter_atom->filter_for1(label,Aa01))	continue;
						#pragma omp for schedule(dynamic) nowait
						for(std::size_t ib01=0; ib01<list_Ab01.size(); ++ib01)
						{
							const TAC &Ab01 = list_Ab01[ib01];
							if(this->filter_atom->filter_for2(label,Aa01,Ab01))	continue;
							// D_mul1 = D_b * D_a0b2
							Tensor<Tdata> D_mul1;
							for(const TAC &Ab2 : list_Ab2)
							{
								if(this->filter_atom->filter_for31(label,Aa01,Ab01,Ab2))	continue;
								const Tensor<Tdata> &D_b = tools.get_Ds_ab(Label::ab::b, Ab01, Ab2);
								if(D_b.empty())	continue;
								const Tensor<Tdata> &D_a0b2 = tools.get_Ds_ab(Label::ab::a0b2, Aa01, Ab2);
								if(D_a0b2.empty())	continue;

								// a0b0b1 = a0b2 * b0b1b2
								Tensor<Tdata> D_tmp1 = Tensor_Multiply::x0y0y1_x0a_y0y1a(D_a0b2, D_b);
								LRI_Cal_Aux::add_Ds(std::move(D_tmp1), D_mul1);
							}
							if(D_mul1.empty())	continue;

							// D_mul2 = D_a2b0 * D_a
							Tensor<Tdata> D_mul2;
							for(const TAC &Aa2 : list_Aa2)
							{
								if(this->filter_atom->filter_for32(label,Aa01,Ab01,Aa2))	continue;
								const Tensor<Tdata> &D_a = Global_Func::find(Ds_a_transpose, Aa01, Aa2);
								if(D_a.empty())	continue;
								const Tensor<Tdata> &D_a2b0 = tools.get_Ds_ab(Label::ab::a2b0, Aa2, Ab01);
								if(D_a2b0.empty())	continue;
								// a1a0b0 = a1a0a2 * a2b0
								Tensor<Tdata> D_tmp2 = Tensor_Multiply::x0x1y1_x0x1a_ay1(D_a, D_a2b0);
								LRI_Cal_Aux::add_Ds(std::move(D_tmp2), D_mul2);
							}
							if(D_mul2.empty())	continue;

							// D_result = D_mul2 * D_mul1
							// b1a1 = a1a0b0 * a0b0b1
							Tensor<Tdata> D_mul3 = Tensor_Multiply::x0y2_x0ab_aby2(D_mul2, D_mul1);
							LRI_Cal_Aux::add_Ds(std::move(D_mul3),
												Ds_result_thread[Aa01][Ab01]);
						} // end for Aa01

						LRI_Cal_Aux::add_Ds_omp_try_map(Ds_result_thread, Ds_result, lock_Ds_result_add_map, fac_add_Ds);
					} // end for Ab01
				} break; // end case a0b2_a2b0

				case Label::ab_ab::a0b2_a2b1:
				{
					const std::vector<TA >  list_Aa01 = LRI_Cal_Aux::filter_list_map(
						list_Aa01_Da,
						data_wrapper(Label::ab::a0b2).Ds_ab );
					const std::vector<TAC> &list_Aa2 = LRI_Cal_Aux::filter_list_map(
						list_Aa2_Da,
						data_wrapper(Label::ab::a2b1).Ds_ab );
					const std::vector<TAC>  list_Ab01 = LRI_Cal_Aux::filter_list_set(
						list_Ab01_Db,
						data_wrapper(Label::ab::a2b1).index_Ds_ab[0]);
					const std::vector<TAC>  list_Ab2 = LRI_Cal_Aux::filter_list_set(
						list_Ab2_Db,
						data_wrapper(Label::ab::a0b2).index_Ds_ab[0]);

					for(const TA &Aa01 : list_Aa01)
					{
						if(this->filter_atom->filter_for1(label,Aa01))	continue;
						#pragma omp for schedule(dynamic) nowait
						for(std::size_t ib01=0; ib01<list_Ab01.size(); ++ib01)
						{
							const TAC &Ab01 = list_Ab01[ib01];
							if(this->filter_atom->filter_for2(label,Aa01,Ab01))	continue;
							// D_mul1 = D_b * D_a0b2
							Tensor<Tdata> D_mul1;
							for(const TAC &Ab2 : list_Ab2)
							{
								if(this->filter_atom->filter_for31(label,Aa01,Ab01,Ab2))	continue;
								const Tensor<Tdata> &D_b = tools.get_Ds_ab(Label::ab::b, Ab01, Ab2);
								if(D_b.empty())	continue;
								const Tensor<Tdata> &D_a0b2 = tools.get_Ds_ab(Label::ab::a0b2, Aa01, Ab2);
								if(D_a0b2.empty())	continue;

								// b0b1a0 = b0b1b2 * a0b2
								Tensor<Tdata> D_tmp1 = Tensor_Multiply::x0x1y0_x0x1a_y0a(D_b, D_a0b2);
								LRI_Cal_Aux::add_Ds(std::move(D_tmp1), D_mul1);
							}
							if(D_mul1.empty())	continue;

							// D_mul2 = D_a2b1 * D_a
							Tensor<Tdata> D_mul2;
							for(const TAC &Aa2 : list_Aa2)
							{
								if(this->filter_atom->filter_for32(label,Aa01,Ab01,Aa2))	continue;
								const Tensor<Tdata> &D_a = tools.get_Ds_ab(Label::ab::a, Aa01, Aa2);
								if(D_a.empty())	continue;
								const Tensor<Tdata> &D_a2b1 = tools.get_Ds_ab(Label::ab::a2b1, Aa2, Ab01);
								if(D_a2b1.empty())	continue;
								// b1a0a1 = a2b1 * a0a1a2
								Tensor<Tdata> D_tmp2 = Tensor_Multiply::x1y0y1_ax1_y0y1a(D_a2b1, D_a);
								LRI_Cal_Aux::add_Ds(std::move(D_tmp2), D_mul2);
							}
							if(D_mul2.empty())	continue;

							// D_result = D_mul2 * D_mul1
							// a1b0 = b1a0a1 * b0b1a0
							Tensor<Tdata> D_mul3 = Tensor_Multiply::x2y0_abx2_y0ab(D_mul2, D_mul1);
							LRI_Cal_Aux::add_Ds(std::move(D_mul3),
												Ds_result_thread[Aa01][Ab01]);
						} // end for Aa01

						LRI_Cal_Aux::add_Ds_omp_try_map(Ds_result_thread, Ds_result, lock_Ds_result_add_map, fac_add_Ds);
					} // end for Ab01
				} break; // end case a0b2_a2b1

				case Label::ab_ab::a1b2_a2b0:
				{
					const std::vector<TA >  list_Aa01 = LRI_Cal_Aux::filter_list_map(
						list_Aa01_Da,
						data_wrapper(Label::ab::a1b2).Ds_ab );
					const std::vector<TAC> &list_Aa2 = LRI_Cal_Aux::filter_list_map(
						list_Aa2_Da,
						data_wrapper(Label::ab::a2b0).Ds_ab );
					const std::vector<TAC>  list_Ab01 = LRI_Cal_Aux::filter_list_set(
						list_Ab01_Db,
						data_wrapper(Label::ab::a2b0).index_Ds_ab[0]);
					const std::vector<TAC>  list_Ab2 = LRI_Cal_Aux::filter_list_set(
						list_Ab2_Db,
						data_wrapper(Label::ab::a1b2).index_Ds_ab[0]);

					for(const TA &Aa01 : list_Aa01)
					{
						if(this->filter_atom->filter_for1(label,Aa01))	continue;
						#pragma omp for schedule(dynamic) nowait
						for(std::size_t ib01=0; ib01<list_Ab01.size(); ++ib01)
						{
							const TAC &Ab01 = list_Ab01[ib01];
							if(this->filter_atom->filter_for2(label,Aa01,Ab01))	continue;
							// D_mul1 = D_b * D_a1b2
							Tensor<Tdata> D_mul1;
							for(const TAC &Ab2 : list_Ab2)
							{
								if(this->filter_atom->filter_for31(label,Aa01,Ab01,Ab2))	continue;
								const Tensor<Tdata> &D_b = tools.get_Ds_ab(Label::ab::b, Ab01, Ab2);
								if(D_b.empty())	continue;
								const Tensor<Tdata> &D_a1b2 = tools.get_Ds_ab(Label::ab::a1b2, Aa01, Ab2);
								if(D_a1b2.empty())	continue;

								// a1b0b1 = a1b2 * b0b1b2
								Tensor<Tdata> D_tmp1 = Tensor_Multiply::x0y0y1_x0a_y0y1a(D_a1b2, D_b);
								LRI_Cal_Aux::add_Ds(std::move(D_tmp1), D_mul1);
							}
							if(D_mul1.empty())	continue;

							// D_mul2 = D_a2b0 * D_a
							Tensor<Tdata> D_mul2;
							for(const TAC &Aa2 : list_Aa2)
							{
								if(this->filter_atom->filter_for32(label,Aa01,Ab01,Aa2))	continue;
								const Tensor<Tdata> &D_a = tools.get_Ds_ab(Label::ab::a, Aa01, Aa2);
								if(D_a.empty())	continue;
								const Tensor<Tdata> &D_a2b0 = tools.get_Ds_ab(Label::ab::a2b0, Aa2, Ab01);
								if(D_a2b0.empty())	continue;
								// a0a1b0 = a0a1a2 * a2b0
								Tensor<Tdata> D_tmp2 = Tensor_Multiply::x0x1y1_x0x1a_ay1(D_a, D_a2b0);
								LRI_Cal_Aux::add_Ds(std::move(D_tmp2), D_mul2);
							}
							if(D_mul2.empty())	continue;

							// D_result = D_mul2 * D_mul1
							// a0b1 = a0a1b0 * a1b0b1
							Tensor<Tdata> D_mul3 = Tensor_Multiply::x0y2_x0ab_aby2(D_mul2, D_mul1);
							LRI_Cal_Aux::add_Ds(std::move(D_mul3),
												Ds_result_thread[Aa01][Ab01]);
						} // end for Aa01

						LRI_Cal_Aux::add_Ds_omp_try_map(Ds_result_thread, Ds_result, lock_Ds_result_add_map, fac_add_Ds);
					} // end for Ab01
				} break; // end case a1b2_a2b0

				default:
					throw std::invalid_argument(std::string(__FILE__)+std::to_string(__LINE__));
			} // end switch(label)
		} // end for label

		LRI_Cal_Aux::add_Ds_omp_wait_map(Ds_result_thread, Ds_result, lock_Ds_result_add_map, fac_add_Ds);
	} // end #pragma omp parallel

	LRI_Cal_Aux::destroy_lock_result(lock_Ds_result_add_map, Ds_result);

  #ifdef __MKL_RI
	mkl_set_num_threads(mkl_threads);
  #endif
}	// end LRI::cal_loop3()

}	// end namespace RI

