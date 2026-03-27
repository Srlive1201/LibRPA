// ===================
//  Author: Peize Lin
//  date: 2022.08.12
// ===================

#pragma once

#include "LRI.h"
#include "LRI_Cal_Aux.h"
#include "../global/Array_Operator.h"

#include <omp.h>
#ifdef __MKL_RI
#include <mkl_service.h>
#endif

namespace RI
{

template<typename TA, typename Tcell, std::size_t Ndim, typename Tdata>
void LRI<TA,Tcell,Ndim,Tdata>::cal_loop4(
	const std::vector<Label::ab_ab> &labels,
	std::vector<std::map<TA, std::map<TAC, Tensor<Tdata>>>> &Ds_result)
{
	using namespace Array_Operator;

	const bool flag_D_b_transpose = [&labels]() -> bool
	{
		for(const Label::ab_ab &label : labels)
			if(LRI_Cal_Aux::judge_x(label)==0)
				return true;
		return false;
	}();
	const std::map<TA, std::map<TAC, Tensor<Tdata>>> Ds_b_transpose
		= flag_D_b_transpose
		? LRI_Cal_Aux::cal_Ds_transpose( this->data_pool.at( this->data_ab_name[Label::ab::b] ).Ds_ab )
		: std::map<TA, std::map<TAC, Tensor<Tdata>>>{};

	assert(!this->coefficients.empty());
	if(Ds_result.empty())
		Ds_result.resize(this->coefficients.size());
	else
		assert(Ds_result.size()==this->coefficients.size());

	omp_lock_t lock_Ds_result_add;
	omp_init_lock(&lock_Ds_result_add);

#ifdef __MKL_RI
	const std::size_t mkl_threads = mkl_get_max_threads();
//	if(!omp_get_nested())
//		mkl_set_num_threads(std::max(1UL,mkl_threads/list_Aa01.size()));
//	else
		mkl_set_num_threads(1);
#endif

	#pragma omp parallel
	{
		std::vector<std::map<TA, std::map<TAC, Tensor<Tdata>>>> Ds_result_thread(this->coefficients.size());
		LRI_Cal_Tools<TA,TC,Tdata> tools(this->period, this->data_pool, this->data_ab_name, Ds_result_thread);

		for(const TA &Aa01 : this->parallel->get_list_Aa01())
		{
			for(const TAC &Aa2 : this->parallel->get_list_Aa2(Aa01))
			{
				const Tensor<Tdata> D_a = tools.get_Ds_ab(Label::ab::a, Aa01, Aa2);
				if(D_a.empty())	continue;
				const Tensor<Tdata> D_a_transpose = LRI_Cal_Aux::tensor3_transpose(D_a);

				const std::vector<TAC> &list_Ab01 = this->parallel->get_list_Ab01(Aa01, Aa2);
				#pragma omp for schedule(dynamic) nowait
				for(std::size_t ib01=0; ib01<list_Ab01.size(); ++ib01)
				{
					const TAC &Ab01 = list_Ab01[ib01];
					const std::vector<Label::ab_ab> labels_filter_Ab01 = tools.filter_Ds_Ab01(labels, Aa01, Aa2, Ab01);
					if(labels_filter_Ab01.empty())	continue;

					std::unordered_map<Label::ab_ab, Tensor<Tdata>> Ds_b01;
					std::unordered_map<Label::ab_ab, Tdata_real> Ds_b01_csm;
					Ds_b01.reserve(Label::array_ab_ab.size());
					Ds_b01_csm.reserve(Label::array_ab_ab.size());

					for(const TAC &Ab2 : this->parallel->get_list_Ab2(Aa01, Aa2, Ab01))
					{
						const Tensor<Tdata> D_b = tools.get_Ds_ab(Label::ab::b, Ab01, Ab2);
						if(D_b.empty())	continue;
						const Tensor<Tdata> D_b_transpose
							= flag_D_b_transpose
							  ? Ds_b_transpose.at(Ab01.first).at({Ab2.first, (Ab2.second-Ab01.second)%this->period})
							  : Tensor<Tdata>{};

						for(const Label::ab_ab &label : labels_filter_Ab01)
						{
							this->cal_funcs[label](
								label,
								Aa01, Aa2, Ab01, Ab2,
								D_a, D_b,
								D_a_transpose, D_b_transpose,
								Ds_b01, Ds_b01_csm,
								tools);
						}
					} // end for Ab2

					if( !LRI_Cal_Aux::judge_Ds_empty(Ds_result_thread) && omp_test_lock(&lock_Ds_result_add) )
					{
						LRI_Cal_Aux::add_Ds(std::move(Ds_result_thread), Ds_result);
						omp_unset_lock(&lock_Ds_result_add);
						Ds_result_thread.clear();
						Ds_result_thread.resize(Ds_result.size());
					}
				} // end for Ab01
			}// end for Aa2
		}// end for Aa01

		if(!LRI_Cal_Aux::judge_Ds_empty(Ds_result_thread))
		{
			omp_set_lock(&lock_Ds_result_add);
			LRI_Cal_Aux::add_Ds(std::move(Ds_result_thread), Ds_result);
			omp_unset_lock(&lock_Ds_result_add);
			Ds_result_thread.clear();
			Ds_result_thread.resize(Ds_result.size());
		}
	} // end #pragma omp parallel

	omp_destroy_lock(&lock_Ds_result_add);
#ifdef __MKL_RI
	mkl_set_num_threads(mkl_threads);
#endif
}

}