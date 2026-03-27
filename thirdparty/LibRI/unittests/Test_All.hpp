// ===================
//  Author: Peize Lin
//  date: 2022.09.20
// ===================

#pragma once

#include "unittests/comm/mix/Communicate_Tensors-test.hpp"
#include "unittests/distribute/Distribute_Equally-test.hpp"
#include "unittests/distribute/Divide_Atoms-test.hpp"
#include "unittests/distribute/Split_Processes-test.hpp"
#include "unittests/global/Blas-test.hpp"
#include "unittests/global/Tensor-test-2.hpp"
#include "unittests/global/Tensor-test-3.hpp"
#include "unittests/global/Tensor_Multiply-test.hpp"
#include "unittests/global/Map_Operator-test.hpp"
//#include "unittests/ri/LRI-loop4-test.hpp"
#include "unittests/ri/LRI-loop3-test.hpp"
#include "unittests/ri/LRI-speed-test.hpp"
#include "unittests/ri/Cell_Nearest-test.hpp"
#include "unittests/physics/Exx-test.hpp"
#include "unittests/physics/RPA-test.hpp"
#include "unittests/physics/GW-test.hpp"

namespace Test_All
{
	static void test_all(int argc, char *argv[])
	{
		Communicate_Tensors_Test::test_comm_judge_map(argc, argv);
		Communicate_Tensors_Test::test_comm_judge_map2(argc, argv);
		Communicate_Tensors_Test::test_comm_judge_map3(argc, argv);
		Communicate_Tensors_Test::test_comm_judge_map2_first(argc, argv);
		Communicate_Tensors_Test::test_comm_judge_map3_first(argc, argv);
		Communicate_Tensors_Test::test_comm_judge_map2_period(argc, argv);
		Communicate_Tensors_Test::test_comm_judge_map3_period(argc, argv);

		Distribute_Equally_Test::test_distribute_atoms(argc, argv);
		Distribute_Equally_Test::test_distribute_atoms_periods(argc, argv);
		Distribute_Equally_Test::test_distribute_atoms_repeatable(argc, argv);
		Distribute_Equally_Test::test_distribute_periods(argc, argv);

		Divide_Atoms_Test::test_divide_atoms();
		Divide_Atoms_Test::test_divide_atoms_with_period();
		Divide_Atoms_Test::test_divide_atoms_periods();

		Split_Processes_Test::test_split_all(argc, argv);

		Blas_Test::test_all();

		Tensor_Test::test_multiply_2();
		Tensor_Test::test_operator_all_3();

		Tensor_Multiply_Test::main<double>();
		Tensor_Multiply_Test::main<std::complex<double>>();

		Map_Operator_Test::test_union_map1();
		Map_Operator_Test::test_union_map2();
		Map_Operator_Test::test_intersection_map1();
		Map_Operator_Test::test_intersection_map2();
		Map_Operator_Test::test_intersection_map3();
		Map_Operator_Test::test_transform_map();
		Map_Operator_Test::test_for_each_map();

//		LRI_Loop4_Test::main<float>(argc, argv);
//		LRI_Loop4_Test::main<double>(argc, argv);
//		LRI_Loop4_Test::main<std::complex<float>>(argc, argv);
//		LRI_Loop4_Test::main<std::complex<double>>(argc, argv);

		LRI_Loop3_Test::main<float>(argc, argv);
		LRI_Loop3_Test::main<double>(argc, argv);
		LRI_Loop3_Test::main<std::complex<float>>(argc, argv);
		LRI_Loop3_Test::main<std::complex<double>>(argc, argv);

		LRI_Speed_Test::test_speed<float>(argc, argv, 1, 1);
		LRI_Speed_Test::test_speed<double>(argc, argv, 1, 1);
		LRI_Speed_Test::test_speed<std::complex<float>>(argc, argv, 1, 1);
		LRI_Speed_Test::test_speed<std::complex<double>>(argc, argv, 1, 1);

		Cell_Nearest_Test::main();

		Exx_Test::main<float>(argc, argv);
		Exx_Test::main<double>(argc, argv);
		Exx_Test::main<std::complex<float>>(argc, argv);
		Exx_Test::main<std::complex<double>>(argc, argv);

		RPA_Test::main<float>(argc, argv);
		RPA_Test::main<double>(argc, argv);
		RPA_Test::main<std::complex<float>>(argc, argv);
		RPA_Test::main<std::complex<double>>(argc, argv);

		GW_Test::main<float>(argc, argv);
		GW_Test::main<double>(argc, argv);
		GW_Test::main<std::complex<float>>(argc, argv);
		GW_Test::main<std::complex<double>>(argc, argv);
	}
}