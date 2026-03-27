// ===================
//  Author: Peize Lin
//  date: 2022.07.11
// ===================

#pragma once

#include "RI/distribute/Split_Processes.h"

#include <fstream>

namespace Split_Processes_Test
{
	static void test_split_all(int argc, char *argv[], const std::vector<std::size_t> &task_sizes)
	{
		MPI_Init(&argc, &argv);

		const std::vector<std::tuple<RI::MPI_Wrapper::mpi_comm, std::size_t, std::size_t>>
			comm_color_sizes = RI::Split_Processes::split_all(MPI_COMM_WORLD, task_sizes);

		std::ofstream ofs("out."+std::to_string(RI::MPI_Wrapper::mpi_get_rank(MPI_COMM_WORLD)));
		for(const auto &comm_color_size : comm_color_sizes)
			ofs<<RI::MPI_Wrapper::mpi_get_rank(std::get<0>(comm_color_size)())<<"\t"
			   <<RI::MPI_Wrapper::mpi_get_size(std::get<0>(comm_color_size)())<<"\t|\t"
			   <<std::get<1>(comm_color_size)<<"\t"
			   <<std::get<2>(comm_color_size)<<std::endl;

		MPI_Finalize();
	}

	static void test_split_all(int argc, char *argv[])
	{
		test_split_all(argc, argv, {8,12});
		/*
		mpirun -n 6								// num_average 恰好为4，完美均分
			rank 0
				0	6	|	0	1
				0	3	|	0	2
				0	1	|	0	3
			rank 1
				1	6	|	0	1
				1	3	|	0	2
				0	1	|	1	3
			rank 2
				2	6	|	0	1
				2	3	|	0	2
				0	1	|	2	3
			rank 3
				3	6	|	0	1
				0	3	|	1	2
				0	1	|	0	3
			rank 4
				4	6	|	0	1
				1	3	|	1	2
				0	1	|	1	3
			rank 5
				5	6	|	0	1
				2	3	|	1	2
				0	1	|	2	3
		*/


		test_split_all(argc, argv, {4,5,6});
		/*
		mpirun -n 7								// num_average 各维度、各进程不同，尽量均分
			rank 0
				0	7	|	0	1
				0	4	|	0	2
				0	2	|	0	2
				0	1	|	0	2
			rank 1
				1	7	|	0	1
				1	4	|	0	2
				1	2	|	0	2
				0	1	|	1	2
			rank 2
				2	7	|	0	1
				2	4	|	0	2
				0	2	|	1	2
				0	1	|	0	2
			rank 3
				3	7	|	0	1
				3	4	|	0	2
				1	2	|	1	2
				0	1	|	1	2
			rank 4
				4	7	|	0	1
				0	3	|	1	2
				0	2	|	0	2
				0	1	|	0	2
			rank 5
				5	7	|	0	1
				1	3	|	1	2
				1	2	|	0	2
				0	1	|	1	2
			rank 6
				6	7	|	0	1
				2	3	|	1	2
				0	1	|	1	2
				0	1	|	0	1
		*/

		test_split_all(argc, argv, {2,10});		// task_sizes[0] << task_sizes[1], 按照第1维划分
		/*
		mpirun -n 3
			rank 0
				0	3	|	0	1
				0	3	|	0	1
				0	1	|	0	3
			rank 1
				1	3	|	0	1
				1	3	|	0	1
				0	1	|	1	3
			rank 2
				2	3	|	0	1
				2	3	|	0	1
				0	1	|	2	3
		*/

		test_split_all(argc, argv, {2,3});
		/*
		mpirun -n 7								// task_sizes > rank_size, 01进程被分到同一task
			rank 0
				0	7	|	0	1
				0	4	|	0	2
				0	2	|	0	3
			rank 1
				1	7	|	0	1
				1	4	|	0	2
				1	2	|	0	3
			rank 2
				2	7	|	0	1
				2	4	|	0	2
				0	1	|	1	3
			rank 3
				3	7	|	0	1
				3	4	|	0	2
				0	1	|	2	3
			rank 4
				4	7	|	0	1
				0	3	|	1	2
				0	1	|	0	3
			rank 5
				5	7	|	0	1
				1	3	|	1	2
				0	1	|	1	3
			rank 6
				6	7	|	0	1
				2	3	|	1	2
				0	1	|	2	3
		*/
	}
}