// ===================
//  Author: Peize Lin
//  date: 2022.09.20
// ===================

#pragma once

#include "unittests/global/Cereal_Test.hpp"
#include "unittests/Comm_Keys/Comm_Keys_3-test.hpp"
#include "unittests/Comm_Trans/Communicate_Map-test-1.hpp"
#include "unittests/Comm_Assemble/Communicate_Map-test-2.hpp"
#include "unittests/Comm_Assemble/Communicate_Map-test-speed.hpp"

namespace Test_All
{
	inline void test_all(int argc, char *argv[])
	{
		Cereal_Test::main(argc,argv);
		Comm_Keys_3_Test::test_all(argc,argv);
		Communicate_Map_Test::test_transmission(argc,argv);
		Communicate_Map_Test::test_assemble(argc,argv);
		Communicate_Map_Test::test_speed(argc,argv);
	}
}