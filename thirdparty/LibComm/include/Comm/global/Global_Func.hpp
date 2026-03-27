//=======================
// AUTHOR : Peize Lin
// DATE :   2023-02-15
//=======================

#pragma once

#include "Global_Func.h"

#include <fstream>
#include <string>
#include <stdexcept>

// BSD like macOS does not have /proc/meminfo
// include unistd.h to get the number of system memory.
#if defined(__MACH__)
#include <unistd.h>
#endif

namespace Comm
{

namespace Global_Func
{
#if defined(__MACH__)
	inline std::size_t memory_available()
	{
		// HACK: always assume half of the system memory free. Only for build and test on macOS
		const std::size_t pages = sysconf(_SC_PHYS_PAGES);
		const std::size_t page_size = sysconf(_SC_PAGE_SIZE);
		return pages * page_size / 2;
	}
#else
	inline std::size_t memory_available()
	{
		constexpr std::size_t kB_to_B = 1024;
		std::ifstream ifs("/proc/meminfo");
		int num = 0;
		std::size_t mem_sum = 0;
		while (ifs.good())
		{
			std::string label, size, kB;
			ifs >> label >> size >> kB;
			if (label == "MemAvailable:")
			{
				return std::stol(size) * kB_to_B;
			}
			else if (label == "MemFree:" || label == "Buffers:" || label == "Cached:")
			{
				mem_sum += std::stol(size);
				++num;
			}

			if(num==3)
			{
				return mem_sum * kB_to_B;
			}
		}
		throw std::runtime_error("read /proc/meminfo error in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));
	}
#endif
}

}