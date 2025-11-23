#include "utils_cmake.h"

#include "envs_cmake.h"
#include "utils_io.h"

#include <cstring>

namespace LIBRPA
{

namespace utils
{

void print_cmake_info()
{
    lib_printf("LibRPA CMake build info:\n");
    lib_printf("| Build type               : %s\n", envs::cmake_build_type);
    lib_printf("| Source code directory    : %s\n", envs::source_dir);
    lib_printf("| Source code Git reference: %s\n", envs::git_ref);
    lib_printf("| Source code Git hash     : %s\n", envs::git_hash);
    lib_printf("| C++ compiler             : %s\n", envs::cxx_compiler);
    lib_printf("| C++ compiler flags       : %s\n", envs::cxx_compiler_flags);
    lib_printf("| Fortran compiler         : %s\n", envs::fortran_compiler);
    lib_printf("| Fortran compiler flags   : %s\n", envs::fortran_compiler_flags);
    lib_printf("| Use GreenX API           : %s\n", envs::use_greenx_api);
    lib_printf("| Use LibRI                : %s\n", envs::use_libri);
    lib_printf("| LibRI include directory  : %s\n", envs::libri_include_dir);
    lib_printf("| LibComm include directory: %s\n", envs::libcomm_include_dir);

}

} /* end of namespace utils */

} /* end of namespace LIBRPA */
