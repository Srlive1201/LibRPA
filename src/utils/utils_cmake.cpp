#include "utils_cmake.h"

#include "envs_cmake.h"
#include "../io/global_io.h"

#include <cstring>

namespace librpa_int
{

void print_cmake_info()
{
    global::lib_printf("LibRPA CMake build info:\n");
    global::lib_printf("| Build type               : %s\n", envs::cmake_build_type);
    global::lib_printf("| Source code directory    : %s\n", envs::source_dir);
    global::lib_printf("| Source code Git reference: %s\n", envs::git_ref);
    global::lib_printf("| Source code Git hash     : %s\n", envs::git_hash);
    global::lib_printf("| C++ compiler             : %s\n", envs::cxx_compiler);
    global::lib_printf("| C++ compiler flags       : %s\n", envs::cxx_compiler_flags);
    global::lib_printf("| Fortran compiler         : %s\n", envs::fortran_compiler);
    global::lib_printf("| Fortran compiler flags   : %s\n", envs::fortran_compiler_flags);
    global::lib_printf("| Use GreenX API           : %s\n", envs::use_greenx_api);
    global::lib_printf("| Use LibRI                : %s\n", envs::use_libri);
    global::lib_printf("| LibRI include directory  : %s\n", envs::libri_include_dir);
    global::lib_printf("| LibComm include directory: %s\n", envs::libcomm_include_dir);

}

} /* end of namespace librpa_int */
