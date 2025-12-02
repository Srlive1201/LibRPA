#include "utils_cmake.h"

#include "envs_cmake.h"
#include "../io/global_io.h"

#include <sstream>

namespace librpa_int
{

void print_cmake_info()
{
    // global::lib_printf("LibRPA CMake build info:\n");
    // global::lib_printf("| Build type               : %s\n", envs::cmake_build_type);
    // global::lib_printf("| Source code directory    : %s\n", envs::source_dir);
    // global::lib_printf("| Source code Git reference: %s\n", envs::git_ref);
    // global::lib_printf("| Source code Git hash     : %s\n", envs::git_hash);
    // global::lib_printf("| C++ compiler             : %s\n", envs::cxx_compiler);
    // global::lib_printf("| C++ compiler flags       : %s\n", envs::cxx_compiler_flags);
    // global::lib_printf("| Fortran compiler         : %s\n", envs::fortran_compiler);
    // global::lib_printf("| Fortran compiler flags   : %s\n", envs::fortran_compiler_flags);
    // global::lib_printf("| Use GreenX API           : %s\n", envs::use_greenx_api);
    // global::lib_printf("| Use LibRI                : %s\n", envs::use_libri);
    // global::lib_printf("| LibRI include directory  : %s\n", envs::libri_include_dir);
    // global::lib_printf("| LibComm include directory: %s\n", envs::libcomm_include_dir);
    global::lib_printf(cmake_info_storage().c_str());
}

const std::string& cmake_info_storage()
{
    static const std::string s = []{
        std::stringstream ss;
        ss << "LibRPA CMake build info:" << "\n"
        << "| Build type               : " << envs::cmake_build_type << "\n"
        << "| Source code directory    : " << envs::source_dir << "\n"
        << "| Source code Git reference: " << envs::git_ref << "\n"
        << "| Source code Git hash     : " << envs::git_hash << "\n"
        << "| C++ compiler             : " << envs::cxx_compiler << "\n"
        << "| C++ compiler flags       : " << envs::cxx_compiler_flags << "\n"
        << "| Fortran compiler         : " << envs::fortran_compiler << "\n"
        << "| Fortran compiler flags   : " << envs::fortran_compiler_flags << "\n"
        << "| Use GreenX API           : " << envs::use_greenx_api << "\n"
        << "| Use LibRI                : " << envs::use_libri << "\n"
        << "| LibRI include directory  : " << envs::libri_include_dir << "\n"
        << "| LibComm include directory: " << envs::libcomm_include_dir << "\n";
        return ss.str();
    }();
    return s;
}

} /* end of namespace librpa_int */
