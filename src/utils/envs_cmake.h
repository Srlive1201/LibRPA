#pragma once

namespace LIBRPA
{

namespace envs
{

//! The path of the source code directory
extern const char * source_dir;

//! Git hash of the source code
extern const char * git_hash;

//! Git reference of the source code
extern const char * git_ref;

//! CMake C++ compiler
extern const char * cxx_compiler;

//! CMake Fortran compiler
extern const char * fortran_compiler;

//! CMake C++ compiler flags
extern const char * cxx_compiler_flags;

//! CMake Fortran compiler flags
extern const char * fortran_compiler_flags;

// CMake Options
//! CMake build type
extern const char * cmake_build_type;

extern const char * use_libri;

extern const char * libri_include_dir;

extern const char * libcomm_include_dir;

extern const char * use_greenx_api;

} /* end of namespace envs */
} /* end of namespace LIBRPA */
