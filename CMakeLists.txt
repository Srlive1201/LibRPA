# main CMake file for LibRPA
cmake_minimum_required(VERSION 3.13)

# extract version number
file(READ "src/Makefile" mkfile)
string(REGEX MATCH " *VERSION *= *[0-9.]*" version ${mkfile})
string(REGEX REPLACE " *VERSION *= *" "" version ${version})
# Project information
project(LibRPA VERSION ${version} LANGUAGES NONE)

# options setup
option(USE_LIBRI "Use LibRI for tensor contraction" OFF)
option(USE_CMAKE_INC "Use cmake.in for configure" ON)
option(DEBUG "Switch on debug mode and facility" OFF)
option(BUILD_LIBRARY "Whether to build library" OFF)

# Add cmake load files
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
include(UseDoxygenDoc)

# C++ standard requirement
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(math_libs "")
set(parallel_libs "")

if(EXISTS "${PROJECT_SOURCE_DIR}/cmake.inc" AND USE_CMAKE_INC)
  include("${PROJECT_SOURCE_DIR}/cmake.inc")
  list(APPEND math_libs ${LAPACK_LIBRARIES})
else()
  if(DEFINED ENV{MKLROOT} AND NOT DEFINED MKLROOT)
    set(MKLROOT "$ENV{MKLROOT}")
  endif()
  if(MKLROOT)
    find_package(IntelMKL REQUIRED)
    include_directories(${MKL_INCLUDE_DIRS} ${MKL_INCLUDE_DIRS}/fftw)
    list(APPEND math_libs IntelMKL::MKL)
  else()
    find_package(LAPACK REQUIRED)
    find_package(ScaLAPACK REQUIRED)
    list(APPEND math_libs
            ${LAPACK_LIBRARIES}
            ${ScaLAPACK_LIBRARIES}
            )
  endif()

  find_package(MPI REQUIRED)
  if (MPI_CXX_FOUND)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_FLAGS}")
  endif()

  find_package(OpenMP REQUIRED)
  if (OpenMP_CXX_FOUND)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  endif()
  list(APPEND parallel_libs
          OpenMP::OpenMP_CXX
          MPI::MPI_CXX
          )
  include_directories(${MPI_CXX_INCLUDE_PATH})
endif()

enable_language(CXX)

include_directories(src)

# use LibRI, which in turn depends on Cereal
if(USE_LIBRI)
  find_package(Cereal)
  if(NOT Cereal_FOUND)
    message(FATAL_ERROR "USE_LIBRI set, but no Cereal headers found as dependency")
  endif()
  find_package(LibComm)
  if(NOT LibComm_FOUND)
    message(FATAL_ERROR "USE_LIBRI set, but no LibComm headers found")
  endif()
  find_package(LibRI)
  if(NOT LibRI_FOUND)
    message(FATAL_ERROR "USE_LIBRI set, but no LibRI headers found")
  endif()
  add_compile_definitions("__USE_LIBRI")
  message(STATUS "The project will be compiled with LibRI")
endif()

if(DEBUG)
  message(STATUS "Debug mode switched on")
  add_compile_definitions("__DEBUG")
endif()

add_subdirectory(src)

add_subdirectory(docs)

enable_testing()
add_subdirectory(tests)
