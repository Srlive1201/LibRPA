# main CMake file for LibRPA
cmake_minimum_required(VERSION 3.13)

project(LibRPA NONE)

# extract version number
file(READ ${PROJECT_SOURCE_DIR}/src/version.h version_file)
string(REGEX MATCH "LIBRPA_MAJOR_VERSION ([0-9]+)" librpa_major_version ${version_file})
set(librpa_major_version ${CMAKE_MATCH_1})
string(REGEX MATCH "LIBRPA_MINOR_VERSION ([0-9]+)" librpa_minor_version ${version_file})
set(librpa_minor_version ${CMAKE_MATCH_1})
string(REGEX MATCH "LIBRPA_MICRO_VERSION ([0-9]+)" librpa_micro_version ${version_file})
set(librpa_micro_version ${CMAKE_MATCH_1})
set(version "${spglib_major_version}.${spglib_minor_version}.${spglib_micro_version}")

set(DRIVER_NAME "chi0_main.exe")
set(LIB_NAME    "rpa")

# options setup
option(USE_LIBRI     "Use LibRI for tensor contraction" OFF)
option(USE_CMAKE_INC "Use cmake.inc for configure" OFF)
option(ENABLE_DEBUG  "Flag of debug mode " OFF)
option(ENABLE_DOC    "Flag to build documentation" OFF)
option(ENABLE_TEST   "Flag to build unit tests" OFF)
option(ENABLE_DRIVER "Flag to build driver executables" ON)

# Add cmake load files
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

set(math_libs "")
set(parallel_libs "")
set(_loaded_cmake_inc OFF)

# TODO: maybe better way to customize environments?
if(EXISTS "${PROJECT_SOURCE_DIR}/cmake.inc" AND USE_CMAKE_INC)
  include("${PROJECT_SOURCE_DIR}/cmake.inc")
  list(APPEND math_libs ${LAPACK_LIBRARIES})
  set(_loaded_cmake_inc ON)
  message(STATUS "Using customize environemnt variables from ${PROJECT_SOURCE_DIR}/cmake.inc")
endif()

# C++ standard requirement
if(USE_LIBRI)
  set(CMAKE_CXX_STANDARD 14)
else()
  set(CMAKE_CXX_STANDARD 11)
endif()
set(CMAKE_CXX_STANDARD_REQUIRED True)

# enable_language later than cmake.inc to take the variables therein into accout
enable_language(CXX)

# for installation
include(GNUInstallDirs)

if(NOT _loaded_cmake_inc)
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

  # message(INFO "${CMAKE_C_COMPILER_LOADED} ${CMAKE_CXX_COMPILER_LOADED} ${CMAKE_Fortran_COMPILER_LOADED}") 
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
  add_compile_definitions("LIBRPA_USE_LIBRI")
  message(STATUS "The project will be compiled with LibRI")
endif()

if(ENABLE_DEBUG)
  message(STATUS "Debug mode switched on")
  add_compile_definitions("LIBRPA_DEBUG")
endif()

if(ENABLE_DRIVER)
  message(STATUS "Build driver")
  add_subdirectory(driver)
endif()

add_subdirectory(src)

if(ENABLE_DOC)
  include(UseDoxygenDoc)
  add_subdirectory(docs)
endif()

if(ENABLE_TEST)
  enable_testing()
  add_subdirectory(tests)
endif()
