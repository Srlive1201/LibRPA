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
option(ENABLE_TEST   "Flag to build unit tests" ON)
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
  message(STATUS "Build with LibRI and dependeicies")
  find_package(Cereal)
  if(NOT Cereal_FOUND)
    message(WARNING "User Cereal not found, use thirdparty/cereal-1.3.0")
    set(Cereal_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/./thirdparty/cereal-1.3.0/include)
  endif()
  find_package(LibComm)
  if(NOT LibComm_FOUND)
    message(WARNING "User LibComm not found, use thirdparty/LibComm")
    set(LibComm_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/./thirdparty/LibComm/include)
  endif()
  find_package(LibRI)
  if(NOT LibRI_FOUND)
    message(WARNING "User LibRI not found, use thirdparty/LibRI")
    set(LibRI_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/./thirdparty/LibRI/include)
  endif()
  add_compile_definitions("LIBRPA_USE_LIBRI")
endif()

if(ENABLE_DEBUG)
  message(STATUS "Debug mode switched on")
  add_compile_definitions("LIBRPA_DEBUG")
endif()

if(ENABLE_TEST)
  # FetchContent_Declare(
  #   googletest
  #   GIT_REPOSITORY https://github.com/google/googletest.git
  #   GIT_TAG release-1.12.1
  # )
  # FetchContent_MakeAvailable(googletest)
  # include(GoogleTest)
  message(STATUS "Test enabled")
  enable_testing()
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
