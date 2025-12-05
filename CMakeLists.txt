# main CMake file for LibRPA
cmake_minimum_required(VERSION 3.13)

project(LibRPA NONE)

# extract version number from version.h file
file(READ ${PROJECT_SOURCE_DIR}/src/version.h version_file)
string(REGEX MATCH "LIBRPA_MAJOR_VERSION ([0-9]+)" librpa_major_version ${version_file})
set(librpa_major_version ${CMAKE_MATCH_1})
string(REGEX MATCH "LIBRPA_MINOR_VERSION ([0-9]+)" librpa_minor_version ${version_file})
set(librpa_minor_version ${CMAKE_MATCH_1})
string(REGEX MATCH "LIBRPA_PATCH_VERSION ([0-9]+)" librpa_patch_version ${version_file})
set(librpa_patch_version ${CMAKE_MATCH_1})
set(version "${librpa_major_version}.${librpa_minor_version}.${librpa_patch_version}")

set(CMAKE_PROJECT_VERSION "${version}")

set(LIBRPA_DRIVER_NAME         "chi0_main.exe")
set(LIBRPA_LIB_NAME            "rpa")
set(LIBRPA_FORTRAN_LIB_NAME    "rpa_f03")

# Public options setup
option(LIBRPA_USE_LIBRI            "Use LibRI for tensor contraction in LiBRPA" OFF)
option(LIBRPA_USE_CMAKE_INC        "Use cmake.inc for configure" OFF)
option(LIBRPA_USE_EXTERNAL_GREENX  "Use external GreenX library rather than the packaged one in LibRPA" OFF)
option(ENABLE_LIBRPA_TEST          "Flag to build LibRPA tests" ON)
option(ENABLE_LIBRPA_DRIVER        "Flag to build LibRPA driver executable" ON)
option(ENABLE_LIBRPA_FORTRAN_BIND  "Flag to build LibRPA Fotran binding" OFF)
option(LIBRPA_VERBOSE_OUTPUT       "Flag to print verbose information in stdout and process output of LibRPA" OFF)
# NOTE: static library not tested
# option(BUILD_LIBRPA_SHARED  "Flag to build shared libraries" ON)

option(ENABLE_LIBRPA_CPP_TEST      "Flag to build tests for C++ components" ON)
option(ENABLE_LIBRPA_FORTRAN_TEST  "Flag to build tests for Fortran bindings" ON)

# Customizable cmake macro
set(LIBRPA_FORTRAN_DP "c_double" CACHE STRING "Kind of double precision that the Fortran binding uses to interface the user code")

if (DEFINED USE_GREENX_MINIMAX)
  message(FATAL_ERROR "USE_GREENX_MINIMAX is deprecated: GreenX is now always used. Please remove USE_GREENX_MINIMAX")
endif()
if (DEFINED USE_GREENX_API)
  if (NOT USE_GREENX_API)
    message(FATAL_ERROR "USE_GREENX_API is deprecated: GreenX is now always used. Please remove USE_GREENX_API")
  endif()
endif()
# Hash number of the bundled greenX. Keep it consistent with thirdparty/greenX
set(_BUNDLED_GREENX_COMMIT_HASH "c90d131")

# Add cmake load files
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

set(math_libs "")
set(parallel_libs "")
set(_loaded_cmake_inc OFF)

# TODO: maybe better way to customize environments?
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/cmake.inc" AND LIBRPA_USE_CMAKE_INC)
  include("${CMAKE_CURRENT_SOURCE_DIR}/cmake.inc")
  list(APPEND math_libs ${LAPACK_LIBRARIES})
  set(_loaded_cmake_inc ON)
  message(STATUS "Using customize environemnt variables from ${CMAKE_CURRENT_SOURCE_DIR}/cmake.inc")
endif()

# C++ standard requirement
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# enable_language later than cmake.inc to take the variables therein into accout
# require Fortran for GreenX library and Fortran binding
enable_language(CXX Fortran)

# bypass the deprecation warning of classic C++ of oneAPI
if(CMAKE_CXX_COMPILER_ID MATCHES Intel)
  if(CMAKE_CXX_COMPILER MATCHES "mpiicpc" OR CMAKE_CXX_COMPILER MATCHES "icpc")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL "2021.7.0.20221019")
      message(STATUS "Intel C++ Compiler Classic detected, adding -diag-disable=10441 flag")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -diag-disable=10441")
    endif()
  endif()
endif()

# Setup installation paths
include(GNUInstallDirs)
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(default_install_prefix "${CMAKE_BINARY_DIR}/install")
  set(CMAKE_INSTALL_PREFIX ${default_install_prefix}
      CACHE STRING "The installation directory is not set, default to ${default_install_prefix}"
      FORCE)
endif()

# Get git commit number
include(GetGitRevisionDescription)
get_git_head_revision(GIT_REF GIT_SHA1)

if(NOT _loaded_cmake_inc)
  if(DEFINED ENV{MKLROOT} AND NOT DEFINED MKLROOT)
    set(MKLROOT "$ENV{MKLROOT}")
  endif()
  if(MKLROOT)
    find_package(IntelMKL REQUIRED)
    # include_directories(${MKL_INCLUDE_DIRS} ${MKL_INCLUDE_DIRS}/fftw)
    include_directories(${MKL_INCLUDE_DIRS})
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

# TODO: always switch on MPI, until test program for MPI is implemented
add_compile_definitions("LIBRPA_USE_MPI")

include_directories(src)

# Switch on LibRI functionalities, which requires Cereal, LibComm and LibRI
# When USE_LIBRI is switched on but the include direcotires of the required dependenccies are not specified
# by command line or environment variables, the internal packages are used.
# If any of the include directory is specified but not found, a fatal error will be raised.
if(LIBRPA_USE_LIBRI)
  message(STATUS "Build with LibRI and dependeicies")

  # Check Cereal
  if(NOT CEREAL_INCLUDE_DIR)
    set(CEREAL_INCLUDE_DIR $ENV{CEREAL_INCLUDE_DIR})
  endif()
  if("${CEREAL_INCLUDE_DIR}" STREQUAL "")
    message(STATUS "Use internal Cereal  : thirdparty/cereal-1.3.0")
    set(Cereal_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/cereal-1.3.0/include)
  else()
    find_package(Cereal)
    if(NOT Cereal_FOUND)
      message(FATAL_ERROR "User Cereal include directory set to ${Cereal_INCLUDE_DIR}, but not found")
    else()
      message(STATUS "Use user Cereal include directory  : ${Cereal_INCLUDE_DIR}")
    endif()
  endif()

  # Check LibComm
  if(NOT LIBCOMM_INCLUDE_DIR)
    set(LIBCOMM_INCLUDE_DIR $ENV{LIBCOMM_INCLUDE_DIR})
  endif()
  if("${LIBCOMM_INCLUDE_DIR}" STREQUAL "")
    message(STATUS "Use internal LibComm : thirdparty/LibComm")
    set(LibComm_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/LibComm/include)
  else()
    find_package(LibComm)
    if(NOT LibComm_FOUND)
      message(FATAL_ERROR "User LibComm include directory set to ${LIBCOMM_INCLUDE_DIR}, but not found")
    else()
      message(STATUS "Use user LibComm include directory : ${LibComm_INCLUDE_DIR}")
    endif()
  endif()

  # Check LibRI
  if(NOT LIBRI_INCLUDE_DIR)
    set(LIBRI_INCLUDE_DIR $ENV{LIBRI_INCLUDE_DIR})
  endif()
  if("${LIBRI_INCLUDE_DIR}" STREQUAL "")
    message(STATUS "Use internal LibRI   : thirdparty/LibRI")
    set(LibRI_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/LibRI/include)
  else()
    find_package(LibRI)
    if(NOT LibRI_FOUND)
      message(FATAL_ERROR "User LibRI include directory set to ${LIBRI_INCLUDE_DIR}, but not found")
    else()
      message(STATUS "Use user LibRI include directory   : ${LibRI_INCLUDE_DIR}")
    endif()
  endif()

  add_compile_definitions("LIBRPA_USE_LIBRI")
endif()

if(CMAKE_BUILD_TYPE MATCHES "Debug")
  message(STATUS "Build type set to Debug, adding LIBRPA_DEBUG preprocessor directive")
  add_compile_definitions("LIBRPA_DEBUG")
  set(LIBRPA_VERBOSE_OUTPUT ON)
endif()

if(LIBRPA_VERBOSE_OUTPUT)
  add_compile_definitions("LIBRPA_VERBOSE")
endif()

# Switch on GreenX
if(NOT LIBRPA_USE_EXTERNAL_GREENX)
  set(LIBRPA_GREENX_COMMIT_HASH "${_BUNDLED_GREENX_COMMIT_HASH}")
  message(STATUS "Generate MiniMax grids by using the bundled GreenX library (collected from commit ${LIBRPA_GREENX_COMMIT_HASH})")
  add_subdirectory(./thirdparty/greenX)
else()
  set(LIBRPA_GREENX_COMMIT_HASH "(EXTERNAL)")
endif()

if(ENABLE_LIBRPA_FORTRAN_BIND)
  message(STATUS "LibRPA Fortran binding enabled")
  add_subdirectory(./binding/fortran)
endif()

if(ENABLE_LIBRPA_TEST)
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

if(ENABLE_LIBRPA_DRIVER)
  message(STATUS "Build driver")
  list(APPEND CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})
  add_subdirectory(driver)
endif()

add_subdirectory(src)
install(FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/include/librpa_compute.h
  ${CMAKE_CURRENT_SOURCE_DIR}/include/librpa_enums.h
  ${CMAKE_CURRENT_SOURCE_DIR}/include/librpa_global.h
  ${CMAKE_CURRENT_SOURCE_DIR}/include/librpa_handler.h
  ${CMAKE_CURRENT_SOURCE_DIR}/include/librpa_input.h
  ${CMAKE_CURRENT_SOURCE_DIR}/include/librpa_func.h
  ${CMAKE_CURRENT_SOURCE_DIR}/include/librpa_options.h
  ${CMAKE_CURRENT_SOURCE_DIR}/include/librpa.h
  ${CMAKE_CURRENT_SOURCE_DIR}/include/librpa.hpp
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

message(STATUS "CMAKE_INSTALL_RPATH: ${CMAKE_INSTALL_RPATH}")
message(STATUS "CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")
message(STATUS "CMAKE_Fortran_FLAGS: ${CMAKE_Fortran_FLAGS}")
