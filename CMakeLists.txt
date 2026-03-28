# main CMake file for LibRPA
cmake_minimum_required(VERSION 3.13)

project(librpa NONE)

# extract version number from version.h file
file(READ ${CMAKE_CURRENT_SOURCE_DIR}/src/version.h version_file)
string(REGEX MATCH "LIBRPA_MAJOR_VERSION ([0-9]+)" librpa_major_version ${version_file})
set(librpa_major_version ${CMAKE_MATCH_1})
string(REGEX MATCH "LIBRPA_MINOR_VERSION ([0-9]+)" librpa_minor_version ${version_file})
set(librpa_minor_version ${CMAKE_MATCH_1})
string(REGEX MATCH "LIBRPA_PATCH_VERSION ([0-9]+)" librpa_patch_version ${version_file})
set(librpa_patch_version ${CMAKE_MATCH_1})
set(version "${librpa_major_version}.${librpa_minor_version}.${librpa_patch_version}")

set(CMAKE_PROJECT_VERSION "${version}")

set(LIBRPA_DRIVER_NAME         "chi0_main.exe")

# Public options setup
option(LIBRPA_USE_LIBRI            "Use LibRI for tensor contraction in LibRPA" OFF)
option(LIBRPA_USE_CMAKE_INC        "Use cmake.inc for configure" OFF)
option(LIBRPA_USE_EXTERNAL_GREENX  "Use an external GreenX target provided by the parent project instead of the bundled copy" OFF)
option(LIBRPA_ENABLE_TEST          "Flag to build LibRPA tests" ON)
option(LIBRPA_ENABLE_DRIVER        "Flag to build LibRPA driver executable" ON)
option(LIBRPA_ENABLE_FORTRAN_BIND  "Flag to build LibRPA Fotran binding" OFF)
option(LIBRPA_VERBOSE_OUTPUT       "Flag to print verbose information in stdout and process output of LibRPA" OFF)
# NOTE: static library not tested
# option(BUILD_LIBRPA_SHARED  "Flag to build shared libraries" ON)

option(LIBRPA_ENABLE_CPP_TEST      "Flag to build tests for C++ components" ON)
option(LIBRPA_ENABLE_FORTRAN_TEST  "Flag to build tests for Fortran bindings" ON)

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
set(_BUNDLED_GREENX_COMMIT_HASH  "c90d131")
set(_BUNDLED_LIBRI_COMMIT_HASH   "1f4200c")
set(_BUNDLED_LIBCOMM_COMMIT_HASH "965bf90")

# Add cmake load files
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

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
# Set to install/ when LibRPA is the root project and
# CMAKE_INSTALL_PREFIX was not manually specified
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT AND (CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR))
  set(default_install_prefix "${CMAKE_SOURCE_DIR}/install")  # default to install at the root project
  set(CMAKE_INSTALL_PREFIX ${default_install_prefix}
      CACHE STRING "The installation directory is not set, default to ${default_install_prefix}"
      FORCE)
  set(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT OFF)
endif()

# Get git commit number when LibRPA is the root project
if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
  include(GetGitRevisionDescription)
  get_git_head_revision(LIBRPA_GIT_REF LIBRPA_GIT_SHA1)
else()
  # Otherwise get the checkpoint recorded in the root project
  if (NOT DEFINED LIBRPA_GIT_REF)
    set(LIBRPA_GIT_REF "GIT-REF-UNKNOWN")
  endif()
  if (NOT DEFINED LIBRPA_GIT_SHA1)
    set(LIBRPA_GIT_SHA1 "GIT-SHA1-UNKNOWN")
  endif()
endif()

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

# Switch on LibRI functionalities, which requires Cereal, LibComm and LibRI
# When USE_LIBRI is switched on but the include direcotires of the required dependenccies are not specified
# by command line or environment variables, the internal packages are used.
# If any of the include directory is specified but not found, a fatal error will be raised.
if(LIBRPA_USE_LIBRI)
  message(STATUS "Build with LibRI and dependencies")

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
    message(STATUS "Use internal LibComm : thirdparty/LibComm (commit hash ${_BUNDLED_LIBCOMM_COMMIT_HASH})")
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
    message(STATUS "Use internal LibRI   : thirdparty/LibRI   (commit hash ${_BUNDLED_LIBRI_COMMIT_HASH})")
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
  add_subdirectory(thirdparty/greenX)
else()
  set(LIBRPA_GREENX_COMMIT_HASH "(EXTERNAL)")
endif()

if(LIBRPA_ENABLE_TEST)
  # FetchContent_Declare(
  #   googletest
  #   GIT_REPOSITORY https://github.com/google/googletest.git
  #   GIT_TAG release-1.12.1
  # )
  # FetchContent_MakeAvailable(googletest)
  # include(GoogleTest)
  message(STATUS "LibRPA test enabled")
  enable_testing()
endif()

if(LIBRPA_ENABLE_FORTRAN_BIND)
  message(STATUS "LibRPA Fortran binding enabled")
  add_subdirectory(binding/fortran)
endif()

if(LIBRPA_ENABLE_DRIVER)
  message(STATUS "LibRPA driver enabled")
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
  ${CMAKE_CURRENT_SOURCE_DIR}/include/librpa_options.h
  ${CMAKE_CURRENT_SOURCE_DIR}/include/librpa.h
  ${CMAKE_CURRENT_SOURCE_DIR}/include/librpa.hpp
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/librpa)

message(STATUS "CMAKE_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX}")
message(STATUS "CMAKE_INSTALL_RPATH: ${CMAKE_INSTALL_RPATH}")
message(STATUS "CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")
message(STATUS "CMAKE_Fortran_FLAGS: ${CMAKE_Fortran_FLAGS}")
