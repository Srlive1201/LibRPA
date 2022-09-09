# main CMake file for LibRPA
cmake_minimum_required(VERSION 3.0)

# extract version number
file(READ "src/Makefile" mkfile)
string(REGEX MATCH " *VERSION *= *[0-9.]*" version ${mkfile})
string(REGEX REPLACE " *VERSION *= *" "" version ${version})
# Project information
project(LibRPA VERSION ${version} LANGUAGES NONE)


# options setup
option(USE_LIBRI "Use LibRI for tensor contraction" OFF)

# Add cmake load files
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
include(UseDoxygenDoc)

# C++ standard requirement
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if(EXISTS "${PROJECT_SOURCE_DIR}/cmake.inc")
  include("${PROJECT_SOURCE_DIR}/cmake.inc")
endif()

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

enable_language(CXX)

add_subdirectory(src)

add_subdirectory(docs)

enable_testing()
add_subdirectory(tests)
