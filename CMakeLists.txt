# main CMake file for LibRPA
cmake_minimum_required(VERSION 3.0)

# extract version number
file(READ "src/Makefile" mkfile)
string(REGEX MATCH " *VERSION *= *[0-9.]*" version ${mkfile})
string(REGEX REPLACE " *VERSION *= *" "" version ${version})
# Project information
project(LibRPA VERSION ${version} LANGUAGES NONE)
# Add cmake load files
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
include(UseDoxygenDoc)

# C++ standard requirement
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if(EXISTS "${PROJECT_SOURCE_DIR}/cmake.inc")
  include("${PROJECT_SOURCE_DIR}/cmake.inc")
endif()

enable_language(CXX)

add_subdirectory(src)

add_subdirectory(docs)
