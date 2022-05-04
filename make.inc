# control file of building LibRPA
# This file will be included when using either make or cmake.

CXX = mpiicpc

# C++ compiler flags
# production
CXXFLAGS = -g -O2 -std=c++11 -qopenmp
# # debug
# CXXFLAGS = -g -O0 -Wall -std=c++11 -qopenmp

# LAPACK library
LAPACK_LIBRARIES = -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core

