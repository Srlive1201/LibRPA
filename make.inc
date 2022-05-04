# control file of building LibRPA
# This file will be included when using either make or cmake.

CXX = icpc
MPICXX = mpiicpc

# C++ compiler flags
# production
CXXFLAGS = -g -O2 -std=c++11 -qopenmp
# # debug
# CXXFLAGS = -g -O0 -Wall -std=c++11 -qopenmp

# Lapack library
LAPACK_DIR = -L$(MKLROOT)/lib/intel64
LAPACK_LIB = -lmkl_intel_lp64 -lmkl_sequential -lmkl_core

