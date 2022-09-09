# control file of building LibRPA
# This file will be included when using either make or cmake.

CXX = mpiicpc

LIBCOMM_DIR=/home/linpz/software/LibComm/include
LIBRI_DIR=/home/linpz/software/libRI/include
CEREAL_DIR=/home/linpz/software/cereal-1.3.0/include

# definitions for preprocessor
#   -D__USE_LIBRI : use libRI for tensor contraction
CPPFLAGS = -D__USE_LIBRI

INCLUDES= -I. -I${CEREAL_DIR} -I$(LIBCOMM_DIR) -I$(LIBRI_DIR)

# C++ compiler flags
# production
CXXFLAGS = -g -O2 -std=c++11 -qopenmp
# # debug
# CXXFLAGS = -g -O0 -Wall -std=c++11 -qopenmp

# LAPACK library
#LAPACK_LIBRARIES = -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
LAPACK_LIBRARIES = -L$(MKLROOT)/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lmkl_blacs_intelmpi_lp64 -lpthread -lm -ldl

# # GNU
# CXX = mpicxx
# CXXFLAGS = -g -O2 -std=c++11 -fopenmp
# # CXXFLAGS = -g -O0 -Wall -std=c++11 -fopenmp
# LAPACK_LIBRARIES = -L$(HOME)/software/lapack-3.9.1 -llapack -lblas -lgfortran
