# Basic C++ compiler setup
# set(CMAKE_C_COMPILER "mpiicc" CACHE STRING "" FORCE)
set(CMAKE_CXX_COMPILER "mpiicpc" CACHE STRING "" FORCE)

set(CMAKE_CXX_FLAGS "-g -O2 -qopenmp" CACHE STRING "" FORCE)

set(LAPACK_LIBRARIES "-L$ENV{MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core")

set(USE_LIBRARY OFF)
