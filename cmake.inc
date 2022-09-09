# Basic C++ compiler setup
# set(CMAKE_C_COMPILER "mpiicc" CACHE STRING "" FORCE)
set(CMAKE_CXX_COMPILER "mpiicpc" CACHE STRING "" FORCE)

set(CMAKE_CXX_FLAGS "-g -O2 -qopenmp" CACHE STRING "" FORCE)

# set(LAPACK_LIBRARIES "-L$ENV{MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core")
set(LAPACK_LIBRARIES "-L$ENV{MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lmkl_blacs_intelmpi_lp64 -lpthread -lm -ldl")

# whether to use LibRI for tensor contraction
set(USE_LIBRI OFF)

# whether to build a library
set(BUILD_LIBRARY OFF)
