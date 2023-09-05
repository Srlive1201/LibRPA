# Basic C++ compiler setup
# set(CMAKE_C_COMPILER "mpiicc" CACHE STRING "" FORCE)
set(CMAKE_CXX_COMPILER "mpiicpc" CACHE STRING "" FORCE)

set(CMAKE_CXX_FLAGS "-g -O2 -fopenmp" CACHE STRING "" FORCE)
if(DEBUG)
  set(CMAKE_CXX_FLAGS "-g -O0 -Wall -fopenmp" CACHE STRING "" FORCE)
endif()


# set(LAPACK_LIBRARIES "-L$ENV{MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core")
set(LAPACK_LIBRARIES "-L$ENV{MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lmkl_blacs_intelmpi_lp64 -lpthread -lm -ldl")

# whether to use LibRI for tensor contraction
# set(USE_LIBRI ON)
