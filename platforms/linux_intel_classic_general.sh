#!/bin/sh

# This script uses Intel classic compilers, Intel MPI and MKL to build LibRPA
# on Linux platform for develop and production use. Tested on Ubuntu, CentOS and SUSE.

# Intel compilers icpc (C++) and ifort (Fortran) needs to be found under directories
# in environment variable PATH, as well as the C++ MPI wrapper mpiicpc.
# Note that Fortran is used only when USE_GREENX_API or ENABLE_FORTRAN_BIND is on.
# Intel MKL will be used as the working math library.

BUILDDIR="${BUILDDIR:=build_intel_classic_mkl}"

# # Ensure environment variables are correctly set.
# # On PC, they can be set by sourcing setup script provided by the vendor.
# # For oneAPI:
# source /opt/intel/oneapi/setvars.sh
# # For old Intel versions
# source /opt/intel/compilers_and_libraries_2020.4.304/linux/bin/compilervars.sh intel64

# # On HPC, usually environment modules are provided.
# # Examples:
# #   MPCDF platforms, oneAPI 2023
# module load intel/2023.1.0.x impi/2021.9 mkl/2023.1

# #   MPCDF platforms, Intel 2020 update 4
# module load intel/19.1.3 impi/2019.9 mkl/2020.4

# # Or you might set it manually, which is not recommended.

# The GreenX API is switched on. In this case, the Fortran compiler needs to be configured.
export CXX=mpiicpc
export FC=ifort

cmake -B $BUILDDIR -DENABLE_TEST=ON -DUSE_LIBRI=OFF -DUSE_GREENX_API=ON
cd $BUILDDIR && make -j 4
