#!/bin/sh

# This script uses Intel LLVM-based compilers, Intel MPI and MKL to build LibRPA
# on Linux platform for develop and production use. Tested on Ubuntu and SUSE.

# Intel compilers icpx (C++) and ifx (Fortran) needs to be found under directories
# in environment variable PATH, as well as the C++ MPI wrapper mpiicpx.
# Note that Fortran is used only when USE_GREENX_API or ENABLE_FORTRAN_BIND is on.
# Intel MKL will be used as the working math library.

BUILDDIR="${BUILDDIR:=build_intel_llvm_mkl}"

# # Ensure environment variables are correctly set.
# # On PC, they can be set by sourcing setup script provided by the vendor.
# # As LLVM-based compilers are only provided in oneAPI toolchain, this can be done by:
# source /opt/intel/oneapi/setvars.sh

# # On HPC, usually environment modules are provided:
# module load intel/2023.1.0.x impi/2021.9 mkl/2023.1

# # Or you might set it manually, which is not recommended.

# The GreenX API is switched on. In this case, the Fortran compiler needs to be configured.
export CXX=mpiicpx
export FC=ifx

cmake -B $BUILDDIR -DENABLE_TEST=ON -DUSE_LIBRI=OFF -DUSE_GREENX_API=ON
cd $BUILDDIR && make -j 4
