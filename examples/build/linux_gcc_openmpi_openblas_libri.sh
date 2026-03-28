#!/bin/sh

# This script uses GCC, OpenMPI, OpenBLAS and ScaLAPACK to build LibRPA on Linux.
# LibRI is enabled using the bundled version. Fortran binding also enabled.
#
# Tested on Fedora 43.
#
# All dependencies (except for those bundled) are installed from dnf package manager.
#

BUILDDIR="${BUILDDIR:=build_gcc_openmpi_scalapack}"

# The following two variables need customize
export OPENBLAS_DIR="/usr/lib64"
export SCALAPACK_DIR="/usr/lib64/openmpi/lib"

# Compilers setup
export CXX=mpicxx
export FC=mpifort
export CC=gcc
export OMPI_CC=gcc
export OMPI_CXX=g++
export OMPI_FC=gfortran

# Compile flags
export FCFLAGS="-fPIC -g -fallow-argument-mismatch -ffree-line-length-none"
export CXXFLAGS="-fPIC -g -rdynamic -Wl,-lgfortran"

cmake -B $BUILDDIR \
  -DLIBRPA_USE_LIBRI=ON \
  -DLIBRPA_ENABLE_TEST=ON \
  -DLIBRPA_ENABLE_FORTRAN_BIND=ON \
  -DLIBRPA_FORTRAN_DP=c_double \
  -DBLAS_LIBRARIES="-L${OPENBLAS_DIR} -lopenblas" \
  -DLAPACK_LIBRARIES="-L${OPENBLAS_DIR} -lopenblas" \
  -DSCALAPACK_DIR="$SCALAPACK_DIR" \
  -DCMAKE_CXX_FLAGS="$CXXFLAGS" \
  -DCMAKE_Fortran_FLAGS="$FCFLAGS"

cmake --build $BUILDDIR -j 4
