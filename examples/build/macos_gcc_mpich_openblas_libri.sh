#!/bin/sh

# This script uses GCC, MPICH, OpenBLAS and ScaLAPACK to build LibRPA
# on macOS for develop and test use. LibRI is enabled using the bundled version.
# Tested for Sequoia, Sonoma and Ventura.

# GCC compilers (g++, gfortran), MPICH and OpenBLAS are installed using homebrew.
# ScaLAPACK needs to be installed manually, because the homebrew ScaLAPACK is
# built on OpenMPI.

BUILDDIR="${BUILDDIR:=build_macos_gcc_mpich}"

# The following two variables need customize
export OPENBLAS_DIR="/opt/homebrew/Cellar/openblas/0.3.28/lib"
export SCALAPACK_DIR="/opt/packages/scalapack/2.2.0/gcc-14.2.0-mpich-4.2.2-openblas"

# Compilers setup
export CXX=mpicxx
export MPICH_CXX=g++-14
export CC=gcc-14
export FC=gfortran-14

# Note: -Wl,-ld_classic may not be necessary for Sequoia.
export FCFLAGS="-Wl,-ld_classic"
export CXXFLAGS="-cxx=g++-14 -g -rdynamic -Wl,-lgfortran -Wl,-ld_classic"

cmake -B $BUILDDIR \
  -DLIBRPA_USE_LIBRI=ON \
  -DLIBRPA_ENABLE_TEST=ON \
  -DBLAS_LIBRARIES="-L${OPENBLAS_DIR} -lblas" \
  -DLAPACK_LIBRARIES="-L${OPENBLAS_DIR} -llapack" \
  -DSCALAPACK_DIR="$SCALAPACK_DIR" \
  -DCMAKE_CXX_FLAGS="$CXXFLAGS" \
  -DCMAKE_Fortran_FLAGS="$FCFLAGS"

cmake --build $BUILDDIR -j 4
