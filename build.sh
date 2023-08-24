#!/usr/bin/env bash

TEST=0
LIBRI=0
DOCS=0
DEBUG=0
BUILDDIR=build
VERBOSE=0
CXX=mpiicpc
FC=mpiifort
GREENX=0

CMAKE=cmake3
# Some build examples
# with LibRI and GreenX minimax
#     CXX=mpiicpc FC=mpiifort LIBRI_INCLUDE_DIR=/home/minyez/projects/LibRI/include cmake -B build -DUSE_GREENX_MINIMAX=ON -DENABLE_TEST=OFF -DUSE_LIBRI=ON

help_info()
{
  echo "Help script to build libRPA"
  echo "Usage: $0 [-t] [--ri] [-d] [-D] [--bd] [-v] [-h] [-c]"
  echo ""
  echo "Options:"
  echo " -h             : print this information and exit"
  echo " -t             : enable test"
  echo " -d             : enable doc"
  echo " -D             : enable debug"
  echo " -v             : verbose"
  echo " --gx           : use GreenX library for minimax grids"
  echo " -c <CXX_COMP>  : set C++ compiler (default mpiicpc)"
  echo " -f <FC_COMP>   : set Fortran compiler, used only with GreenX (default mpiifort)"
  echo " --bd <builddir>: set build directory (default 'build')"
  echo " --ri           : use libRI"
}

while [ ${#} -gt 0 ]; do
  case "$1" in
  -c   ) CXX=$2; shift 1;;
  -t   ) TEST=1;;
  -f   ) FC=$2; shift 1;;
  --ri ) LIBRI=1;;
  --gx ) GREENX=1;;
  -d   ) DOCS=1;;
  -D   ) DEBUG=1;;
  --bd ) BUILDDIR=$2; shift 1;;
  -v   ) VERBOSE=1;;
  -h | --help | help | h ) help_info; exit 0;;
  * ) echo "Unknown option $0"; exit 1;;
  esac
  shift 1
done

options=""
if (( TEST )); then
  options="-DENABLE_TEST=ON $options"
else
  options="-DENABLE_TEST=OFF $options"
fi

if (( DOCS )); then
  options="-DENABLE_DOC=ON $options"
else
  options="-DENABLE_DOC=OFF $options"
fi

if (( DEBUG )); then
  options="-DCMAKE_BUILD_TYPE=Debug $options"
fi

if (( LIBRI )); then
  options="-DUSE_LIBRI=ON $options"
fi

if (( GREENX )); then
  options="-DUSE_GREENX_MINIMAX=ON $options"
fi

echo "CXX=$CXX cmake -B $BUILDDIR $options"
CXX="$CXX" $CMAKE -B "$BUILDDIR" $options
if (( VERBOSE )); then
  $CMAKE --build "$BUILDDIR" -v
else
  $CMAKE --build "$BUILDDIR"
fi
