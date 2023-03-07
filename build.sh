#!/usr/bin/env bash

TEST=0
LIBRI=0
DOCS=0
DEBUG=0
BUILDDIR=build
VERBOSE=0
CXX=mpiicpc

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
  echo " -c <CXX_COMP>  : set C++ compiler (default mpiicpc)"
  echo " --bd <builddir>: set build directory"
  echo " --ri           : use libRI"
}

while [ ${#} -gt 0 ]; do
  case "$1" in
  -c   ) CXX=$1; shift 1;;
  -t   ) TEST=1;;
  --ri ) LIBRI=1;;
  -d   ) DOCS=1;;
  -D   ) DEBUG=1;;
  --bd ) BUILDDIR=$1; shift 1;;
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
  options="-DENABLE_DEBUG=ON $options"
else
  options="-DENABLE_DEBUG=OFF $options"
fi

CXX=$CXX cmake -B "$BUILDDIR" $options
if (( VERBOSE )); then
  cmake --build "$BUILDDIR" -v
else
  cmake --build "$BUILDDIR"
fi
