# Compile Options

## Overview

| Option                                        | Type | Default |
|-----------------------------------------------|------|---------|
| [`USE_LIBRI`](#use-libri)                     | Bool | `OFF`   |
| [`USE_CMAKE_INC`](#use_cmake_inc)             | Bool | `OFF`   |
| [`USE_GREENX_API`](#use-greenx-api)           | Bool | `OFF`   |
| [`ENABLE_FORTRAN_BIND`](#enable-fortran-bind) | Bool | `OFF`   |

## `USE_LIBRI`

When switching on, the code will be compiled with [LibRI](https://github.com/abacusmodeling/LibRI)
to handle contraction of RI tensors.

Note that the *GW* and EXX functionality requires the code compiled with LibRI, i.e. `-DUSE_LIBRI=ON`.
RPA correlation energy can be computed without this flag on.

## `USE_CMAKE_INC`

When switched on, the `cmake.inc` file will be used to initialize the compilers and other options.

This option would be deprecated in the future because it can be replaced by parsing `cmake.inc` to
the command line option `-C` of cmake.

## `USE_GREENX_API`

The minimax grids are part of the [Green X](https://nomad-coe.github.io/greenX/) library.
When `OFF`, the plain-text minimax grids stored under `src/minimax_grid/GreenX` will be used.
The transform coefficients are then calculated by calling an Python script.
Switching to `ON` will make the code link to the GreenX library and call the API to generate the minimax grids and transform matrices.

In principle these two ways to get the minimax grids should be essentially the same.
However, the plain-text grids were extracted at the early stage of the Green X library,
and the grids can be missing for certain energy range and number of grid points.
Thus it is recommended to use the API.

## `ENABLE_FORTRAN_BIND`

When swicthed on, the Fortran binding of LibRPA will be built.
