# Compile Options

## Overview

| Option                                        | Type   | Default |
|-----------------------------------------------|--------|---------|
| [`USE_LIBRI`](#use-libri)                     | Bool   | `OFF`   |
| [`USE_CMAKE_INC`](#use-cmake-inc)             | Bool   | `OFF`   |
| [`USE_GREENX_API`](#use-greenx-api)           | Bool   | `OFF`   |
| [`ENABLE_FORTRAN_BIND`](#enable-fortran-bind) | Bool   | `OFF`   |
| [`ENABLE_DRIVER`](#enable-driver)             | Bool   | `ON`    |
| [`ENABLE_TEST`](#enable-test)                 | Bool   | `ON`    |
| [`LIBRI_INCLUDE_DIR`](#libri-include-dir)     | String | empty   |
| [`LIBCOMM_INCLUDE_DIR`](#libcomm-include-dir) | String | empty   |

The options should be parsed as cmake command line options, for example, `-DUSE_LIBRI=ON`.

(use-libri)=
## `USE_LIBRI`

When switching on, the code will be compiled with [LibRI](https://github.com/abacusmodeling/LibRI)
to handle contraction of RI tensors.

Note that the *GW* and EXX functionality requires the code compiled with LibRI, i.e. `-DUSE_LIBRI=ON`.
RPA correlation energy can be computed without this flag on.

(use-cmake-inc)=
## `USE_CMAKE_INC`

When switched on, the `cmake.inc` file will be used to initialize the compilers and other options.

This option would be deprecated in the future because it can be replaced by parsing `cmake.inc` to
the command line option `-C` of cmake.

(use-greenx-api)=
## `USE_GREENX_API`

The minimax grids are part of the [Green X](https://nomad-coe.github.io/greenX/) library.
When `OFF`, the plain-text minimax grids stored under `src/minimax_grid/GreenX` will be used.
The transform coefficients are then calculated by calling an Python script.
Switching to `ON` will make the code link to the GreenX library and call the API to generate the minimax grids and transform matrices.

In principle these two ways to get the minimax grids should be essentially the same.
However, the plain-text grids were extracted at the early stage of the Green X library,
and the grids can be missing for certain energy range and number of grid points.
Thus it is recommended to use the API.

(enable-fortran-bind)=
## `ENABLE_FORTRAN_BIND`

When swicthed on, the Fortran binding for LibRPA will be built.

(enable-driver)=
## `ENABLE_DRIVER`

When swicthed on, the driver executable of LibRPA will be built.

(enable-test)=
## `ENABLE_TEST`

When swicthed on, the unit tests of LibRPA will be built.
After successful compile of LibRPA, one can issue `make test` under the build directory
to perform the unit tests.

```{note}
At present the unit tests do not cover the whole code base.
We are still working on it.
```

(libri-include-dir)=
## `LIBRI_INCLUDE_DIR`

The path to the include directory of LibRI.
When empty, the internal LibRI will be used.
Otherwise it will search for `RI/ri/RI_Tools.h` under the specified directory.
Error will be raised if the search fails.

(libcomm-include-dir)=
## `LIBCOMM_INCLUDE_DIR`

The path to the include directory of LibComm.
When empty, the internal LibComm will be used.
Otherwise it will search for `Comm/Comm_Tools.h` under the specified directory.
Error will be raised if the search fails.
