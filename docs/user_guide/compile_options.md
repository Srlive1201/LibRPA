# Compile Options

## Overview

| Option                                                      | Type              | Default    |
|-------------------------------------------------------------|-------------------|------------|
| [`LIBRPA_USE_LIBRI`](#librpa-use-libri)                     | bool              | `OFF`      |
| [`LIBRPA_USE_CMAKE_INC`](#librpa-use-cmake-inc)             | bool              | `OFF`      |
| [`LIBRPA_USE_EXTERNAL_GREENX`](#librpa-use-external-greenx) | bool              | `OFF`      |
| [`LIBRPA_ENABLE_FORTRAN_BIND`](#librpa-enable-fortran-bind) | bool              | `OFF`      |
| [`LIBRPA_FORTRAN_DP`](#librpa-fortran-dp)                   | string or integer | `c_double` |
| [`LIBRPA_ENABLE_DRIVER`](#librpa-enable-driver)             | bool              | `ON`       |
| [`LIBRPA_ENABLE_TEST`](#librpa-enable-test)                 | bool              | `ON`       |
| [`LIBRPA_ENABLE_CPP_TEST`](#librpa-enable-cpp-test)         | bool              | `ON`       |
| [`LIBRPA_ENABLE_FORTRAN_TEST`](#librpa-enable-fortran-test) | bool              | `ON`       |
| [`LIBRI_INCLUDE_DIR`](#libri-include-dir)                   | string            | empty      |
| [`LIBCOMM_INCLUDE_DIR`](#libcomm-include-dir)               | string            | empty      |
| [`CEREAL_INCLUDE_DIR`](#cereal-include-dir)                 | string            | empty      |
| [`SCALAPACK_DIR`](#scalapack-dir)                           | string            | empty      |

These options can be parsed on the CMake command line, for example:

```sh
cmake -DLIBRPA_USE_LIBRI=ON
```

(librpa-use-libri)=
## `LIBRPA_USE_LIBRI`

When enabled, LibRPA is compiled with [LibRI](https://github.com/abacusmodeling/LibRI)
for RI tensor contractions.

The *GW* and EXX functionalities require LibRPA to be compiled with LibRI, i.e. `-DLIBRPA_USE_LIBRI=ON`.
By contrast, the RPA correlation energy can also be computed without this option.

(librpa-use-cmake-inc)=
## `LIBRPA_USE_CMAKE_INC`

When enabled, the `cmake.inc` file is used to initialize compilers and other build options.

**Deprecated**. It is recommended to use standard CMake command-line options such as `-C` or `-D` to specify custom variables.

(librpa-use-external-greenx)=
## `LIBRPA_USE_EXTERNAL_GREENX`

Controls whether LibRPA uses the bundled GreenX library or an external one.

The minimax grids used by LibRPA are provided through the
[GreenX](https://nomad-coe.github.io/greenX/) library.

When this option is `OFF` (default), LibRPA builds and links against the bundled GreenX source distributed with LibRPA under `thirdparty/greenX`.

When this option is `ON`, LibRPA does not build the bundled GreenX copy.
Instead, it expects an external GreenX library to be provided by the parent or higher-level CMake project.
In particular, the CMake target `LibGXMiniMax` must already be defined and available for linking.

This option is mainly intended for developer workflows or project setups in which GreenX is managed outside LibRPA.

(librpa-enable-fortran-bind)=
## `LIBRPA_ENABLE_FORTRAN_BIND`

When enabled, the Fortran bindings of LibRPA are built.

(librpa-fortran-dp)=
## `LIBRPA_FORTRAN_DP`

Specifies the Fortran kind used for double-precision real and complex data in the Fortran bindings.

The default value is `c_double`, which is suitable when interoperability with C is desired.
This option may also be set to an integer kind value if needed by the calling code.

This option is meaningful only if `LIBRPA_ENABLE_FORTRAN_BIND=ON`.

(librpa-enable-driver)=
## `LIBRPA_ENABLE_DRIVER`

When enabled, the LibRPA driver executable is built.

(librpa-enable-test)=
## `LIBRPA_ENABLE_TEST`

When enabled, the unit tests of LibRPA are built.

After LibRPA has been compiled successfully, the tests can be run from the build directory with:
```sh
ctest
```
or equivalently
```sh
make test
```

```{note}
At present, the unit tests do not cover the entire code base.
Test coverage is still being expanded.
```

(librpa-enable-cpp-test)=
## `LIBRPA_ENABLE_CPP_TEST`

When enabled, the C++ unit tests are built.

This option is meaningful only if `LIBRPA_ENABLE_TEST=ON`.

(librpa-enable-fortran-test)=
## `LIBRPA_ENABLE_FORTRAN_TEST`

When enabled, the Fortran unit tests are built.

This option is meaningful only if both `LIBRPA_ENABLE_TEST=ON` and `LIBRPA_ENABLE_FORTRAN_BIND=ON`.

(libri-include-dir)=
## `LIBRI_INCLUDE_DIR`

Specifies the path to the LibRI include directory.

If this variable is empty, the internal LibRI copy is used.
Otherwise, CMake searches for `RI/ri/RI_Tools.h` under the specified directory.
An error is raised if the file cannot be found.

Example:
```sh
cmake -DLIBRI_INCLUDE_DIR=/path/to/LibRI/include
```

(libcomm-include-dir)=
## `LIBCOMM_INCLUDE_DIR`

Specifies the path to the LibComm include directory.

If this variable is empty, the internal LibComm copy is used.
Otherwise, CMake searches for `Comm/Comm_Tools.h` under the specified directory.
An error is raised if the file cannot be found.

Example:
```sh
cmake -DLIBCOMM_INCLUDE_DIR=/path/to/LibComm/include
```

(cereal-include-dir)=
## `CEREAL_INCLUDE_DIR`

Specifies the path to the cereal include directory.

If this variable is empty, the bundled cereal copy is used.
Otherwise, CMake searches for `cereal/cereal.hpp` under the specified directory.
An error is raised if the file cannot be found.

Example:
```sh
cmake -DCEREAL_INCLUDE_DIR=/path/to/cereal/include
```

(scalapack-dir)=
## `SCALAPACK_DIR`

`SCALAPACK_DIR` specifies the installation path of ScaLAPACK and is used to
locate the ScaLAPACK libraries when `MKLROOT` is not defined.

This variable can be provided in two ways:

- as a CMake option:

  ```bash
  cmake -DSCALAPACK_DIR=/path/to/scalapack
  ```

- or as an environment variable:

  ```bash
  export SCALAPACK_DIR=/path/to/scalapack
  cmake
  ```

This option is intended for environments where ScaLAPACK is provided as a
standalone installation rather than through Intel MKL.
