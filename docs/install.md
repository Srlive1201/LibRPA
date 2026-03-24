# Installation

## Dependencies

LibRPA depends on the following core software components:

- a C++ compiler and an MPI library
- a Fortran compiler (and MPI support if [`LIBRPA_ENABLE_FORTRAN_BIND`](<librpa-enable-fortran-bind>) is enabled)
- BLAS and LAPACK libraries
- a ScaLAPACK library
- the [GreenX](https://github.com/nomad-coe/greenX) library for minimax
  time-frequency grids

For *GW*, the following packages are additionally required:

- [LibRI](https://github.com/abacusmodeling/LibRI) for tensor contractions
- [LibComm](https://github.com/abacusmodeling/LibComm), which is required by LibRI for communication of tensor data between processes
- [cereal](https://uscilab.github.io/cereal), which is required by LibRI for data serialization

Some of these dependencies are located under the `thirdparty/` directory.
Depending on the package, they are included either as Git submodules or as
bundled source code distributed with LibRPA.

## Download

You can obtain the LibRPA code by cloning the GitHub repository:

```bash
git clone https://github.com/Srlive1201/LibRPA
```

Then enter the repository and initialize the submodules:
```bash
cd LibRPA
git submodule update --init --recursive
```

After this step, the source tree is ready for compilation.

## Compile

To compile LibRPA, you need working compiler and library toolchains for C++,
Fortran, MPI, BLAS/LAPACK, and ScaLAPACK.

The Intel compilers and Intel MPI together with MKL from Intel oneAPI are often
the most straightforward choice. Alternatively, LibRPA can also be built with
GCC/GFortran together with an open-source MPI implementation such as
[MPICH](https://www.mpich.org) and an open-source ScaLAPACK library [Netlib ScaLAPACK](https://www.netlib.org/scalapack>).

LibRPA uses CMake as its build system.
Ensure that the compilers and required libraries can be found under directories specified
by relevant environment variables, and under the root directory of the source tree:
```bash
mkdir build
cd build
cmake ..
make -j 4
```
This searches for the paths of required dependencies, and builds the LibRPA library and the driver executable.

After a successful build, the driver executable `chi0_main.exe` and the shared
library `src/librpa.so` (`src/librpa.dylib` on macOS, or `src/librpa.a` if `BUILD_SHARED_LIBS` is disabled)
can be found in the **build directory**.

You can specify the compilers through environment variables when invoking
CMake. For example, to use the Intel classic C++ and Fortran compilers:

```bash
CXX=mpiicpc FC=mpiifort cmake ..
```

To help CMake find the correct BLAS/LAPACK and ScaLAPACK libraries at link
time, you may need to ensure that the corresponding library directories are
visible through `LIBRARY_PATH` or `LD_LIBRARY_PATH`. For example, when using
MKL:

```bash
export LD_LIBRARY_PATH="$MKLROOT/lib/intel64:$LD_LIBRARY_PATH"
CXX=mpiicpc FC=mpiifort cmake ..
```

By default, LibRPA builds and links against the bundled GreenX source
distributed under thirdparty/greenX.
If you want to use an external GreenX instead, you should enable the CMake
option [`LIBRPA_USE_EXTERNAL_GREENX`](<librpa-use-external-greenx>):

```bash
cmake -DLIBRPA_USE_EXTERNAL_GREENX=ON ..
```

In this case, LibRPA does not build the bundled GreenX copy. Instead, the
parent or higher-level CMake project must provide the external GreenX target
`LibGXMiniMax`.

Several build scripts are provided on the [`Build Examples`](examples/build/index)
page to help users build LibRPA on different platforms and with different toolchains.
You may use them as starting points and adapt them to your local environment.

For a complete list of compile options, please refer to the
[Compile Options](user_guide/compile_options) page.
