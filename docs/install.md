# Installation

## Dependencies

LibRPA is built on top of several external libraries:

* [LibRI](https://github.com/abacusmodeling/LibRI) for tensor contraction
* [LibComm](https://github.com/abacusmodeling/LibComm)
  as a dependency of LibRI for communication of tensor data between processes
* [cereal](https://uscilab.github.io/cereal)
  as a dependency of LibRI for data serialization.
* Minimax time-frequency grids: The original grids are obtained from the CP2K
  code, but now they are available as a component of the [GreenX](https://github.com/nomad-coe/greenX)library.
* The Python package `scipy` is required, if the original CP2K Minimax grids are
  used, to calculate the transformation matrix.

The dependencies can be found under the `thirdparty/` directory. Each
dependency is included either as a submodule (LibRI and LibComm), or by storing
the code of a release version (possibly with minor modification). One exception
is the original Minimax grids from CP2K: the grids are stored as plain text
files in the `src/minimax_grid/GreenX` directory.

## Download

You can obtain the LibRPA code by cloning the GitHub repository:

```bash
git clone https://github.com/Srlive1201/LibRPA
```

Then go to the cloned repository and download dependencies which are submodules
of the project

```bash
cd LibRPA && git submodule update --init --recursive
```

When compiling LibRPA without linking the GreenX library, the `scipy` package
is required. Its installation can be checked by running

```bash
python -c "import scipy"
```

If you get a non-zero return code, you may need to install it by following the
[official documentation](https://scipy.org/install).

That's all and you are ready to compile.

## Compile

To compile LibRPA, you need a C++ compiler supporting MPI and a ScaLAPACK library.
The Intel MPI compiler and MKL library from Intel oneAPI tools (both base and
hpc toolkits) appear to be the most straightforward choice. Alternatively, you
can use GNU Compiler Collection (GCC) along with open source MPI
implementation (e.g. [MPICH](https://www.mpich.org)) and ScaLAPACK library ([Netlib](https://www.netlib.org/scalapack>)).
A Fortran compiler is also needed if you want to use the GreenX library API or build the Fortran binding.

CMake build system is used to build LibRPA library and driver executable.
Under the root directory of LibRPA, run the following commands

```bash
mkdir build
cd build
cmake -DUSE_LIBRI=ON ..
make -j 4
```

This will build both the library and drivers of LibRPA.
When the build process finished, you can find the driver `librpa_driver.exe`
and the shared library `src/librpa.so` in the `build` directory.

You can specify the compiler by prefixing the cmake command.
For example, to use the Intel C++ classic compiler

```bash
CXX=mpiicpc cmake -DUSE_LIBRI=ON ..
```

LibRPA compiled from the above commands will use the original Minimax grids
from CP2K during calculations. To use the updated version of Minimax grids in
the GreenX library, the following cmake command should be used

```bash
cmake -DUSE_LIBRI=ON -DUSE_GREENX_API=ON ..
```

Note that for CMake to find the correct ScaLAPACK libraries for linking, you
may need to add the directory of the libraries to `LIBRARY_PATH` or
`LD_LIBRARY_PATH` before the cmake command.
For example, to use the MKL libraries

```bash
export LD_LIBRARY_PATH="$MKLROOT/lib/intel64:$LD_LIBRARY_PATH"
CXX=mpiicpc FC=ifort cmake -DUSE_LIBRI=ON -DUSE_GREENX_API=ON ..
```

For a comprehensive list of compile options, please refer to the [user guide](user_guide/compile_options).
