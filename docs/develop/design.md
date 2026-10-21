# Code Design

LibRPA is designed to provide efficient and scalable implementations of various beyond-DFT methods, including but not limited to random-phase approximation (RPA) correlation energies. The library is developed in C++ and Python, utilizing a hybrid parallelism model with MPI and OpenMP to maximize performance on modern high-performance computing systems.

## Framework of RPA calculation

![Framework of RPA calculation](../assets/librpa_flowchat.jpg)

## RPA interface to FHI-aims

![RPA interface to FHI-aims](../assets/aims_interface.drawio.jpg)

## RPA interface to ABACUS

![RPA interface to ABACUS](../assets/abacus_interface.drawio.jpg)

## Numerical G0W0 task and the driver

The G0W0 numerical workflow lives in `src/core/g0w0_task`; the core `Dataset`
owns its mean fields, response, screening state and parallel layouts. The API
registry resolves public handlers to these states. Numerical tasks assemble
lower-level math, MPI and I/O services. Atomic-basis adapters and physics-specific
output assembly live in core; generic serializers accept value buffers, and
lower modules do not include core physics headers.

The G0W0 driver parses host input, invokes the public C/C++ APIs and writes
requested output. Ordinary input filenames remain in the driver. Numerical
velocity preparation and analytic head/wing construction belong to the library.
Read-only spectrum queries supply the driver with energies, occupations and
k-point metadata for presentation.
