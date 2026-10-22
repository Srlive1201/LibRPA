# Code Design

LibRPA is designed to provide efficient and scalable implementations of various beyond-DFT methods, including but not limited to random-phase approximation (RPA) correlation energies. The library is developed in C++ and Python, utilizing a hybrid parallelism model with MPI and OpenMP to maximize performance on modern high-performance computing systems.

## Framework of RPA calculation

![Framework of RPA calculation](../assets/librpa_flowchat.jpg)

## RPA interface to FHI-aims

![RPA interface to FHI-aims](../assets/aims_interface.drawio.jpg)

## RPA interface to ABACUS

![RPA interface to ABACUS](../assets/abacus_interface.drawio.jpg)

## QSGW task boundary

The new QSGW iteration, restart state, live optical updates and diagnostic results
are implemented in `src/core/qsgw/`. The public functions in `librpa_qsgw.h`
marshal numerical inputs and output callbacks; the driver reads host files and
writes the returned output events.

The existing `src/api/Dataset` ownership and G0W0 workflow remain in place.
QSGW borrows only the numerical objects it uses through its own `TaskContext`;
the API task retains the original dataset for the borrowed references' lifetime.
Initialization and the existing G0W0 Sigma build are supplied by API adapters.
The optional response observer runs at the resident full-response boundary so a
companion task can consume that response before screening releases it. Core QSGW
code includes no API or driver header. Its entering-state QP diagnostic and table
formatters belong to QSGW and preserve the existing continuation options and
output conventions.
