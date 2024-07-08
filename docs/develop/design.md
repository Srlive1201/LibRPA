# Code Design

LibRPA is designed to provide efficient and scalable implementations of various beyond-DFT methods, including but not limited to random-phase approximation (RPA) correlation energies. The library is developed in C++ and Python, utilizing a hybrid parallelism model with MPI and OpenMP to maximize performance on modern high-performance computing systems.

## Framework of RPA calculation

![Framework of RPA calculation](../assets/librpa_flowchat.jpg)

## RPA interface to FHI-aims

![RPA interface to FHI-aims](../assets/aims_interface.drawio.jpg)

## RPA interface to ABACUS

![RPA interface to ABACUS](../assets/abacus_interface.drawio.jpg)

## TODO

- [ ] Adapt RPA force work by Mohammad in the [backup branch](https://github.com/Srlive1201/LibRPA/tree/master-backup-240416)
