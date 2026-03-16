# Quasi-particle band structure using LibRPA driver with FHI-aims dataset

This tutorial guides you through the process of setting up and running a calculation using FHI-aims for geometry and SCF computations and LibRPA for quasi-particle energies by one-shot GW method.

The silicon unit cell is used as test case.

## 1. **Prerequisites**
Before starting, ensure that you have:
- [FHI-aims](https://fhi-aims.org/get-the-code-menu/get-the-code) installed with `USE_GREENX` set to `ON`.
- LibRPA [installed](../../install.md) with LibRI enabled.

## 2. **FHI-aims Input Files**
These files are required to run the initial SCF and geometry optimization using FHI-aims:
- **`control.in`**: This file contains control parameters for FHI-aims. Below is an example:
  ```text
    # Basic model
    xc               pbe
    k_grid           4 4 4
    occupation_type  gaussian 0.001

    # GW switches
    qpe_calc         gw_expt
    freq_grid_type   minimax
    frequency_points 16
    anacon_type      1

    # Output flags
    output librpa binary develop fold_C
    output band   0.50000  0.50000  0.50000   0.00000  0.00000  0.00000 13 L G
    output band   0.00000  0.00000  0.00000   0.50000  0.00000  0.50000 13 G X

    [light species default for Si]
  ```
- **`geometry.in`**: Contains the geometry of the system. For silicon:
  ```text
    lattice_vector   3.8301668167   0.0000000000   0.0000000000
    lattice_vector   1.9150834084   3.3170217640   0.0000000000
    lattice_vector   1.9150834084   1.1056739213   3.1273181102
    atom_frac        0.0000000000   0.0000000000   0.0000000000  Si
    atom_frac        0.2500000000   0.2500000000   0.2500000000  Si
  ```

Then run FHI-aims to generate dataset files.

## 3. **LibRPA Input File**

- **`librpa.in`**: This file specifies the parameters for the GW calculation. Example:
  ```text
    task = g0w0_band
    nfreq = 16
    option_dielect_func = 0
    replace_w_head = t
    use_scalapack_gw_wc = t
    parallel_routing = libri
  ```

  Here `replace_w_head` is switched on, so that LibRPA uses the dielectric function
  directly computed in FHI-aims for the correction of dielectric matrix to speed up
  k-point convergence.

## 4. **Run the GW Calculation with LibRPA**

After obtaining the output files from FHI-aims and setting up the parameters in `librpa.in`, you can use LibRPA to calculate the quasi-particle band structure:
```bash
mpirun -np 4 /path/to/LibRPA/build/chi0_main.exe  >  LibRPA.out
```

## 5. **LibRPA Output**

After successful run, LibRPA writes the band structures to the corresponding files:
- `KS_band_spin_<ispin>.dat`: input Kohn-Sham band structure, for reference purpose
- `EXX_band_spin_<ispin>.dat`: non-self-consistent exact-exchange-only band structure
- `GW_band_spin_<ispin>.dat`: GW quasi-particle band structure

All band structure files share the same format as those in FHI-aims, and hence
can be post-processed similarly using existing scripts for FHI-aims.
