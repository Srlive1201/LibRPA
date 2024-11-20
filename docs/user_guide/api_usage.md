# API Usage

## Summary of user API functions

The user APIs that can be called from within DFT code are defined in {librpa}`librpa.h`.
There are three groups of API functions

- Environment setup
- Input parsing
- Computation

### Environment setup

Runtime environment for LibRPA should be properly setup to parse the input data and call
computation function. They are managed by the environment setup APIs

- {librpa}`initialize_librpa_environment`
- {librpa}`finalize_librpa_environment`
- {librpa}`get_default_librpa_params`
- {librpa}`set_librpa_params`

Any calls to input parsing and computation APIs should be surrounded by calls to `initialize_librpa_environment` and `finalize_librpa_environment`.
To configure runtime parameters for input parsing and computation, an {librpa}`LibRPAParams` object should be created, modified and parsed to `set_librpa_params`.
Parameters are attributes of `LibRPAParams` object, and can be initialized using `get_default_librpa_params`.
For example, to configure the number of frequency grids and Green's function threshold,
one can add the following code after the runtime environment is initialized
```c
LibRPAParams params;
get_default_librpa_params(&params);
params.nfreq = 16;
params.gf_R_threshold = 1.e-3;
set_librpa_params(&params);
```

### Input parsing

Input parsing APIs common to all tasks include
- {librpa}`set_dimension`
- {librpa}`set_wg_ekb_efermi`
- {librpa}`set_ao_basis_wfc`
- {librpa}`set_latvec_and_G`
- {librpa}`set_kgrids_kvec_tot`
- {librpa}`set_ibz2bz_index_and_weight`
- {librpa}`set_ao_basis_aux`
- {librpa}`set_aux_bare_coulomb_k_atom_pair`
- {librpa}`set_aux_bare_coulomb_k_2D_block`

These functions parse information from the hosting DFT code to LibRPA, such as
dimension of lattice vectors, orbital basis and auxiliary basis, k-mesh points, eigenvalues and
wave functions of mean-field calculation.
They have to be called after the runtime environment is setup.

### Computation

After data have been correctly set up by input parsing functions, computation functions
can be called to get quantities of interest.
For example, the RPA correlation energy, including the contribution from each irreducible k-point, can be obtained
by calling {librpa}`get_rpa_correlation_energy`.

## Schematic example

The following C++ code gives an impression of how LibRPA APIs should be called
in the host DFT code.

```cpp
#include <librpa.h>

int main()
{
    // Before: setup and SCF calculation of hosting program

    // *** Initialize LibRPA runtime environment ***
    initialize_librpa_environment(MPI_COMM_WORLD, 0, 1, "LibRPA.out");

    // *** Set up runtime parameters ***
    LibRPAParams params_librpa;
    get_default_librpa_params(&params_librpa);
    params_librpa.nfreq = 16;
    set_librpa_params(&params_librpa);

    // *** Input parsing ***
    set_dimension(nspins, nkpts, nstates, nbasis, natoms);
    set_wg_ekb_efermi(nspins, nkpts, nstates, wg, ekb, efermi);
    set_latvec_and_G(lat_mat, G_mat);
    set_kgrids_kvec_tot(nk1, nk2, nk3, kvecs);
    set_ibz2bz_index_and_weight(nk_irk, ibz2bz_index, wk_irk);

    // Loop over ia1, ia2 and R
    set_ao_basis_aux(ia1, ia2, n_basis_1, n_basis_2, n_auxbas_1, R, lri_c, 1);
    // Loop over ia1, ia2 and iq
    set_aux_cut_coulomb_k_atom_pair(iq, ia1, ia2, n_auxbas_1, n_auxbas_2, coul_real, coul_imag)

    // *** Compute RPA correlation energy ***
    double rpa_corr, rpa_corr_irk_contrib[nk_irk];
    get_rpa_correlation_energy(&rpa_corr, rpa_corr_irk_contrib);

    // *** Clean up and finalize LibRPA runtime ***
    finalize_librpa_environment();

    // After: post-processing in hosting program and other tasks
}
```
