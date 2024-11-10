# API Usage

The user APIs that can be called from within DFT code are defined in {librpa}`librpa.h`.
There are three groups of API functions

- Environment setup
- Input parsing
- Computation

## Environment setup

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

## Input parsing

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

## Computation
After data have been correctly set up by input parsing functions, computation functions
can be called to get quantities of interest.
For example, the RPA correlation energy, including the contribution from each irreducible k-point, can be obtained
by calling {librpa}`get_rpa_correlation_energy`.
