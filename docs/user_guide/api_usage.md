# API Usage

## Summary of user API functions

The user-facing APIs that can be called from within a DFT code are declared in
{librpa}`librpa.h` (C), {librpa}`librpa.hpp` (C++), and {librpa}`librpa_f03.f90` (Fortran 2003).

These APIs fall into three groups:

- global environment management (initialization and finalization)
- input parsing
- computation

They are centered around two management structures:

- {librpa}`LibrpaHandler`, which manages input parsing and computations
- {librpa}`LibrpaOptions`, which stores runtime parameters

## Handler object

The {librpa}`LibrpaHandler` (C) / {librpa}`librpa::Handler` (C++) / `type(LibrpaHandler)` (Fortran) is the central object that manages the LibRPA computation. It encapsulates all input data and computation state.

The typical workflow is:
1. Create a handler using the MPI communicator
2. Set input data via setter methods
3. Call computation functions
4. Free the handler when done

**C:**
```c
LibrpaHandler *h = librpa_create_handler(MPI_COMM_WORLD);
// ... set input data and compute ...
librpa_destroy_handler(h);
```

**C++:**
```cpp
librpa::Handler h(MPI_COMM_WORLD);
// ... set input data and compute ...
h.free();
```

**Fortran:**
```fortran
type(LibrpaHandler) :: h
call h%init(MPI_COMM_WORLD)
! ... set input data and compute ...
call h%free()
```

(environment-setup)=
## Environment setup

Runtime environment for LibRPA should be properly setup to parse the input data and call
computation functions. They are managed by the environment setup APIs (C/C++/Fortran):

- {librpa}`librpa_init_global` / {librpa}`librpa::init_global`
- {librpa}`librpa_finalize_global` / {librpa}`librpa::finalize_global`
- {librpa}`librpa_init_options` / {librpa}`librpa::Options`
- {librpa}`librpa_set_output_dir` / {librpa}`librpa::Options::set_output_dir`

Any calls to input parsing and computation APIs should be surrounded by calls to `librpa_init_global` and `librpa_finalize_global`.

To configure runtime parameters for computation, an {librpa}`LibrpaOptions` object should be created, modified and passed to computation functions.
Parameters are attributes of `LibrpaOptions` object, and can be initialized using `librpa_init_options`.
For example, to configure the number of frequency grids and Green's function threshold:
```c
// C API
LibrpaOptions opts;
librpa_init_options(&opts);
opts.nfreq = 16;
opts.gf_threshold = 1.e-3;
opts.output_level = LIBRPA_VERBOSE_INFO;
```

```cpp
// C++ API
librpa::Options opts;
opts.nfreq = 16;
opts.gf_threshold = 1.e-3;
opts.output_level = librpa::Verbose::LIBRPA_VERBOSE_INFO;
```

```fortran
! Fortran API
type(LibrpaOptions) :: opts
call opts%init()
opts%nfreq = 16
opts%gf_threshold = 1.0d-3
opts%output_level = LIBRPA_VERBOSE_INFO
```

## Input parsing

Input parsing APIs common to all tasks include (C/C++/Fortran):

- {librpa}`librpa_set_scf_dimension` / {librpa}`Handler::set_scf_dimension` / {librpa}`handler%set_scf_dimension <librpahandler::set_scf_dimension>`
- {librpa}`librpa_set_wg_ekb_efermi` / {librpa}`Handler::set_wg_ekb_efermi` / {librpa}`handler%set_wg_ekb_efermi <librpahandler::set_wg_ekb_efermi>`
- {librpa}`librpa_set_wfc` / {librpa}`Handler::set_wfc` / {librpa}`handler%set_wfc <librpahandler::set_wfc>`
- {librpa}`librpa_set_ao_basis_wfc` / {librpa}`Handler::set_ao_basis_wfc` / {librpa}`handler%set_ao_basis_wfc <librpahandler::set_ao_basis_wfc>`
- {librpa}`librpa_set_ao_basis_aux` / {librpa}`Handler::set_ao_basis_aux` / {librpa}`handler%set_ao_basis_aux <librpahandler::set_ao_basis_aux>`
- {librpa}`librpa_set_latvec_and_G` / {librpa}`Handler::set_latvec_and_G` / {librpa}`handler%set_latvec_and_G <librpahandler::set_latvec_and_g>`
- {librpa}`librpa_set_atoms` / {librpa}`Handler::set_atoms` / {librpa}`handler%set_atoms <librpahandler::set_atoms>`
- {librpa}`librpa_set_kgrids_kvec` / {librpa}`Handler::set_kgrids_kvec` / {librpa}`handler%set_kgrids_kvec <librpahandler::set_kgrids_kvec>`
- {librpa}`librpa_set_ibz_mapping` / {librpa}`Handler::set_ibz_mapping` / {librpa}`handler%set_ibz_mapping <librpahandler::set_ibz_mapping>`
- {librpa}`librpa_set_lri_coeff` / {librpa}`Handler::set_lri_coeff` / {librpa}`handler%set_lri_coeff <librpahandler::set_lri_coeff>`
- {librpa}`librpa_set_aux_bare_coulomb_k_atom_pair` / {librpa}`Handler::set_aux_bare_coulomb_k_atom_pair` / {librpa}`handler%set_aux_bare_coulomb_k_atom_pair <librpahandler::set_aux_bare_coulomb_k_atom_pair>`
- {librpa}`librpa_set_aux_cut_coulomb_k_atom_pair` / {librpa}`Handler::set_aux_cut_coulomb_k_atom_pair` / {librpa}`handler%set_aux_cut_coulomb_k_atom_pair <librpahandler::set_aux_cut_coulomb_k_atom_pair>`
- {librpa}`librpa_set_aux_bare_coulomb_k_2d_block` / {librpa}`Handler::set_aux_bare_coulomb_k_2d_block` / {librpa}`handler%set_aux_bare_coulomb_k_2d_block <librpahandler::set_aux_bare_coulomb_k_2d_block>`
- {librpa}`librpa_set_aux_cut_coulomb_k_2d_block` / {librpa}`Handler::set_aux_cut_coulomb_k_2d_block` / {librpa}`handler%set_aux_cut_coulomb_k_2d_block <librpahandler::set_aux_cut_coulomb_k_2d_block>`
- {librpa}`librpa_set_band_kvec` / {librpa}`Handler::set_band_kvec` / {librpa}`handler%set_band_kvec <librpahandler::set_band_kvec>`
- {librpa}`librpa_set_band_occ_eigval` / {librpa}`Handler::set_band_occ_eigval` / {librpa}`handler%set_band_occ_eigval <librpahandler::set_band_occ_eigval>`
- {librpa}`librpa_set_wfc_band` / {librpa}`Handler::set_wfc_band` / {librpa}`handler%set_wfc_band <librpahandler::set_wfc_band>`

These functions parse information from the hosting DFT code to LibRPA, such as:
- dimension of lattice vectors, orbital basis and auxiliary basis
- k-mesh points, eigenvalues and occupation numbers
- wave functions of mean-field calculation
- local RI coefficients
- Coulomb matrices in auxiliary basis
- data for band structure calculation

They have to be called after the runtime environment is initialized,
see [Environment setup](#environment-setup) above.

## Computation

After data have been correctly set up by input parsing functions, computation functions can be called to get quantities of interest:

**General**
- {librpa}`librpa_get_imaginary_frequency_grids`
  / {librpa}`Handler::get_imaginary_frequency_grids`
  / {librpa}`handler%get_imaginary_frequency_grids <librpahandler::get_imaginary_frequency_grids>`

**RPA correlation energy:**
- {librpa}`librpa_get_rpa_correlation_energy`
  / {librpa}`Handler::get_rpa_correlation_energy`
  / {librpa}`handler%get_rpa_correlation_energy <librpahandler::get_rpa_correlation_energy>`

**Exact exchange:**
- {librpa}`librpa_build_exx` / {librpa}`Handler::build_exx` / {librpa}`handler%build_exx <librpahandler::build_exx>`
- {librpa}`librpa_get_exx_pot_kgrid` / {librpa}`Handler::get_exx_pot_kgrid` / {librpa}`handler%get_exx_pot_kgrid <librpahandler::get_exx_pot_kgrid>`
- {librpa}`librpa_get_exx_pot_band_k` / {librpa}`Handler::get_exx_pot_band_k` / {librpa}`handler%get_exx_pot_band_k <librpahandler::get_exx_pot_band_k>`

**G0W0 self-energy:**
- {librpa}`librpa_build_g0w0_sigma` / {librpa}`Handler::build_g0w0_sigma` / {librpa}`handler%build_g0w0_sigma <librpahandler::build_g0w0_sigma>`
- {librpa}`librpa_get_g0w0_sigc_kgrid` / {librpa}`Handler::get_g0w0_sigc_kgrid` / {librpa}`handler%get_g0w0_sigc_kgrid <librpahandler::get_g0w0_sigc_kgrid>`
- {librpa}`librpa_get_g0w0_sigc_band_k` / {librpa}`Handler::get_g0w0_sigc_band_k` / {librpa}`handler%get_g0w0_sigc_band_k <librpahandler::get_g0w0_sigc_band_k>`

## Schematic examples

### C++

The following C++ code gives an impression of how LibRPA APIs should be called in the host DFT code:

```cpp
#include <librpa.hpp>
#include <mpi.h>

int main()
{
    MPI_Init(...);
    // If LibRI is enabled, MPI thread support should be set to MPI_THREAD_MULTIPLE
    // MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

    // Before: setup and SCF calculation of hosting program

    // *** Initialize LibRPA runtime environment ***
    librpa::init_global();

    // *** Set up runtime parameters ***
    librpa::Options opts;
    opts.nfreq = 16;
    opts.gf_threshold = 1e-3;
    opts.output_level = librpa::Verbose::LIBRPA_VERBOSE_INFO;

    // *** Create handler ***
    librpa::Handler h(MPI_COMM_WORLD);

    // *** Input parsing ***
    h.set_scf_dimension(nspins, nkpts, nstates, nbasis);
    h.set_wg_ekb_efermi(nspins, nkpts, nstates, wg, ekb, efermi);
    h.set_latvec_and_G(lat_mat, G_mat);
    h.set_kgrids_kvec(nk1, nk2, nk3, kvecs);
    h.set_ibz_mapping(ibz_map);
    h.set_ao_basis_wfc(nbs_wfc);
    h.set_ao_basis_aux(nbs_aux);

    // Loop over atoms and R cells to set LRI coefficients
    for (int i = 0; i < natoms; ++i) {
        for (int j = 0; j < natoms; ++j) {
            for (auto& R : r_cells) {
                h.set_lri_coeff(LIBRPA_ROUTING_AUTO, i, j, nbasis_i, nbasis_j, naux_mu, R.data(), Cs.data());
            }
        }
    }

    // Loop over k-points to set Coulomb matrices
    for (int ik = 0; ik < nkpts; ++ik) {
        h.set_aux_bare_coulomb_k_atom_pair(ik, i, j, naux_i, naux_j, Vq_real.data(), Vq_imag.data(), vq_threshold);
    }

    // *** Compute RPA correlation energy ***
    std::vector<std::complex<double>> rpa_contrib_ibzk(nkpts_ibz);
    double Ec = h.get_rpa_correlation_energy(opts, rpa_contrib_ibzk);

    // *** Compute exact exchange ***
    h.build_exx(opts);
    std::vector<double> vexx = h.get_exx_pot_kgrid(opts, nspins, iks_this, i_state_low, i_state_high);

    // *** Compute G0W0 self-energy ***
    h.build_g0w0_sigma(opts);
    std::vector<std::complex<double>> sigc = h.get_g0w0_sigc_kgrid(opts, nspins, iks_this, i_state_low, i_state_high, vxc, vexx);

    // *** Clean up ***
    h.free();
    librpa::finalize_global();

    MPI_Finalize();
}
```

### Fortran

```fortran
use librpa_f03
use mpi

implicit none

type(LibrpaOptions) :: opts
type(LibrpaHandler) :: h
integer :: ierr
integer :: comm, nkpts_ibz
complex(dp), allocatable :: contrib_ibzk(:)

call MPI_Init(ierr)
comm = MPI_COMM_WORLD

! Initialize LibRPA
call librpa_init_global()

! Set up options
call opts%init()
opts%nfreq = 16
opts%gf_threshold = 1.0d-3

! Initialize handler
call h%init(comm)

! Input parsing
call h%set_scf_dimension(nspins, nkpts, nstates, nbasis)
call h%set_wg_ekb_efermi(nspins, nkpts, nstates, wg, ekb, efermi)
call h%set_latvec_and_G(lat, recplatt)
! ... more input calls ...

! Compute RPA correlation energy
allocate(contrib_ibzk(nkpts_ibz))
Ec = h%get_rpa_correlation_energy(opts, nkpts_ibz, contrib_ibzk)

! Clean up
call h%free()

! Finalize LibRPA
call librpa_finalize_global()

call MPI_Finalize(ierr)
```
