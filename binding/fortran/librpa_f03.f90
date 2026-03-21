module librpa_f03

   use iso_c_binding, only: c_char, c_ptr, c_int, c_double, c_null_ptr, c_size_t
   implicit none

   private

   !=======================================================================
   ! Public types, constants, and functions
   !=======================================================================
   public :: LibrpaOptions
   public :: LibrpaHandler

   ! Double precision to connect to the user code. Control it through -DLIBRPA_FORTRAN_DP
   integer, parameter, public :: dp = ${LIBRPA_FORTRAN_DP}

   integer, parameter, public :: LIBRPA_MAX_STRLEN = 200

   integer, parameter, public :: LIBRPA_VERBOSE_DEBUG = 4
   integer, parameter, public :: LIBRPA_VERBOSE_WARN = 3
   integer, parameter, public :: LIBRPA_VERBOSE_INFO = 2
   integer, parameter, public :: LIBRPA_VERBOSE_CRITICAL = 1
   integer, parameter, public :: LIBRPA_VERBOSE_SILENT = 0

   integer, parameter :: LIBRPA_UNSET = -101
   integer, parameter :: LIBRPA_AUTO = -51

   ! Parallel routing
   integer, parameter, public :: LIBRPA_ROUTING_UNSET = LIBRPA_UNSET
   integer, parameter, public :: LIBRPA_ROUTING_AUTO = LIBRPA_AUTO
   integer, parameter, public :: LIBRPA_ROUTING_RTAU = 0
   integer, parameter, public :: LIBRPA_ROUTING_ATOMPAIR = 1
   integer, parameter, public :: LIBRPA_ROUTING_LIBRI = 2

   ! Time-frequency grid types
   integer, parameter, public :: LIBRPA_TFGRID_UNSET = LIBRPA_UNSET
   integer, parameter, public :: LIBRPA_TFGRID_GL = 0
   integer, parameter, public :: LIBRPA_TFGRID_GCI= 1
   integer, parameter, public :: LIBRPA_TFGRID_GCII = 2
   integer, parameter, public :: LIBRPA_TFGRID_MINIMAX = 3
   integer, parameter, public :: LIBRPA_TFGRID_EVENSPACED = 4
   integer, parameter, public :: LIBRPA_TFGRID_EVENSPACED_TF = 5

   public :: librpa_init_global
   public :: librpa_finalize_global
   public :: librpa_get_major_version
   public :: librpa_get_minor_version
   public :: librpa_get_patch_version
   public :: librpa_test
   public :: librpa_print_profile

   !=======================================================================
   ! Switch as bool type, defined as in include/librpa_enums.h
   integer(c_int), parameter :: LIBRPA_SWITCH_OFF = 0
   integer(c_int), parameter :: LIBRPA_SWITCH_ON = 1
   character(kind=c_char), allocatable, target, save :: redirect_path_buf(:)

   complex(dp), parameter :: CIMAG = (0.0d0, 1.0d0)

   !===== C-side options type =====
   ! Must have the same data layout as the struct defined in include/librpa_options.h
   type, bind(c) :: LibrpaOptions_c
      ! Common runtime control
      character(kind=c_char, len=1) :: output_dir(LIBRPA_MAX_STRLEN)
      integer(c_int) :: parallel_routing
      integer(c_int) :: output_level
      real(c_double) :: cs_threshold
      real(c_double) :: vq_threshold
      integer(c_int) :: use_soc
      integer(c_int) :: use_kpara_scf_eigvec
      integer(c_int) :: tfgrids_type
      integer(c_int) :: nfreq
      real(c_double) :: tfgrids_freq_min
      real(c_double) :: tfgrids_freq_interval
      real(c_double) :: tfgrids_freq_max
      real(c_double) :: tfgrids_time_min
      real(c_double) :: tfgrids_time_interval

      ! RPA specific
      real(c_double) :: gf_threshold
      integer(c_int) :: use_scalapack_ecrpa

      ! GW specific
      integer(c_int) :: n_params_anacon
      integer(c_int) :: use_scalapack_gw_wc
      integer(c_int) :: replace_w_head
      integer(c_int) :: option_dielect_func
      real(c_double) :: sqrt_coulomb_threshold
      real(c_double) :: libri_chi0_threshold_C
      real(c_double) :: libri_chi0_threshold_G
      real(c_double) :: libri_exx_threshold_C
      real(c_double) :: libri_exx_threshold_D
      real(c_double) :: libri_exx_threshold_V
      real(c_double) :: libri_g0w0_threshold_C
      real(c_double) :: libri_g0w0_threshold_G
      real(c_double) :: libri_g0w0_threshold_Wc

      ! Output controls
      integer(c_int) :: output_gw_sigc_mat
      integer(c_int) :: output_gw_sigc_mat_rt
      integer(c_int) :: output_gw_sigc_mat_rf
   end type LibrpaOptions_c

   ! High-level Fortran wrapper
   type :: LibrpaOptions
      type(LibrpaOptions_c), private :: opts_c

      character(len=LIBRPA_MAX_STRLEN) :: output_dir
      integer :: parallel_routing
      integer :: output_level
      real(dp) :: cs_threshold
      real(dp) :: vq_threshold
      logical :: use_soc
      logical :: use_kpara_scf_eigvec
      integer :: tfgrids_type
      integer :: nfreq
      real(dp) :: tfgrids_freq_min
      real(dp) :: tfgrids_freq_interval
      real(dp) :: tfgrids_freq_max
      real(dp) :: tfgrids_time_min
      real(dp) :: tfgrids_time_interval
      real(dp) :: gf_threshold
      logical :: use_scalapack_ecrpa
      integer :: n_params_anacon
      logical :: use_scalapack_gw_wc
      logical :: replace_w_head
      integer :: option_dielect_func
      real(dp) :: sqrt_coulomb_threshold
      real(dp) :: libri_chi0_threshold_C
      real(dp) :: libri_chi0_threshold_G
      real(dp) :: libri_exx_threshold_C
      real(dp) :: libri_exx_threshold_D
      real(dp) :: libri_exx_threshold_V
      real(dp) :: libri_g0w0_threshold_C
      real(dp) :: libri_g0w0_threshold_G
      real(dp) :: libri_g0w0_threshold_Wc
      logical :: output_gw_sigc_mat
      logical :: output_gw_sigc_mat_rt
      logical :: output_gw_sigc_mat_rf

      contains
         procedure :: init => librpa_init_options
         procedure :: set_output_dir => librpa_set_output_dir
   end type LibrpaOptions

   interface
      ! void librpa_init_options(LibrpaOptions *opts);
      subroutine librpa_init_options_c(opts_c) bind(c, name="librpa_init_options")
         import :: LibrpaOptions_c
         type(LibrpaOptions_c) :: opts_c
      end subroutine librpa_init_options_c

      subroutine librpa_set_output_dir_c(opts_c, output_dir) bind(c, name="librpa_set_output_dir")
         import :: LibrpaOptions_c, c_char
         type(LibrpaOptions_c) :: opts_c
         character(kind=c_char, len=1), dimension(*), intent(in) :: output_dir
      end subroutine librpa_set_output_dir_c
   end interface

   integer, parameter :: SYNC_OPTS_C2F = 1
   integer, parameter :: SYNC_OPTS_F2C = -1

   ! Global environment
   interface
      subroutine librpa_init_global_c(f_comm, sw_redirect, path, sw_process) &
            bind(c, name="librpa_init_global_fortran")
         import :: c_int, c_char
         integer(c_int) :: f_comm
         integer(c_int), value :: sw_redirect
         character(kind=c_char), dimension(*), intent(in) :: path
         integer(c_int), value :: sw_process
      end subroutine librpa_init_global_c

      subroutine librpa_finalize_global_c() bind(c, name="librpa_finalize_global")
      end subroutine librpa_finalize_global_c

      subroutine librpa_test_c() bind(c, name="librpa_test")
      end subroutine librpa_test_c

      subroutine librpa_print_profile_c() bind(c, name="librpa_print_profile")
      end subroutine librpa_print_profile_c
   end interface

   ! Version information
   interface
      function librpa_get_major_version_c() bind(c, name="librpa_get_major_version")
         import :: c_int
         integer(c_int) :: librpa_get_major_version_c
      end function librpa_get_major_version_c

      function librpa_get_minor_version_c() bind(c, name="librpa_get_minor_version")
         import :: c_int
         integer(c_int) :: librpa_get_minor_version_c
      end function librpa_get_minor_version_c

      function librpa_get_patch_version_c() bind(c, name="librpa_get_patch_version")
         import :: c_int
         integer(c_int) :: librpa_get_patch_version_c
      end function librpa_get_patch_version_c
   end interface

   ! High-level Fortran wrapper
   type :: LibrpaHandler
      type(c_ptr) :: ptr_c_handle = c_null_ptr
      contains
         ! Initialization and destruction
         procedure :: create  => librpa_create_handler
         procedure :: destroy => librpa_destroy_handler
         ! Input
         procedure :: set_scf_dimension => librpa_set_scf_dimension
         procedure :: set_wg_ekb_efermi => librpa_set_wg_ekb_efermi
         procedure :: set_wfc => librpa_set_wfc
         procedure :: set_ao_basis_wfc => librpa_set_ao_basis_wfc
         procedure :: set_ao_basis_aux => librpa_set_ao_basis_aux
         procedure :: set_latvec_and_G => librpa_set_latvec_and_G
         procedure :: set_atoms => librpa_set_atoms
         procedure :: set_kgrids_kvec => librpa_set_kgrids_kvec
         procedure :: set_ibz_mapping => librpa_set_ibz_mapping
         procedure :: set_lri_coeff => librpa_set_lri_coeff
         procedure :: set_aux_bare_coulomb_k_atom_pair => librpa_set_aux_bare_coulomb_k_atom_pair
         procedure :: set_aux_cut_coulomb_k_atom_pair => librpa_set_aux_cut_coulomb_k_atom_pair
         procedure :: set_aux_bare_coulomb_k_2d_block => librpa_set_aux_bare_coulomb_k_2d_block
         procedure :: set_aux_cut_coulomb_k_2d_block => librpa_set_aux_cut_coulomb_k_2d_block
         procedure :: set_dielect_func_imagfreq => librpa_set_dielect_func_imagfreq
         procedure :: set_band_kvec => librpa_set_band_kvec
         procedure :: set_wfc_band => librpa_set_wfc_band
         procedure :: set_band_occ_eigval => librpa_set_band_occ_eigval
         procedure :: reset_band_data => librpa_reset_band_data
         ! Compute
         procedure :: get_imaginary_frequency_grids => librpa_get_imaginary_frequency_grids
         procedure :: get_rpa_correlation_energy => librpa_get_rpa_correlation_energy
         procedure :: build_exx => librpa_build_exx
         procedure :: get_exx_pot_kgrid => librpa_get_exx_pot_kgrid
         procedure :: get_exx_pot_band_k => librpa_get_exx_pot_band_k
         procedure :: build_g0w0_sigma => librpa_build_g0w0_sigma
         procedure :: get_g0w0_qpe_kgrid => librpa_get_g0w0_qpe_kgrid
   end type LibrpaHandler

   interface
      function librpa_create_handler_c(f_comm) bind(c, name="librpa_create_handler_fortran")
         import :: c_ptr, c_int
         integer(c_int) :: f_comm
         type(c_ptr) :: librpa_create_handler_c
      end function librpa_create_handler_c

      subroutine librpa_destroy_handler_c(h) bind(c, name="librpa_destroy_handler")
         import :: c_ptr
         type(c_ptr), value :: h
      end subroutine librpa_destroy_handler_c
   end interface

   ! Input functions interface
   interface
      subroutine librpa_set_scf_dimension_c(h, nspins, nkpts, nstates, nbasis) &
                                            bind(c, name="librpa_set_scf_dimension")
         import :: c_ptr, c_int
         type(c_ptr), value :: h
         integer(c_int), value :: nspins, nkpts, nstates, nbasis
      end subroutine librpa_set_scf_dimension_c

      subroutine librpa_set_wg_ekb_efermi_c(h, nspins, nkpts, nstates, wg, ekb, efermi) bind(c, name="librpa_set_wg_ekb_efermi")
         import :: c_ptr, c_int, c_double
         type(c_ptr), value :: h
         integer(c_int), value :: nspins, nkpts, nstates
         real(c_double), dimension(*), intent(in) :: wg
         real(c_double), dimension(*), intent(in) :: ekb
         real(c_double), value :: efermi
      end subroutine librpa_set_wg_ekb_efermi_c

      subroutine librpa_set_wfc_c(h, ispin, ik, nstates_local, nbasis_local, wfc_real, wfc_imag) bind(c, name="librpa_set_wfc")
         import :: c_ptr, c_int, c_double
         type(c_ptr), value :: h
         integer(c_int), value :: ispin, ik, nstates_local, nbasis_local
         real(c_double), dimension(*), intent(in) :: wfc_real
         real(c_double), dimension(*), intent(in) :: wfc_imag
      end subroutine librpa_set_wfc_c

      subroutine librpa_set_wfc_packed_c(h, ispin, ik, nstates_local, nbasis_local, wfc) &
            bind(c, name="librpa_set_wfc_packed")
         import :: c_ptr, c_int
         type(c_ptr), value :: h
         integer(c_int), value :: ispin, ik, nstates_local, nbasis_local
         type(c_ptr), value :: wfc
      end subroutine librpa_set_wfc_packed_c

      subroutine librpa_set_ao_basis_wfc_c(h, natoms, nbs_wfc) bind(c, name="librpa_set_ao_basis_wfc")
         import :: c_ptr, c_int, c_size_t
         type(c_ptr), value :: h
         integer(c_int), value :: natoms
         integer(c_size_t), dimension(*), intent(in) :: nbs_wfc
      end subroutine librpa_set_ao_basis_wfc_c

      subroutine librpa_set_ao_basis_aux_c(h, natoms, nbs_aux) bind(c, name="librpa_set_ao_basis_aux")
         import :: c_ptr, c_int, c_size_t
         type(c_ptr), value :: h
         integer(c_int), value :: natoms
         integer(c_size_t), dimension(*), intent(in) :: nbs_aux
      end subroutine librpa_set_ao_basis_aux_c

      subroutine librpa_set_latvec_and_G_c(h, latt, recplatt) bind(c, name="librpa_set_latvec_and_G")
         import :: c_ptr, c_double
         type(c_ptr), value :: h
         real(c_double), dimension(9), intent(in) :: latt, recplatt
      end subroutine librpa_set_latvec_and_G_c

      subroutine librpa_set_atoms_c(h, natoms, types, posi_cart) bind(c, name="librpa_set_atoms")
         import :: c_ptr, c_int, c_double
         type(c_ptr), value :: h
         integer(c_int), value :: natoms
         integer(c_int), dimension(*), intent(in) :: types
         real(c_double), dimension(*), intent(in) :: posi_cart
      end subroutine librpa_set_atoms_c

      subroutine librpa_set_kgrids_kvec_c(h, nk1, nk2, nk3, kvecs) bind(c, name="librpa_set_kgrids_kvec")
         import :: c_ptr, c_int, c_double
         type(c_ptr), value :: h
         integer(c_int), value :: nk1, nk2, nk3
         real(c_double), dimension(*), intent(in) :: kvecs
      end subroutine librpa_set_kgrids_kvec_c

      subroutine librpa_set_ibz_mapping_c(h, nkpts, map_ibzk) bind(c, name="librpa_set_ibz_mapping")
         import :: c_ptr, c_int
         type(c_ptr), value :: h
         integer(c_int), value :: nkpts
         integer(c_int), dimension(*), intent(in) :: map_ibzk
      end subroutine librpa_set_ibz_mapping_c

      subroutine librpa_set_lri_coeff_c(h, routing, i_atom, j_atom, nao_i, nao_j, naux_i, &
                                        r, coeff) &
            bind(c, name="librpa_set_lri_coeff")
         import :: c_ptr, c_int, c_double
         type(c_ptr), value :: h
         integer(c_int), value :: routing, i_atom, j_atom, nao_i, nao_j, naux_i
         integer(c_int), dimension(3), intent(in) :: r
         real(c_double), dimension(*), intent(in) :: coeff
      end subroutine librpa_set_lri_coeff_c

      subroutine librpa_set_aux_bare_coulomb_k_atom_pair_c &
            (h, ik, i_atom, j_atom, naux_i, naux_j, vq_real, vq_imag, vq_threshold) &
            bind(c, name="librpa_set_aux_bare_coulomb_k_atom_pair")
         import :: c_ptr, c_int, c_double
         type(c_ptr), value :: h
         integer(c_int), value :: ik, i_atom, j_atom, naux_i, naux_j
         real(c_double), dimension(*), intent(in) :: vq_real, vq_imag
         real(c_double), value :: vq_threshold
      end subroutine librpa_set_aux_bare_coulomb_k_atom_pair_c

      subroutine librpa_set_aux_cut_coulomb_k_atom_pair_c &
            (h, ik, i_atom, j_atom, naux_i, naux_j, vq_real, vq_imag, vq_threshold) &
            bind(c, name="librpa_set_aux_cut_coulomb_k_atom_pair")
         import :: c_ptr, c_int, c_double
         type(c_ptr), value :: h
         integer(c_int), value :: ik, i_atom, j_atom, naux_i, naux_j
         real(c_double), dimension(*), intent(in) :: vq_real, vq_imag
         real(c_double), value :: vq_threshold
      end subroutine librpa_set_aux_cut_coulomb_k_atom_pair_c

      subroutine librpa_set_aux_bare_coulomb_k_2d_block_c &
            (h, ik, mu_begin, mu_end, nu_begin, nu_end, vq_real, vq_imag) &
            bind(c, name="librpa_set_aux_bare_coulomb_k_2d_block")
         import :: c_ptr, c_int, c_double
         type(c_ptr), value :: h
         integer(c_int), value :: ik, mu_begin, mu_end, nu_begin, nu_end
         real(c_double), dimension(*), intent(in) :: vq_real, vq_imag
      end subroutine librpa_set_aux_bare_coulomb_k_2d_block_c

      subroutine librpa_set_aux_cut_coulomb_k_2d_block_c &
            (h, ik, mu_begin, mu_end, nu_begin, nu_end, vq_real, vq_imag) &
            bind(c, name="librpa_set_aux_cut_coulomb_k_2d_block")
         import :: c_ptr, c_int, c_double
         type(c_ptr), value :: h
         integer(c_int), value :: ik, mu_begin, mu_end, nu_begin, nu_end
         real(c_double), dimension(*), intent(in) :: vq_real, vq_imag
      end subroutine librpa_set_aux_cut_coulomb_k_2d_block_c

      subroutine librpa_set_dielect_func_imagfreq_c(h, nfreq, omegas_imag, dielect_func) &
            bind(c, name="librpa_set_dielect_func_imagfreq")
         import :: c_ptr, c_int, c_double
         type(c_ptr), value :: h
         integer(c_int), value :: nfreq
         real(c_double), dimension(*), intent(in) :: omegas_imag
         real(c_double), dimension(*), intent(in) :: dielect_func
      end subroutine librpa_set_dielect_func_imagfreq_c

      subroutine librpa_set_band_kvec_c(h, nkpts_band, kfrac_band) &
            bind(c, name="librpa_set_band_kvec")
         import :: c_ptr, c_int, c_double
         type(c_ptr), value :: h
         integer(c_int), value :: nkpts_band
         real(c_double), dimension(*), intent(in) :: kfrac_band
      end subroutine librpa_set_band_kvec_c

      subroutine librpa_set_band_occ_eigval_c(h, nspins, nkpts_band, nstates, occ, eig) &
            bind(c, name="librpa_set_band_occ_eigval")
         import :: c_ptr, c_int, c_double
         type(c_ptr), value :: h
         integer(c_int), value :: nspins, nkpts_band, nstates
         real(c_double), dimension(*), intent(in) :: occ
         real(c_double), dimension(*), intent(in) :: eig
      end subroutine librpa_set_band_occ_eigval_c

      subroutine librpa_set_wfc_band_c(h, ispin, ik_band, nstates_local, nbasis_local, wfc_real, wfc_imag) &
            bind(c, name="librpa_set_wfc_band")
         import :: c_ptr, c_int, c_double
         type(c_ptr), value :: h
         integer(c_int), value :: ispin, ik_band, nstates_local, nbasis_local
         real(c_double), dimension(*), intent(in) :: wfc_real
         real(c_double), dimension(*), intent(in) :: wfc_imag
      end subroutine librpa_set_wfc_band_c

      subroutine librpa_set_wfc_band_packed_c(h, ispin, ik_band, nstates_local, nbasis_local, wfc) &
            bind(c, name="librpa_set_wfc_band_packed")
         import :: c_ptr, c_int
         type(c_ptr), value :: h
         integer(c_int), value :: ispin, ik_band, nstates_local, nbasis_local
         type(c_ptr), value :: wfc
      end subroutine librpa_set_wfc_band_packed_c

      subroutine librpa_reset_band_data_c(h) bind(c, name="librpa_reset_band_data")
         import :: c_ptr
         type(c_ptr), value :: h
      end subroutine librpa_reset_band_data_c
   end interface

   ! Compute functions interface
   interface
      subroutine librpa_get_imaginary_frequency_grids_c(h, opts, omegas, weights) &
            bind(c, name="librpa_get_imaginary_frequency_grids")
         import :: LibrpaOptions_c, c_ptr, c_double
         type(c_ptr), value :: h
         type(LibrpaOptions_c), intent(in) :: opts
         real(c_double), dimension(*), intent(inout) :: omegas, weights
      end subroutine librpa_get_imaginary_frequency_grids_c

      function librpa_get_rpa_correlation_energy_c(h, opts, nkpts_ibz, contrib_ibzk_re, contrib_ibzk_im) &
            bind(c, name="librpa_get_rpa_correlation_energy")
         import :: LibrpaOptions_c, c_ptr, c_int, c_double
         type(c_ptr), value :: h
         type(LibrpaOptions_c), intent(in) :: opts
         integer(c_int), intent(in), value :: nkpts_ibz
         real(c_double), dimension(*), intent(inout) :: contrib_ibzk_re, contrib_ibzk_im
         real(c_double) :: librpa_get_rpa_correlation_energy_c
      end function librpa_get_rpa_correlation_energy_c

      subroutine librpa_build_exx_c(h, opts) bind(c, name="librpa_build_exx")
         import :: LibrpaOptions_c, c_ptr
         type(c_ptr), value :: h
         type(LibrpaOptions_c), intent(in) :: opts
      end subroutine librpa_build_exx_c

      subroutine librpa_get_exx_pot_kgrid_c(h, opts, n_spins, n_kpts_this, iks_this, &
                                            i_state_low, i_state_high, vexx) &
            bind(c, name="librpa_get_exx_pot_kgrid")
         import :: LibrpaOptions_c, c_ptr, c_int, c_double
         type(c_ptr), value :: h
         type(LibrpaOptions_c), intent(in) :: opts
         integer(c_int), value :: n_spins, n_kpts_this, i_state_low, i_state_high
         integer(c_int), dimension(*), intent(in) :: iks_this
         real(c_double), dimension(*), intent(inout) :: vexx
      end subroutine librpa_get_exx_pot_kgrid_c

      subroutine librpa_get_exx_pot_band_k_c(h, opts, n_spins, n_kpts_band_this, iks_band_this, &
                                             i_state_low, i_state_high, vexx_band) &
            bind(c, name="librpa_get_exx_pot_band_k")
         import :: LibrpaOptions_c, c_ptr, c_int, c_double
         type(c_ptr), value :: h
         type(LibrpaOptions_c), intent(in) :: opts
         integer(c_int), value :: n_spins, n_kpts_band_this, i_state_low, i_state_high
         integer(c_int), dimension(*), intent(in) :: iks_band_this
         real(c_double), dimension(*), intent(inout) :: vexx_band
      end subroutine librpa_get_exx_pot_band_k_c

      subroutine librpa_build_g0w0_sigma_c(h, opts) bind(c, name="librpa_build_g0w0_sigma")
         import :: LibrpaOptions_c, c_ptr
         type(c_ptr), value :: h
         type(LibrpaOptions_c), intent(in) :: opts
      end subroutine librpa_build_g0w0_sigma_c

      subroutine librpa_get_g0w0_qpe_kgrid_c(h, opts, n_spins, n_kpoints_local, iks_local, &
                                             i_state_low, i_state_high, vxc, vexx, sigc_re, sigc_im) &
            bind(c, name="librpa_get_g0w0_qpe_kgrid")
         import :: LibrpaOptions_c, c_ptr, c_int, c_double
         type(c_ptr), value :: h
         type(LibrpaOptions_c), intent(in) :: opts
         integer(c_int), value :: n_spins, n_kpoints_local, i_state_low, i_state_high
         integer(c_int), dimension(*), intent(in) :: iks_local
         real(c_double), dimension(*), intent(in) :: vxc, vexx
         real(c_double), dimension(*), intent(inout) :: sigc_re, sigc_im
      end subroutine librpa_get_g0w0_qpe_kgrid_c
   end interface

   ! Helper to communicate runtime options between C and Fortran types
   interface sync_opt
      module procedure sync_opt_string
      module procedure sync_opt_switch
      module procedure sync_opt_int
      module procedure sync_opt_dp
   end interface

contains

   ! Synchronize C/C++ and Fortran strings
   subroutine sync_opt_string(f_string, c_string, direction)
      implicit none
      character(len=*), intent(inout) :: f_string
      character(len=1, kind=c_char), dimension(*), intent(inout) :: c_string
      integer, intent(in) :: direction

      if (direction .eq. SYNC_OPTS_C2F) then
         call c_f_string_chars(c_string, f_string)
      else if (direction .eq. SYNC_OPTS_F2C) then
         call f_c_string_chars(f_string, c_string, trim_f=.true.)
      end if
   end subroutine sync_opt_string

   ! Synchronize LibrpaSwitch with Fortran logical
   subroutine sync_opt_switch(f_logical, c_switch, direction)
      implicit none
      integer(kind=c_int), intent(inout) :: c_switch
      logical, intent(inout) :: f_logical
      integer, intent(in) :: direction

      if (direction .eq. SYNC_OPTS_C2F) then
         f_logical = (c_switch .eq. LIBRPA_SWITCH_ON)
      else if (direction .eq. SYNC_OPTS_F2C) then
         c_switch = LIBRPA_SWITCH_OFF
         if (f_logical) c_switch = LIBRPA_SWITCH_ON
      end if
   end subroutine sync_opt_switch

   ! Synchronize C/C++ int with Fortran integer
   subroutine sync_opt_int(f_integer, c_int_value, direction)
      implicit none
      integer, intent(inout) :: f_integer
      integer(kind=c_int), intent(inout) :: c_int_value
      integer, intent(in) :: direction

      if (direction .eq. SYNC_OPTS_C2F) then
         f_integer = int(c_int_value)
      else if (direction .eq. SYNC_OPTS_F2C) then
         c_int_value = int(f_integer, kind=c_int)
      end if
   end subroutine sync_opt_int

   ! Synchronize C/C++ and Fortran double precision numbers
   subroutine sync_opt_dp(f_dp, c_double_value, direction)
      implicit none
      real(dp), intent(inout) :: f_dp
      real(c_double), intent(inout) :: c_double_value
      integer, intent(in) :: direction

      if (direction .eq. SYNC_OPTS_C2F) then
         f_dp = real(c_double_value, kind=dp)
      else if (direction .eq. SYNC_OPTS_F2C) then
         c_double_value = real(f_dp, kind=c_double)
      end if
   end subroutine sync_opt_dp

   !> @brief Copy a Fortran character varaible to C char array
   !!
   !> Adapted from https://fortranwiki.org/fortran/show/c_interface_module
   subroutine f_c_string_chars(f_string, c_string, c_string_len, trim_f)
      use iso_c_binding, only: c_null_char
      implicit none

      character(len=*), intent(in) :: f_string
      character(len=1, kind=c_char), dimension(*), intent(out) :: c_string
      ! Max string length, INCLUDING THE TERMINAL NUL
      integer, intent(in), optional :: c_string_len
      logical, intent(in), optional :: trim_f

      integer :: i, strlen

      if (present(trim_f)) then
         if (trim_f) then
            strlen = len(trim(f_string))
         else
            strlen = len(f_string)
         end if
      else
         strlen = len(f_string)
      end if
      ! print*, "strlen ", strlen
      if (present(c_string_len)) then
         if (c_string_len <= 0) return
         strlen = min(strlen, c_string_len - 1)
      end if

      do i = 1, strlen
         c_string(i) = f_string(i:i)
      end do

      c_string(strlen + 1) = c_null_char
   end subroutine f_c_string_chars

   !> @brief Copy a C string, passed as a char-array reference, to a Fortran string.
   !!
   !> copied from https://fortranwiki.org/fortran/show/c_interface_module
   subroutine c_f_string_chars(c_string, f_string)
      use iso_c_binding, only: c_null_char
      implicit none

      character(len=1, kind=c_char), intent(in) :: c_string(*)
      character(len=*), intent(out) :: f_string
      integer :: i
      i=1
      do while (c_string(i) /= c_null_char .and. i <= len(f_string))
         f_string(i:i) = c_string(i)
         i = i + 1
      end do
      if (i < len(f_string)) f_string(i:) = ' '
   end subroutine c_f_string_chars


   !> @brief Convert Fortran logical to C integer as boolean
   subroutine f_c_bool(f_logical, c_bi)
      implicit none
      logical, intent(in) :: f_logical
      integer(kind=c_int), intent(out) :: c_bi

      if (f_logical) then
         c_bi = LIBRPA_SWITCH_ON
      else
         c_bi = LIBRPA_SWITCH_OFF
      endif
   end subroutine f_c_bool

   !> @brief Convert C integer as boolean to Fortran logical
   subroutine c_f_bool(c_bi, f_logical)
      implicit none
      integer(kind=c_int), intent(in) :: c_bi
      logical, intent(out) :: f_logical

      f_logical = (c_bi .eq. LIBRPA_SWITCH_ON)
   end subroutine c_f_bool

   ! Synchronize option values between the Fortran object and the containing C object
   ! Everytime opts_c used through any C interface, its value should be synchronized from opts
   !   call sync_opts(opts, SYNC_OPTS_F2C)
   subroutine sync_opts(opts, direction)
      type(LibrpaOptions), intent(inout) :: opts
      integer, intent(in) :: direction

      if ((direction .ne. SYNC_OPTS_C2F) .and. (direction .ne. SYNC_OPTS_F2C)) then
         stop "internal error - illegal direction value"
      end if

      call sync_opt(opts%output_dir,              opts%opts_c%output_dir,              direction)
      call sync_opt(opts%parallel_routing,        opts%opts_c%parallel_routing,        direction)
      call sync_opt(opts%output_level,            opts%opts_c%output_level,            direction)
      call sync_opt(opts%cs_threshold,            opts%opts_c%cs_threshold,            direction)
      call sync_opt(opts%vq_threshold,            opts%opts_c%vq_threshold,            direction)
      call sync_opt(opts%use_soc,                 opts%opts_c%use_soc,                 direction)
      call sync_opt(opts%use_kpara_scf_eigvec,    opts%opts_c%use_kpara_scf_eigvec,    direction)
      call sync_opt(opts%tfgrids_type,            opts%opts_c%tfgrids_type,            direction)
      call sync_opt(opts%nfreq,                   opts%opts_c%nfreq,                   direction)
      call sync_opt(opts%tfgrids_freq_min,        opts%opts_c%tfgrids_freq_min,        direction)
      call sync_opt(opts%tfgrids_freq_interval,   opts%opts_c%tfgrids_freq_interval,   direction)
      call sync_opt(opts%tfgrids_freq_max,        opts%opts_c%tfgrids_freq_max,        direction)
      call sync_opt(opts%tfgrids_time_min,        opts%opts_c%tfgrids_time_min,        direction)
      call sync_opt(opts%tfgrids_time_interval,   opts%opts_c%tfgrids_time_interval,   direction)
      call sync_opt(opts%gf_threshold,            opts%opts_c%gf_threshold,            direction)
      call sync_opt(opts%use_scalapack_ecrpa,     opts%opts_c%use_scalapack_ecrpa,     direction)
      call sync_opt(opts%n_params_anacon,         opts%opts_c%n_params_anacon,         direction)
      call sync_opt(opts%option_dielect_func,     opts%opts_c%option_dielect_func,     direction)
      call sync_opt(opts%use_scalapack_gw_wc,     opts%opts_c%use_scalapack_gw_wc,     direction)
      call sync_opt(opts%sqrt_coulomb_threshold,  opts%opts_c%sqrt_coulomb_threshold,  direction)
      call sync_opt(opts%replace_w_head,          opts%opts_c%replace_w_head,          direction)
      call sync_opt(opts%libri_chi0_threshold_C,  opts%opts_c%libri_chi0_threshold_C,  direction)
      call sync_opt(opts%libri_chi0_threshold_G,  opts%opts_c%libri_chi0_threshold_G,  direction)
      call sync_opt(opts%libri_exx_threshold_C,   opts%opts_c%libri_exx_threshold_C,   direction)
      call sync_opt(opts%libri_exx_threshold_D,   opts%opts_c%libri_exx_threshold_D,   direction)
      call sync_opt(opts%libri_exx_threshold_V,   opts%opts_c%libri_exx_threshold_V,   direction)
      call sync_opt(opts%libri_g0w0_threshold_C,  opts%opts_c%libri_g0w0_threshold_C,  direction)
      call sync_opt(opts%libri_g0w0_threshold_G,  opts%opts_c%libri_g0w0_threshold_G,  direction)
      call sync_opt(opts%libri_g0w0_threshold_Wc, opts%opts_c%libri_g0w0_threshold_Wc, direction)
      call sync_opt(opts%output_gw_sigc_mat,      opts%opts_c%output_gw_sigc_mat,      direction)
      call sync_opt(opts%output_gw_sigc_mat_rt,   opts%opts_c%output_gw_sigc_mat_rt,   direction)
      call sync_opt(opts%output_gw_sigc_mat_rf,   opts%opts_c%output_gw_sigc_mat_rf,   direction)
   end subroutine

   subroutine librpa_init_options(opts)
      implicit none
      class(LibrpaOptions), intent(inout) :: opts
      call librpa_init_options_c(opts%opts_c)
      call sync_opts(opts, SYNC_OPTS_C2F)
   end subroutine librpa_init_options

   subroutine librpa_set_output_dir(opts, output_dir)
      implicit none
      class(LibrpaOptions), intent(inout) :: opts
      character(len=*), intent(in) :: output_dir
      opts%output_dir = trim(output_dir)
      call sync_opt(opts%output_dir, opts%opts_c%output_dir, SYNC_OPTS_F2C)
   end subroutine librpa_set_output_dir

   !> @brief Initialize the global computing environment of LibRPA
   !>
   !> It should be called after MPI initialization and before other LibRPA functions.
   !>
   !> @param  sw_redirect    Switch of redirecting standard output (default false)
   !> @param  redirect_path  Path of redirected output, only used when `sw_redirect` is true
   !> @param  sw_process     Switch of writing per-process output (default true)
   subroutine librpa_init_global(sw_redirect, redirect_path, sw_process)
      use iso_c_binding, only: c_null_char
      use mpi, only: MPI_COMM_WORLD
      implicit none

      logical, intent(in), optional :: sw_redirect, sw_process
      character(len=*), intent(in), optional :: redirect_path

      character(len=*), parameter :: def = "stdout"
      integer(c_int) :: s1, s2, f_comm
      character(kind=c_char), allocatable, target :: path_c(:)
      character(len=:), allocatable :: tmp
      integer :: n, i

      s1 = LIBRPA_SWITCH_OFF
      if (present(sw_redirect)) then
         if (sw_redirect) s1 = LIBRPA_SWITCH_ON
      end if

      s2 = LIBRPA_SWITCH_ON
      if (present(sw_process)) then
         if (.not. sw_process) s2 = LIBRPA_SWITCH_OFF
      end if

      if (present(redirect_path)) then
        tmp = trim(redirect_path)
      else
        tmp = trim(def)
      end if

      n = len(tmp)
      if (allocated(redirect_path_buf)) deallocate(redirect_path_buf)
      allocate(redirect_path_buf(n+1))
      do i = 1, n
         redirect_path_buf(i) = tmp(i:i)
      end do
      redirect_path_buf(n+1) = c_null_char

      f_comm = int(MPI_COMM_WORLD, kind=c_int)

      call librpa_init_global_c(f_comm, s1, redirect_path_buf, s2)

      !call librpa_init_global_c(s1, redirect_path_buf, s2)
      !if (allocated(path_c)) deallocate(path_c)
   end subroutine librpa_init_global

   !> @brief Release all internal data and finalize the global computing environment of LibRPA
   !>
   !> It should be called after all LibRPA operations are finished.
   subroutine librpa_finalize_global()
      implicit none
      call librpa_finalize_global_c()
      if (allocated(redirect_path_buf)) deallocate(redirect_path_buf)
   end subroutine librpa_finalize_global

   integer function librpa_get_major_version() result(v)
      implicit none
      v = librpa_get_major_version_c()
   end function librpa_get_major_version

   integer function librpa_get_minor_version() result(v)
      implicit none
      v = librpa_get_minor_version_c()
   end function librpa_get_minor_version

   integer function librpa_get_patch_version() result(v)
      implicit none
      v = librpa_get_patch_version_c()
   end function librpa_get_patch_version

   subroutine librpa_test()
      implicit none
      call librpa_test_c()
   end subroutine librpa_test

   subroutine librpa_print_profile()
      implicit none
      call librpa_print_profile_c()
   end subroutine librpa_print_profile

   subroutine librpa_create_handler(this, comm)
      use iso_c_binding, only: c_associated
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      integer, intent(in) :: comm
      integer(c_int) :: f_comm

      if (c_associated(this%ptr_c_handle)) call this%destroy()

      f_comm = int(comm, kind=c_int)
      this%ptr_c_handle = librpa_create_handler_c(f_comm)
      ! this%ptr_c_handle = librpa_create_handler_c(mpi_comm_f2c(comm))
      ! this%ptr_c_handle = librpa_create_handler_c(comm)
   end subroutine librpa_create_handler

   subroutine librpa_destroy_handler(this)
      use iso_c_binding, only: c_associated
      implicit none
      class(LibrpaHandler), intent(inout) :: this

      if (c_associated(this%ptr_c_handle)) then
         call librpa_destroy_handler_c(this%ptr_c_handle)
         this%ptr_c_handle = c_null_ptr
      end if
   end subroutine librpa_destroy_handler

   ! Input functions
   subroutine librpa_set_scf_dimension(this, nspins, nkpts, nstates, nbasis)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      integer, intent(in) :: nspins, nkpts, nstates, nbasis

      integer(c_int) :: nspins_c, nkpts_c, nstates_c, nbasis_c
      ! integer(c_int) :: st_istate_c, nstates_local_c, st_ibasis_c, nbasis_local_c

      nspins_c = int(nspins, kind=c_int)
      nkpts_c = int(nkpts, kind=c_int)
      nstates_c = int(nstates, kind=c_int)
      nbasis_c = int(nbasis, kind=c_int)
      ! st_istate_c = int(st_istate, kind=c_int) - 1
      ! nstates_local_c = int(nstates_local, kind=c_int)
      ! st_ibasis_c = int(st_ibasis, kind=c_int) - 1
      ! nbasis_local_c = int(nbasis_local, kind=c_int)

      call librpa_set_scf_dimension_c(this%ptr_c_handle, nspins_c, nkpts_c, nstates_c, nbasis_c)
   end subroutine librpa_set_scf_dimension

   subroutine librpa_set_wg_ekb_efermi(this, nspins, nkpts, nstates, wg, ekb, efermi)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      integer, intent(in) :: nspins, nkpts, nstates
      real(dp), intent(in) :: wg(nstates, nkpts, nspins)
      real(dp), intent(in) :: ekb(nstates, nkpts, nspins)
      real(dp), intent(in) :: efermi

      integer(c_int) :: nspins_c, nkpts_c, nstates_c
      real(c_double), allocatable :: wg_c(:,:,:), ekb_c(:,:,:)

      nspins_c = int(nspins, kind=c_int)
      nkpts_c = int(nkpts, kind=c_int)
      nstates_c = int(nstates, kind=c_int)

      if (dp == c_double) then
         call librpa_set_wg_ekb_efermi_c(this%ptr_c_handle, nspins_c, nkpts_c, nstates_c, wg, ekb, real(efermi, kind=c_double))
      else
         allocate(wg_c(nstates, nkpts, nspins))
         allocate(ekb_c(nstates, nkpts, nspins))
         wg_c = real(wg, kind=c_double)
         ekb_c = real(ekb, kind=c_double)
         call librpa_set_wg_ekb_efermi_c(this%ptr_c_handle, nspins_c, nkpts_c, nstates_c, wg_c, ekb_c, real(efermi, kind=c_double))
         deallocate(wg_c, ekb_c)
      end if
   end subroutine librpa_set_wg_ekb_efermi

   !> @brief Set the wave-function expansion coefficients
   !> @param ispin          spin index (starting from 1) of the wave function
   !> @param ik             (global) k-point index (starting from 1) of the wave function
   !> @param nstates_local  local dimenstion (number of states) of the parsed wave-function
   !> @param nbasis_local   local dimenstion (number of basis functions) of the parsed wave-function
   !> @param wfc_cplx       Complex-valued wave function to parse
   subroutine librpa_set_wfc(this, ispin, ik, nstates_local, nbasis_local, wfc_cplx)
      use iso_c_binding, only: c_int, c_double, c_loc
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      integer, intent(in) :: ispin, ik, nstates_local, nbasis_local
      complex(dp), intent(in), target :: wfc_cplx(nbasis_local, nstates_local)

      real(c_double), allocatable :: wfc_real(:,:), wfc_imag(:,:)
      integer(c_int) :: ispin_c, ik_c, nstates_local_c, nbasis_local_c

      ispin_c = int(ispin-1, kind=c_int)
      ik_c = int(ik-1, kind=c_int)
      nstates_local_c = int(nstates_local, kind=c_int)
      nbasis_local_c = int(nbasis_local, kind=c_int)

      if (dp == c_double) then
         ! Fast path without create intermediate Fortran arrays
         call librpa_set_wfc_packed_c(&
            this%ptr_c_handle, ispin_c, ik_c, &
            nstates_local_c, nbasis_local_c, c_loc(wfc_cplx))
      else
         allocate(wfc_real(nbasis_local, nstates_local))
         allocate(wfc_imag(nbasis_local, nstates_local))
         wfc_real = real(wfc_cplx, kind=c_double)
         wfc_imag = real(aimag(wfc_cplx), kind=c_double)
         call librpa_set_wfc_c(this%ptr_c_handle, ispin_c, ik_c, &
            nstates_local_c, nbasis_local_c, wfc_real, wfc_imag)
         deallocate(wfc_real, wfc_imag)
      end if
   end subroutine librpa_set_wfc

   subroutine set_ao_basis(h, natoms, nbs, is_aux)
      implicit none
      type(LibrpaHandler), intent(inout) :: h
      integer, intent(in) :: natoms
      integer, intent(in) :: nbs(natoms)
      logical, intent(in) :: is_aux

      integer(c_size_t), allocatable :: nbs_c(:)
      integer(c_int) :: natoms_c

      allocate(nbs_c(natoms))
      nbs_c = int(nbs, kind=c_size_t)
      natoms_c = int(natoms, kind=c_int)
      if (is_aux) then
         call librpa_set_ao_basis_aux_c(h%ptr_c_handle, natoms_c, nbs_c)
      else
         call librpa_set_ao_basis_wfc_c(h%ptr_c_handle, natoms_c, nbs_c)
      end if
      deallocate(nbs_c)
   end subroutine set_ao_basis

   !> @brief Set the wave-function atomic basis
   !>
   !> @param natoms   number of atoms
   !> @param nbs_wfc  number of wave-function basis on each atom
   subroutine librpa_set_ao_basis_wfc(this, natoms, nbs_wfc)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      integer, intent(in) :: natoms
      integer, intent(in) :: nbs_wfc(natoms)

      call set_ao_basis(this, natoms, nbs_wfc, .false.)
   end subroutine librpa_set_ao_basis_wfc

   !> @brief Set the auxiliary atomic basis
   !>
   !> @param natoms   number of atoms
   !> @param nbs_aux  number of auxiliary basis functions on each atom
   subroutine librpa_set_ao_basis_aux(this, natoms, nbs_aux)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      integer, intent(in) :: natoms
      integer, intent(in) :: nbs_aux(natoms)

      call set_ao_basis(this, natoms, nbs_aux, .true.)
   end subroutine librpa_set_ao_basis_aux

   !> @brief Set the direct and reciprocal lattice vectors
   !>
   !> Each column is a lattice/reciprocal lattice vector.
   !>
   !> @param latt     lattice vectors (in Bohr)
   !> @param recplatt reciprocal lattice vectors (in Bohr^-1)
   subroutine librpa_set_latvec_and_G(this, latt, recplatt)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      real(dp), dimension(3, 3), intent(in) :: latt, recplatt

      real(c_double) :: latt_c(3,3), recplatt_c(3,3)

      latt_c = real(latt, kind=c_double)
      recplatt_c = real(recplatt, kind=c_double)
      call librpa_set_latvec_and_G_c(this%ptr_c_handle, latt_c, recplatt_c)
   end subroutine librpa_set_latvec_and_G

   !> @brief Set types and coordinates of the atoms in the model
   !>
   !> @param natoms     number of atoms
   !> @param types      species type of each atom
   !> @param pos_cart   Cartesian coordinates of each atom
   subroutine librpa_set_atoms(this, natoms, types, posi_cart)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      integer, intent(in) :: natoms
      integer, dimension(natoms), intent(in) :: types
      real(dp), dimension(3, natoms), intent(in) :: posi_cart

      integer(c_int) :: natoms_c
      integer(c_int), allocatable :: types_c(:)
      real(c_double), allocatable :: posi_cart_c(:,:)

      allocate(types_c(natoms))

      types_c = int(types, kind=c_int)
      natoms_c = int(natoms, kind=c_int)

      if (dp == c_double) then
         call librpa_set_atoms_c(this%ptr_c_handle, natoms_c, types_c, posi_cart)
      else
         allocate(posi_cart_c(3, natoms))
         posi_cart_c = real(posi_cart, kind=c_double)
         call librpa_set_atoms_c(this%ptr_c_handle, natoms_c, types_c, posi_cart_c)
         deallocate(posi_cart_c)
      end if
      deallocate(types_c)
   end subroutine librpa_set_atoms

   subroutine librpa_set_kgrids_kvec(this, nk1, nk2, nk3, kvecs)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      integer, intent(in) :: nk1, nk2, nk3
      real(dp), intent(in) :: kvecs(3, nk1*nk2*nk3)

      integer(c_int) :: nk1_c, nk2_c, nk3_c
      real(c_double), allocatable :: kvecs_c(:,:)

      nk1_c = int(nk1, kind=c_int)
      nk2_c = int(nk2, kind=c_int)
      nk3_c = int(nk3, kind=c_int)

      if (dp == c_double) then
         call librpa_set_kgrids_kvec_c(this%ptr_c_handle, nk1_c, nk2_c, nk3_c, kvecs)
      else
         allocate(kvecs_c(3, nk1*nk2*nk3))
         kvecs_c = real(kvecs, kind=c_double)
         call librpa_set_kgrids_kvec_c(this%ptr_c_handle, nk1_c, nk2_c, nk3_c, kvecs_c)
         deallocate(kvecs_c)
      end if
   end subroutine librpa_set_kgrids_kvec

   !> @brief Set the mapping from full k-point list to the irreducbile sector
   !>
   !> Example: four-k-point case where the first two and last points are in the irreducbile sector,
   !>          and the third point is mapped to the second, then map_ibzk should be (1, 2, 2, 4)
   !>
   !> @param nkpts     number of k-points in the full Brillouin zone
   !> @param map_ibzk  mapping to the k-point in the irreducible sector
   subroutine librpa_set_ibz_mapping(this, nkpts, map_ibzk)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      integer, intent(in) :: nkpts
      integer, dimension(nkpts), intent(in) :: map_ibzk

      integer :: ik
      integer(c_int) :: nkpts_c
      integer(c_int), allocatable :: map_ibzk_c(:)

      allocate(map_ibzk_c(nkpts))
      map_ibzk_c = int(map_ibzk, kind=c_int) - 1
      nkpts_c = int(nkpts, kind=c_int)
      call librpa_set_ibz_mapping_c(this%ptr_c_handle, nkpts_c, map_ibzk_c)
      deallocate(map_ibzk_c)
   end subroutine librpa_set_ibz_mapping

   !> @brief Set the local RI coefficients
   !>
   !> @param routing  Parallel routing, should be one of the `LIBRPA_ROUTING_*` parameters
   !> @param i_atom   Index of atom I (starting from 1)
   !> @param j_atom   Index of atom J (starting from 1)
   !> @param nao_i    Number of wave-functions basis on atom I
   !> @param nao_j    Number of wave-functions basis on atom J
   !> @param naux_i   Number of auxiliary basis on atom I
   !> @param r        Index of unit cell in the crystal, with (0,0,0) at the origin
   !> @param coeff    Local RI coefficients associated with atom pair I-J, with auxiliary basis on I.
   !>
   subroutine librpa_set_lri_coeff(this, routing, i_atom, j_atom, nao_i, nao_j, naux_i, r, coeff)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      integer, intent(in) :: routing, i_atom, j_atom, nao_i, nao_j, naux_i
      integer, dimension(3), intent(in) :: r
      real(dp), contiguous, intent(in) :: coeff(:, :, :)

      integer(c_int) :: r_c(3)
      integer(c_int) :: routing_c, i_atom_c, j_atom_c, nao_i_c, nao_j_c, naux_i_c
      real(c_double), allocatable :: coeff_c(:,:,:)

      ! Sanity check
      if (size(coeff,1) /= naux_i .or. size(coeff,2) /= nao_j .or. size(coeff,3) /= nao_i) then
         write(*,*) "wrong coeff shape: input (", naux_i, nao_i , nao_j, ") | ", &
                    "internal (", size(coeff,1), size(coeff,2), size(coeff,3), ")"
         error stop "librpa_set_lri_coeff: coeff has wrong shape"
      end if

      ! Check if the routing is a valid LIBRPA_ROUTING_ parameter
      select case (routing)
         case (LIBRPA_ROUTING_AUTO)
         case (LIBRPA_ROUTING_RTAU)
         case (LIBRPA_ROUTING_ATOMPAIR)
         case (LIBRPA_ROUTING_LIBRI)
         case default
            write(*,*) "Invalid routing parameter:", routing
            error stop
      end select

      r_c = int(r, kind=c_int)
      routing_c = int(routing, c_int)
      i_atom_c = int(i_atom-1, c_int)
      j_atom_c = int(j_atom-1, c_int)
      nao_i_c = int(nao_i, c_int)
      nao_j_c = int(nao_j, c_int)
      naux_i_c = int(naux_i, c_int)
      if (dp == c_double) then
         call librpa_set_lri_coeff_c(this%ptr_c_handle, &
               routing_c, i_atom_c, j_atom_c, nao_i_c, nao_j_c, naux_i_c, r_c, coeff)
      else
         allocate(coeff_c(naux_i, nao_j, nao_i))
         coeff_c = real(coeff, kind=c_double)
         call librpa_set_lri_coeff_c(this%ptr_c_handle, &
               routing_c, i_atom_c, j_atom_c, nao_i_c, nao_j_c, naux_i_c, r_c, coeff_c)
         deallocate(coeff_c)
      end if
   end subroutine librpa_set_lri_coeff

   subroutine set_aux_coulomb_k_atom_pair(h, ik, i_atom, j_atom, naux_i, naux_j, vq, vq_threshold, is_cut)
      type(LibrpaHandler), intent(inout) :: h
      integer, intent(in) :: ik, i_atom, j_atom, naux_i, naux_j
      complex(dp), intent(in) :: vq(naux_i, naux_j)
      real(dp), intent(in) :: vq_threshold
      logical, intent(in) :: is_cut

      integer(c_int) :: ik_c, i_atom_c, j_atom_c, naux_i_c, naux_j_c
      real(c_double), allocatable :: vq_real(:,:), vq_imag(:,:)
      real(c_double) :: thres_c

      allocate(vq_real(naux_j, naux_i), vq_imag(naux_j, naux_i))

      ik_c = int(ik-1, kind=c_int)
      i_atom_c = int(i_atom-1, kind=c_int)
      j_atom_c = int(j_atom-1, kind=c_int)
      naux_i_c = int(naux_i, kind=c_int)
      naux_j_c = int(naux_j, kind=c_int)
      vq_real = transpose(real(vq, kind=c_double))
      vq_imag = transpose(real(aimag(vq), kind=c_double))
      thres_c = real(vq_threshold, kind=c_double)
      if (is_cut) then
         call librpa_set_aux_cut_coulomb_k_atom_pair_c(h%ptr_c_handle, &
            ik_c, i_atom_c, j_atom_c, naux_i_c, naux_j_c, vq_real, vq_imag, thres_c)
      else
         call librpa_set_aux_bare_coulomb_k_atom_pair_c(h%ptr_c_handle, &
            ik_c, i_atom_c, j_atom_c, naux_i_c, naux_j_c, vq_real, vq_imag, thres_c)
      end if

      deallocate(vq_real, vq_imag)
   end subroutine set_aux_coulomb_k_atom_pair 

   subroutine librpa_set_aux_bare_coulomb_k_atom_pair &
         (this, ik, i_atom, j_atom, naux_i, naux_j, vq, vq_threshold)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      integer, intent(in) :: ik, i_atom, j_atom, naux_i, naux_j
      complex(dp), intent(in) :: vq(naux_i, naux_j)
      real(dp), intent(in) :: vq_threshold

      call set_aux_coulomb_k_atom_pair(this, ik, i_atom, j_atom, naux_i, naux_j, vq, vq_threshold, .false.)
   end subroutine librpa_set_aux_bare_coulomb_k_atom_pair

   subroutine librpa_set_aux_cut_coulomb_k_atom_pair &
         (this, ik, i_atom, j_atom, naux_i, naux_j, vq, vq_threshold)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      integer, intent(in) :: ik, i_atom, j_atom, naux_i, naux_j
      complex(dp), intent(in) :: vq(naux_i, naux_j)
      real(dp), intent(in) :: vq_threshold

      call set_aux_coulomb_k_atom_pair(this, ik, i_atom, j_atom, naux_i, naux_j, vq, vq_threshold, .true.)
   end subroutine librpa_set_aux_cut_coulomb_k_atom_pair

   subroutine set_aux_coulomb_k_2d_block(h, ik, mu_begin, mu_end, nu_begin, nu_end, vq, is_cut)
      implicit none
      type(LibrpaHandler), intent(inout) :: h
      integer, intent(in) :: ik, mu_begin, mu_end, nu_begin, nu_end
      complex(dp), intent(in) :: vq(mu_end-mu_begin+1, nu_end-nu_begin+1)
      logical, intent(in) :: is_cut

      integer(c_int) :: ik_c, mb, me, nb, ne
      real(c_double), allocatable :: vq_real(:,:), vq_imag(:,:)

      ik_c = int(ik-1, kind=c_int)
      ! The C interface uses included beginning and excluded ending.
      ! In Fortran, it is more common that both of two ends are included
      mb = int(mu_begin-1, kind=c_int)
      me = int(mu_end, kind=c_int)
      nb = int(nu_begin-1, kind=c_int)
      ne = int(nu_end, kind=c_int)

      allocate(vq_real(ne-nb, me-mb))
      allocate(vq_imag(ne-nb, me-mb))
      vq_real = transpose(real(vq, kind=c_double))
      vq_imag = transpose(real(aimag(vq), kind=c_double))
      if (is_cut) then
         call librpa_set_aux_cut_coulomb_k_2d_block_c(h%ptr_c_handle, ik_c, mb, me, nb, ne, vq_real, vq_imag)
      else
         call librpa_set_aux_bare_coulomb_k_2d_block_c(h%ptr_c_handle, ik_c, mb, me, nb, ne, vq_real, vq_imag)
      end if
      deallocate(vq_real)
      deallocate(vq_imag)
   end subroutine set_aux_coulomb_k_2d_block

   subroutine librpa_set_aux_bare_coulomb_k_2d_block &
         (this, ik, mu_begin, mu_end, nu_begin, nu_end, vq)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      integer, intent(in) :: ik, mu_begin, mu_end, nu_begin, nu_end
      complex(dp), intent(in) :: vq(mu_end-mu_begin+1, nu_end-nu_begin+1)

      call set_aux_coulomb_k_2d_block(this, ik, mu_begin, mu_end, nu_begin, nu_end, vq, .false.)
   end subroutine librpa_set_aux_bare_coulomb_k_2d_block

   subroutine librpa_set_aux_cut_coulomb_k_2d_block &
         (this, ik, mu_begin, mu_end, nu_begin, nu_end, vq)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      integer, intent(in) :: ik, mu_begin, mu_end, nu_begin, nu_end
      complex(dp), intent(in) :: vq(mu_end-mu_begin+1, nu_end-nu_begin+1)

      call set_aux_coulomb_k_2d_block(this, ik, mu_begin, mu_end, nu_begin, nu_end, vq, .true.)
   end subroutine librpa_set_aux_cut_coulomb_k_2d_block

   subroutine librpa_set_dielect_func_imagfreq(this, nfreq, omegas_imag, dielect_func)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      integer, intent(in) :: nfreq
      real(dp), dimension(nfreq), intent(in) :: omegas_imag
      real(dp), dimension(nfreq), intent(in) :: dielect_func

      integer(c_int) :: nfreq_c
      real(c_double), allocatable :: omegas_c(:), df_c(:)

      nfreq_c = int(nfreq, kind=c_int)
      if (dp == c_double) then
         call librpa_set_dielect_func_imagfreq_c(this%ptr_c_handle, nfreq_c, omegas_imag, dielect_func)
      else
         allocate(omegas_c(nfreq))
         allocate(df_c(nfreq))
         omegas_c(:) = real(omegas_imag(:), kind=c_double)
         df_c(:) = real(dielect_func(:), kind=c_double)
         call librpa_set_dielect_func_imagfreq_c(this%ptr_c_handle, nfreq_c, omegas_c, df_c)
         deallocate(omegas_c, df_c)
      end if
   end subroutine librpa_set_dielect_func_imagfreq

   subroutine librpa_set_band_kvec(this, nkpts_band, kfrac_band)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      integer, intent(in) :: nkpts_band
      real(dp), intent(in) :: kfrac_band(3, nkpts_band)

      integer(c_int) :: nkb_c
      real(c_double), allocatable :: kfb_c(:,:)

      nkb_c = int(nkpts_band, kind=c_int)

      if (dp == c_double) then
         call librpa_set_band_kvec_c(this%ptr_c_handle, nkb_c, kfrac_band)
      else
         allocate(kfb_c(3, nkpts_band))
         kfb_c = real(kfrac_band, kind=c_double)
         call librpa_set_band_kvec_c(this%ptr_c_handle, nkb_c, kfb_c)
         deallocate(kfb_c)
      end if
   end subroutine librpa_set_band_kvec

   subroutine librpa_set_band_occ_eigval(this, nspins, nkpts_band, nstates, occ, eig)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      integer, intent(in) :: nspins, nkpts_band, nstates
      real(dp), intent(in) :: occ(nstates, nkpts_band, nspins)
      real(dp), intent(in) :: eig(nstates, nkpts_band, nspins)

      integer(c_int) :: nspins_c, nkpts_band_c, nstates_c
      real(c_double), allocatable :: occ_c(:,:,:), eig_c(:,:,:)

      nspins_c = int(nspins, kind=c_int)
      nkpts_band_c = int(nkpts_band, kind=c_int)
      nstates_c = int(nstates, kind=c_int)

      if (dp == c_double) then
         call librpa_set_band_occ_eigval_c(this%ptr_c_handle, nspins_c, nkpts_band_c, nstates_c, occ, eig)
      else
         allocate(occ_c(nstates, nkpts_band, nspins))
         allocate(eig_c(nstates, nkpts_band, nspins))
         occ_c = real(occ, kind=c_double)
         eig_c = real(eig, kind=c_double)
         call librpa_set_band_occ_eigval_c(this%ptr_c_handle, nspins_c, nkpts_band_c, nstates_c, occ_c, eig_c)
         deallocate(occ_c, eig_c)
      end if
   end subroutine librpa_set_band_occ_eigval

   !> @brief Set the wave-function expansion coefficients for band calculation
   !> @param ispin          spin index (starting from 1) of the wave function
   !> @param ik_band        (global) k-point index (starting from 1) of the wave function
   !> @param nstates_local  local dimenstion (number of states) of the parsed wave-function
   !> @param nbasis_local   local dimenstion (number of basis functions) of the parsed wave-function
   !> @param wfc_cplx       Complex-valued wave function to parse
   subroutine librpa_set_wfc_band(this, ispin, ik_band, nstates_local, nbasis_local, wfc_cplx)
      use iso_c_binding, only: c_int, c_double, c_loc
      implicit none

      class(LibrpaHandler), intent(inout) :: this
      integer, intent(in) :: ispin, ik_band, nstates_local, nbasis_local
      complex(dp), intent(in), target :: wfc_cplx(nbasis_local, nstates_local)

      real(c_double), allocatable :: wfc_real(:,:), wfc_imag(:,:)
      integer(c_int) :: ispin_c, ikb_c, nstates_local_c, nbasis_local_c

      ispin_c = int(ispin-1, kind=c_int)
      ikb_c = int(ik_band-1, kind=c_int)
      nstates_local_c = int(nstates_local, kind=c_int)
      nbasis_local_c = int(nbasis_local, kind=c_int)

      if (dp == c_double) then
         ! Fast path without create intermediate Fortran arrays
         call librpa_set_wfc_band_packed_c(&
            this%ptr_c_handle, ispin_c, ikb_c, &
            nstates_local_c, nbasis_local_c, c_loc(wfc_cplx))
      else
         allocate(wfc_real(nbasis_local, nstates_local))
         allocate(wfc_imag(nbasis_local, nstates_local))
         wfc_real = real(wfc_cplx, kind=c_double)
         wfc_imag = real(aimag(wfc_cplx), kind=c_double)
         call librpa_set_wfc_band_c(this%ptr_c_handle, ispin_c, ikb_c, &
            nstates_local_c, nbasis_local_c, wfc_real, wfc_imag)
         deallocate(wfc_real, wfc_imag)
      end if
   end subroutine librpa_set_wfc_band

   subroutine librpa_reset_band_data(this)
      implicit none
      class(LibrpaHandler), intent(inout) :: this

      call librpa_reset_band_data_c(this%ptr_c_handle)
   end subroutine librpa_reset_band_data

   ! Compute functions
   subroutine librpa_get_imaginary_frequency_grids(this, opts, omegas, weights)
      class(LibrpaHandler), intent(inout) :: this
      type(LibrpaOptions), intent(inout) :: opts
      real(dp), allocatable, intent(inout) :: omegas(:), weights(:)

      real(c_double), allocatable :: omegas_c(:), weights_c(:)

      if (allocated(omegas)) then
         if (size(omegas) .ne. opts%nfreq) deallocate(omegas)
      end if
      if (.not. allocated(omegas)) allocate(omegas(opts%nfreq))

      if (allocated(weights)) then
         if (size(weights) .ne. opts%nfreq) deallocate(weights)
      end if
      if (.not. allocated(weights)) allocate(weights(opts%nfreq))

      ! write(*, *) "size(omegas): ", size(omegas), " size(weights): ", size(weights)

      call sync_opts(opts, SYNC_OPTS_F2C)
      if (c_double == dp) then
         call librpa_get_imaginary_frequency_grids_c(this%ptr_c_handle, opts%opts_c, omegas, weights)
      else
         allocate(omegas_c(opts%nfreq))
         allocate(weights_c(opts%nfreq))
         call librpa_get_imaginary_frequency_grids_c(this%ptr_c_handle, opts%opts_c, omegas_c, weights_c)
         omegas(:) = real(omegas_c(:), kind=dp)
         weights(:) = real(weights_c(:), kind=dp)
         deallocate(omegas_c, weights_c)
      end if
   end subroutine librpa_get_imaginary_frequency_grids

   real(dp) function librpa_get_rpa_correlation_energy(this, opts, nkpts_ibz, contrib_ibzk) result(e)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      type(LibrpaOptions), intent(inout) :: opts
      integer, intent(in) :: nkpts_ibz
      complex(dp), dimension(nkpts_ibz), intent(inout) :: contrib_ibzk

      integer(c_int) :: nkpts_ibz_c
      real(c_double), allocatable :: contrib_ibzk_re(:), contrib_ibzk_im(:)

      call sync_opts(opts, SYNC_OPTS_F2C)

      nkpts_ibz_c = int(nkpts_ibz, kind=c_int)
      allocate(contrib_ibzk_re(nkpts_ibz), contrib_ibzk_im(nkpts_ibz))
      e = librpa_get_rpa_correlation_energy_c(this%ptr_c_handle, opts%opts_c, &
                                              nkpts_ibz_c, contrib_ibzk_re, contrib_ibzk_im)
      if (dp == c_double) then
         contrib_ibzk = contrib_ibzk_re + contrib_ibzk_im * CIMAG
      else
         contrib_ibzk = real(contrib_ibzk_re, kind=dp) + real(contrib_ibzk_im, kind=dp) * CIMAG
      end if
      deallocate(contrib_ibzk_re, contrib_ibzk_im)
   end function librpa_get_rpa_correlation_energy

   subroutine librpa_build_exx(this, opts)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      type(LibrpaOptions), intent(inout) :: opts

      call sync_opts(opts, SYNC_OPTS_F2C)
      call librpa_build_exx_c(this%ptr_c_handle, opts%opts_c)
   end subroutine librpa_build_exx

   subroutine librpa_get_exx_pot_kgrid(this, opts, n_spins, n_kpts_this, iks_this, &
                                       i_state_low, i_state_high, vexx)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      type(LibrpaOptions), intent(inout) :: opts
      integer, contiguous, dimension(:), intent(in) :: iks_this
      integer, intent(in) :: n_spins, n_kpts_this, i_state_low, i_state_high
      real(dp), dimension(i_state_high - i_state_low + 1, n_kpts_this, n_spins), intent(inout) :: vexx

      integer(c_int), allocatable :: iks_this_c(:)
      real(c_double), allocatable :: vexx_c(:,:,:)
      integer(c_int) :: n_spins_c, n_kpts_this_c, i_state_low_c, i_state_high_c
      integer :: n_states_calc

      n_spins_c = int(n_spins, kind=c_int)
      i_state_low_c = int(i_state_low - 1, kind=c_int)
      ! i_state_high is not included in the C interface, so no minus 1 here
      i_state_high_c = int(i_state_high, kind=c_int)

      n_kpts_this_c = int(n_kpts_this, kind=c_int)
      allocate(iks_this_c(max(1, n_kpts_this)))
      if (n_kpts_this > 0) then
         iks_this_c = int(iks_this(1:n_kpts_this), kind=c_int) - 1
         ! write(*,*) "size(iks_local) ", size(iks_local)
         ! write(*,*) "size(iks_local_c) ", size(iks_local_c)
         ! write(*,*) "iks_local_c ", iks_local_c
      end if

      call sync_opts(opts, SYNC_OPTS_F2C)
      if (dp == c_double) then
         call librpa_get_exx_pot_kgrid_c(this%ptr_c_handle, opts%opts_c, n_spins_c, n_kpts_this_c, &
                                         iks_this_c, i_state_low_c, i_state_high_c, vexx)
      else
         n_states_calc = i_state_high - i_state_low + 1
         allocate(vexx_c(n_states_calc, max(n_kpts_this, 1), n_spins))
         call librpa_get_exx_pot_kgrid_c(this%ptr_c_handle, opts%opts_c, n_spins_c, n_kpts_this_c, &
                                         iks_this_c, i_state_low_c, i_state_high_c, vexx_c)
         if (n_kpts_this > 0) vexx = real(vexx_c, kind=dp)
         deallocate(vexx_c)
      end if

      deallocate(iks_this_c)
   end subroutine librpa_get_exx_pot_kgrid

   subroutine librpa_get_exx_pot_band_k(this, opts, n_spins, n_kpts_band_this, iks_band_this, &
                                        i_state_low, i_state_high, vexx_band)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      type(LibrpaOptions), intent(inout) :: opts
      integer, contiguous, dimension(:), intent(in) :: iks_band_this
      integer, intent(in) :: n_spins, n_kpts_band_this, i_state_low, i_state_high
      real(dp), dimension(i_state_high - i_state_low + 1, n_kpts_band_this, n_spins), intent(inout) :: vexx_band

      integer(c_int), allocatable :: iks_band_this_c(:)
      real(c_double), allocatable :: vexx_band_c(:,:,:)
      integer(c_int) :: n_spins_c, n_kpts_this_c, i_state_low_c, i_state_high_c
      integer :: n_states_calc

      n_spins_c = int(n_spins, kind=c_int)
      i_state_low_c = int(i_state_low - 1, kind=c_int)
      ! i_state_high is not included in the C interface, so no minus 1 here
      i_state_high_c = int(i_state_high, kind=c_int)

      n_kpts_this_c = int(n_kpts_band_this, kind=c_int)
      allocate(iks_band_this_c(max(1, n_kpts_band_this)))
      if (n_kpts_band_this > 0) then
         iks_band_this_c = int(iks_band_this(1:n_kpts_band_this), kind=c_int) - 1
         ! write(*,*) "size(iks_local) ", size(iks_local)
         ! write(*,*) "size(iks_local_c) ", size(iks_local_c)
         ! write(*,*) "iks_local_c ", iks_local_c
      end if

      call sync_opts(opts, SYNC_OPTS_F2C)
      if (dp == c_double) then
         call librpa_get_exx_pot_band_k_c(this%ptr_c_handle, opts%opts_c, n_spins_c, n_kpts_this_c, &
                                          iks_band_this_c, i_state_low_c, i_state_high_c, vexx_band)
      else
         n_states_calc = i_state_high - i_state_low + 1
         allocate(vexx_band_c(n_states_calc, max(n_kpts_band_this, 1), n_spins))
         call librpa_get_exx_pot_band_k_c(this%ptr_c_handle, opts%opts_c, n_spins_c, n_kpts_this_c, &
                                          iks_band_this_c, i_state_low_c, i_state_high_c, vexx_band_c)
         if (n_kpts_band_this > 0) vexx_band = real(vexx_band_c, kind=dp)
         deallocate(vexx_band_c)
      end if

      deallocate(iks_band_this_c)
   end subroutine librpa_get_exx_pot_band_k

   subroutine librpa_build_g0w0_sigma(this, opts)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      type(LibrpaOptions), intent(inout) :: opts

      call sync_opts(opts, SYNC_OPTS_F2C)
      call librpa_build_g0w0_sigma_c(this%ptr_c_handle, opts%opts_c)
   end subroutine librpa_build_g0w0_sigma

   subroutine librpa_get_g0w0_qpe_kgrid(this, opts, n_spins, n_kpoints_local, iks_local, &
                                        i_state_low, i_state_high, vxc, vexx, sigc)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      type(LibrpaOptions), intent(inout) :: opts
      integer, contiguous, dimension(:), intent(in) :: iks_local
      integer, intent(in) :: n_spins, n_kpoints_local, i_state_low, i_state_high
      real(dp), dimension(i_state_high - i_state_low + 1, n_kpoints_local, n_spins), intent(in) :: vxc
      real(dp), dimension(i_state_high - i_state_low + 1, n_kpoints_local, n_spins), intent(in) :: vexx
      complex(dp), dimension(i_state_high - i_state_low + 1, n_kpoints_local, n_spins), intent(inout) :: sigc

      integer(c_int), allocatable :: iks_local_c(:)
      real(c_double), allocatable :: vxc_c(:,:,:), vexx_c(:,:,:)
      real(c_double), allocatable :: sigc_re_c(:,:,:), sigc_im_c(:,:,:)
      integer(c_int) :: n_spins_c, n_kpoints_local_c, i_state_low_c, i_state_high_c
      integer :: n_states_calc

      n_spins_c = int(n_spins, kind=c_int)
      i_state_low_c = int(i_state_low - 1, kind=c_int)
      i_state_high_c = int(i_state_high, kind=c_int)

      n_kpoints_local_c = int(n_kpoints_local, kind=c_int)
      if (n_kpoints_local > 0) then
         allocate(iks_local_c(n_kpoints_local))
         iks_local_c = int(iks_local(1:n_kpoints_local), kind=c_int) - 1
      else
         allocate(iks_local_c(1))
      end if

      n_states_calc = i_state_high - i_state_low + 1
      allocate(sigc_re_c(n_states_calc, max(n_kpoints_local, 1), n_spins))
      allocate(sigc_im_c(n_states_calc, max(n_kpoints_local, 1), n_spins))

      call sync_opts(opts, SYNC_OPTS_F2C)
      if (dp == c_double) then
         call librpa_get_g0w0_qpe_kgrid_c(this%ptr_c_handle, opts%opts_c, n_spins_c, n_kpoints_local_c, &
                                          iks_local_c, i_state_low_c, i_state_high_c, vxc, vexx, sigc_re_c, sigc_im_c)
      else
         allocate(vexx_c(n_states_calc, max(n_kpoints_local, 1), n_spins))
         allocate(vxc_c(n_states_calc, max(n_kpoints_local, 1), n_spins))
         if (n_kpoints_local > 0) then
            vxc_c = real(vxc, kind=c_double)
            vexx_c = real(vexx, kind=c_double)
         end if
         call librpa_get_g0w0_qpe_kgrid_c(this%ptr_c_handle, opts%opts_c, n_spins_c, n_kpoints_local_c, &
                                          iks_local_c, i_state_low_c, i_state_high_c, vxc_c, vexx_c, sigc_re_c, sigc_im_c)
         if (allocated(vxc_c)) deallocate(vxc_c)
         if (allocated(vexx_c)) deallocate(vexx_c)
      end if

      if (n_kpoints_local > 0) then
         sigc(:,:,:) = cmplx(sigc_re_c, sigc_im_c, kind=dp)
      end if

      deallocate(sigc_re_c)
      deallocate(sigc_im_c)
      deallocate(iks_local_c)
   end subroutine librpa_get_g0w0_qpe_kgrid

end module librpa_f03
