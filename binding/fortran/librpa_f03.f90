module librpa_f03

   use iso_c_binding, only: c_char, c_ptr, c_int, c_double, c_null_ptr, c_size_t
   implicit none

   private

   !=======================================================================
   ! Public types, constants, and functions
   !=======================================================================
   public :: LibrpaOptions
   public :: LibrpaHandler

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

   !=======================================================================
   ! Switch as bool type, defined as in include/librpa_enums.h
   integer(c_int), parameter :: LIBRPA_SWITCH_OFF = 0
   integer(c_int), parameter :: LIBRPA_SWITCH_ON = 1
   character(kind=c_char), allocatable, target, save :: redirect_path_buf(:)

   complex(8), parameter :: CIMAG = (0.0d0, 1.0d0)

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
      real(8) :: cs_threshold
      real(8) :: vq_threshold
      logical :: use_soc
      integer :: tfgrids_type
      integer :: nfreq
      real(8) :: tfgrids_freq_min
      real(8) :: tfgrids_freq_interval
      real(8) :: tfgrids_freq_max
      real(8) :: tfgrids_time_min
      real(8) :: tfgrids_time_interval
      real(8) :: gf_threshold
      logical :: use_scalapack_ecrpa
      integer :: n_params_anacon
      logical :: use_scalapack_gw_wc
      logical :: replace_w_head
      integer :: option_dielect_func
      real(8) :: sqrt_coulomb_threshold
      real(8) :: libri_chi0_threshold_C
      real(8) :: libri_chi0_threshold_G
      real(8) :: libri_exx_threshold_C
      real(8) :: libri_exx_threshold_D
      real(8) :: libri_exx_threshold_V
      real(8) :: libri_g0w0_threshold_C
      real(8) :: libri_g0w0_threshold_G
      real(8) :: libri_g0w0_threshold_Wc
      logical :: output_gw_sigc_mat
      logical :: output_gw_sigc_mat_rt
      logical :: output_gw_sigc_mat_rf

      contains
         procedure :: init => librpa_init_options
   end type LibrpaOptions

   interface
      ! void librpa_init_options(LibrpaOptions *opts);
      subroutine librpa_init_options_c(opts_c) bind(c, name="librpa_init_options")
         import :: LibrpaOptions_c
         type(LibrpaOptions_c) :: opts_c
      end subroutine librpa_init_options_c
   end interface

   integer, parameter :: SYNC_OPTS_C2F = 1
   integer, parameter :: SYNC_OPTS_F2C = -1

   ! Global environment
   interface
      subroutine librpa_init_global_c(sw_redirect, path, sw_process) bind(c, name="librpa_init_global")
         import :: c_int, c_char
         integer(c_int), value :: sw_redirect
         character(kind=c_char), dimension(*), intent(in) :: path
         integer(c_int), value :: sw_process
      end subroutine librpa_init_global_c

      subroutine librpa_finalize_global_c() bind(c, name="librpa_finalize_global")
      end subroutine librpa_finalize_global_c
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
         ! Compute
         procedure :: get_rpa_correlation_energy => librpa_get_rpa_correlation_energy
   end type LibrpaHandler

   interface
      function librpa_create_handler_c(comm_c) bind(c, name="librpa_create_handler")
         import :: c_ptr, c_int
         integer(c_int), value :: comm_c
         type(c_ptr) :: librpa_create_handler_c
      end function librpa_create_handler_c

      subroutine librpa_destroy_handler_c(h) bind(c, name="librpa_destroy_handler")
         import :: c_ptr
         type(c_ptr), value :: h
      end subroutine librpa_destroy_handler_c
   end interface

   ! Input functions interface
   interface
      subroutine librpa_set_scf_dimension_c(h, nspins, nkpts, nstates, nbasis, &
                                            st_istate, nstates_local, st_ibasis, nbasis_local) &
                                            bind(c, name="librpa_set_scf_dimension")
         import :: c_ptr, c_int
         type(c_ptr), value :: h
         integer(c_int), value :: nspins, nkpts, nstates, nbasis, st_istate, nstates_local, st_ibasis, nbasis_local
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

      function librpa_get_rpa_correlation_energy_c(h, opts, nkpts_ibz, contrib_ibzk_re, contrib_ibzk_im) &
            bind(c, name="librpa_get_rpa_correlation_energy")
         import :: LibrpaOptions_c, c_ptr, c_int, c_double
         type(c_ptr), value :: h
         type(LibrpaOptions_c), intent(in) :: opts
         integer(c_int), intent(in) :: nkpts_ibz
         real(c_double), dimension(*), intent(inout) :: contrib_ibzk_re, contrib_ibzk_im
         real(c_double) :: librpa_get_rpa_correlation_energy_c
      end function librpa_get_rpa_correlation_energy_c
   end interface

contains

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

      if (direction .eq. SYNC_OPTS_C2F) then
         ! From C object to Fortran. This case is rare, usually only for initialization
         call c_f_string_chars(opts%opts_c%output_dir, opts%output_dir)
         opts%parallel_routing = opts%opts_c%parallel_routing
         opts%output_level = opts%opts_c%output_level
         opts%cs_threshold = opts%opts_c%cs_threshold
         opts%vq_threshold = opts%opts_c%vq_threshold
         call c_f_bool(opts%opts_c%use_soc, opts%use_soc)
         opts%tfgrids_type = opts%opts_c%tfgrids_type
         opts%nfreq = opts%opts_c%nfreq
         opts%tfgrids_freq_min = opts%opts_c%tfgrids_freq_min
         opts%tfgrids_freq_interval = opts%opts_c%tfgrids_freq_interval
         opts%tfgrids_freq_max = opts%opts_c%tfgrids_freq_max
         opts%tfgrids_time_min = opts%opts_c%tfgrids_time_min
         opts%tfgrids_time_interval = opts%opts_c%tfgrids_time_interval
         opts%gf_threshold = opts%opts_c%gf_threshold
         call c_f_bool(opts%opts_c%use_scalapack_ecrpa, opts%use_scalapack_ecrpa)
         opts%n_params_anacon = opts%opts_c%n_params_anacon
         call c_f_bool(opts%opts_c%use_scalapack_gw_wc, opts%use_scalapack_gw_wc)
         call c_f_bool(opts%opts_c%replace_w_head, opts%replace_w_head)
         opts%option_dielect_func = opts%opts_c%option_dielect_func
         opts%sqrt_coulomb_threshold = opts%opts_c%sqrt_coulomb_threshold
         opts%libri_chi0_threshold_C = opts%opts_c%libri_chi0_threshold_C
         opts%libri_chi0_threshold_G = opts%opts_c%libri_chi0_threshold_G
         opts%libri_exx_threshold_C = opts%opts_c%libri_exx_threshold_C
         opts%libri_exx_threshold_D = opts%opts_c%libri_exx_threshold_D
         opts%libri_exx_threshold_V = opts%opts_c%libri_exx_threshold_V
         opts%libri_g0w0_threshold_C = opts%opts_c%libri_g0w0_threshold_C
         opts%libri_g0w0_threshold_G = opts%opts_c%libri_g0w0_threshold_G
         opts%libri_g0w0_threshold_Wc = opts%opts_c%libri_g0w0_threshold_Wc
         call c_f_bool(opts%opts_c%output_gw_sigc_mat, opts%output_gw_sigc_mat)
         call c_f_bool(opts%opts_c%output_gw_sigc_mat_rt, opts%output_gw_sigc_mat_rt)
         call c_f_bool(opts%opts_c%output_gw_sigc_mat_rf, opts%output_gw_sigc_mat_rf)
      else if (direction .eq. SYNC_OPTS_F2C) then
         ! From Fortran object to C, should be called at the beginning of each compute API function
         call f_c_string_chars(opts%output_dir, opts%opts_c%output_dir, trim_f=.true.)
         opts%opts_c%parallel_routing = opts%parallel_routing
         opts%opts_c%output_level = opts%output_level
         opts%opts_c%cs_threshold = opts%cs_threshold
         opts%opts_c%vq_threshold = opts%vq_threshold
         call f_c_bool(opts%use_soc, opts%opts_c%use_soc)
         opts%opts_c%tfgrids_type = opts%tfgrids_type
         opts%opts_c%nfreq = opts%nfreq
         opts%opts_c%tfgrids_freq_min = opts%tfgrids_freq_min
         opts%opts_c%tfgrids_freq_interval = opts%tfgrids_freq_interval
         opts%opts_c%tfgrids_freq_max = opts%tfgrids_freq_max
         opts%opts_c%tfgrids_time_min = opts%tfgrids_time_min
         opts%opts_c%tfgrids_time_interval = opts%tfgrids_time_interval
         opts%opts_c%gf_threshold = opts%gf_threshold
         call f_c_bool(opts%use_scalapack_ecrpa, opts%opts_c%use_scalapack_ecrpa)
         opts%opts_c%n_params_anacon = opts%n_params_anacon
         call f_c_bool(opts%use_scalapack_gw_wc, opts%opts_c%use_scalapack_gw_wc)
         call f_c_bool(opts%replace_w_head, opts%opts_c%replace_w_head)
         opts%opts_c%option_dielect_func = opts%option_dielect_func
         opts%opts_c%sqrt_coulomb_threshold = opts%sqrt_coulomb_threshold
         opts%opts_c%libri_chi0_threshold_C = opts%libri_chi0_threshold_C
         opts%opts_c%libri_chi0_threshold_G = opts%libri_chi0_threshold_G
         opts%opts_c%libri_exx_threshold_C = opts%libri_exx_threshold_C
         opts%opts_c%libri_exx_threshold_D = opts%libri_exx_threshold_D
         opts%opts_c%libri_exx_threshold_V = opts%libri_exx_threshold_V
         opts%opts_c%libri_g0w0_threshold_C = opts%libri_g0w0_threshold_C
         opts%opts_c%libri_g0w0_threshold_G = opts%libri_g0w0_threshold_G
         opts%opts_c%libri_g0w0_threshold_Wc = opts%libri_g0w0_threshold_Wc
         call f_c_bool(opts%output_gw_sigc_mat, opts%opts_c%output_gw_sigc_mat)
         call f_c_bool(opts%output_gw_sigc_mat_rt, opts%opts_c%output_gw_sigc_mat_rt)
         call f_c_bool(opts%output_gw_sigc_mat_rf, opts%opts_c%output_gw_sigc_mat_rf)
      else
         stop "internal error - illegal direction value"
      end if
   end subroutine

   subroutine librpa_init_options(opts)
      implicit none
      class(LibrpaOptions), intent(inout) :: opts
      call librpa_init_options_c(opts%opts_c)
      call sync_opts(opts, SYNC_OPTS_C2F)
   end subroutine librpa_init_options

   subroutine librpa_init_global(sw_redirect, redirect_path, sw_process)
      use iso_c_binding, only: c_null_char
      implicit none

      logical, intent(in), optional :: sw_redirect, sw_process
      character(len=*), intent(in), optional :: redirect_path

      character(len=*), parameter :: def = "stdout"
      integer(c_int) :: s1, s2
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

      call librpa_init_global_c(s1, redirect_path_buf, s2)
      !if (allocated(path_c)) deallocate(path_c)
   end subroutine librpa_init_global

   subroutine librpa_finalize_global()
      implicit none
      call librpa_finalize_global_c()
      if (allocated(redirect_path_buf)) deallocate(redirect_path_buf)
   end subroutine librpa_finalize_global

   integer function librpa_get_major_version() result(v)
      integer(c_int) :: v_c
      v_c = librpa_get_major_version_c()
      v = int(v_c)
   end function librpa_get_major_version

   integer function librpa_get_minor_version() result(v)
      integer(c_int) :: v_c
      v_c = librpa_get_minor_version_c()
      v = int(v_c)
   end function librpa_get_minor_version

   integer function librpa_get_patch_version() result(v)
      integer(c_int) :: v_c
      v_c = librpa_get_patch_version_c()
      v = int(v_c)
   end function librpa_get_patch_version

   subroutine librpa_create_handler(this, comm)
      use iso_c_binding, only: c_associated
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      integer, intent(in) :: comm
      integer(c_int) :: comm_c = 0

      if (c_associated(this%ptr_c_handle)) call this%destroy()

      comm_c = int(comm, kind=c_int)
      this%ptr_c_handle = librpa_create_handler_c(comm_c)
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
   subroutine librpa_set_scf_dimension(this, nspins, nkpts, nstates, nbasis, st_istate, nstates_local, st_ibasis, nbasis_local)
      use iso_c_binding, only: c_associated
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      integer, intent(in) :: nspins, nkpts, nstates, nbasis
      integer, intent(in) :: st_istate, nstates_local, st_ibasis, nbasis_local

      integer(c_int) :: nspins_c, nkpts_c, nstates_c, nbasis_c
      integer(c_int) :: st_istate_c, nstates_local_c, st_ibasis_c, nbasis_local_c

      nspins_c = int(nspins, kind=c_int)
      nkpts_c = int(nkpts, kind=c_int)
      nstates_c = int(nstates, kind=c_int)
      nbasis_c = int(nbasis, kind=c_int)
      st_istate_c = int(st_istate, kind=c_int) - 1
      nstates_local_c = int(nstates_local, kind=c_int)
      st_ibasis_c = int(st_ibasis, kind=c_int) - 1
      nbasis_local_c = int(nbasis_local, kind=c_int)

      call librpa_set_scf_dimension_c(this%ptr_c_handle, nspins_c, nkpts_c, nstates_c, nbasis_c, &
                                      st_istate_c, nstates_local_c, st_ibasis_c, nbasis_local_c)
   end subroutine librpa_set_scf_dimension

   subroutine librpa_set_wg_ekb_efermi(this, nspins, nkpts, nstates, wg, ekb, efermi)
      use iso_c_binding, only: c_associated
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      integer, intent(in) :: nspins, nkpts, nstates
      real(8), intent(in) :: wg(nstates, nkpts, nspins)
      real(8), intent(in) :: ekb(nstates, nkpts, nspins)
      real(8), intent(in) :: efermi

      integer(c_int) :: nspins_c, nkpts_c, nstates_c
      integer :: ispin, ik, istate, idx, ierr
      integer :: n
      real(c_double), allocatable :: wg_c(:), ekb_c(:)
      real(c_double) :: efermi_c

      n = nspins * nkpts * nstates
      allocate(wg_c(n))
      allocate(ekb_c(n))

      idx = 0
      do ispin = 1, nspins
         do ik = 1, nkpts
            do istate = 1, nstates
               idx = idx + 1
               wg_c(idx) = wg(istate, ik, ispin)
               ekb_c(idx) = ekb(istate, ik, ispin)
            end do
         end do
      end do

      nspins_c = int(nspins, kind=c_int)
      nkpts_c = int(nkpts, kind=c_int)
      nstates_c = int(nstates, kind=c_int)
      efermi_c = real(efermi, kind=c_double)

      call librpa_set_wg_ekb_efermi_c(this%ptr_c_handle, nspins_c, nkpts_c, nstates_c, wg_c, ekb_c, efermi_c)
      deallocate(wg_c, ekb_c)
   end subroutine librpa_set_wg_ekb_efermi

   !> @brief Set the wave-function expansion coefficients
   !> @param ispin          spin index (starting from 1) of the wave function
   !> @param ik             k-point index (starting from 1) of the wave function
   !> @param nstates_local  local dimenstion (number of states) of the parsed wave-function
   !> @param nbasis_local   local dimenstion (number of basis functions) of the parsed wave-function
   !> @param wfc_cplx       Complex-valued wave function to parse
   subroutine librpa_set_wfc(this, ispin, ik, nstates_local, nbasis_local, wfc_cplx)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      integer, intent(in) :: ispin, ik, nstates_local, nbasis_local
      complex(kind=8), intent(in) :: wfc_cplx(nbasis_local, nstates_local)

      integer(c_int) :: ispin_c, ik_c, nstates_local_c, nbasis_local_c
      real(c_double), allocatable :: wfc_real(:,:), wfc_imag(:,:)

      ispin_c = int(ispin, kind=c_int) - 1
      ik_c = int(ik, kind=c_int) - 1

      allocate(wfc_real(nbasis_local, nstates_local))
      allocate(wfc_imag(nbasis_local, nstates_local))
      wfc_real(:,:) = real(wfc_cplx, kind=c_double)
      wfc_imag(:,:) = real(aimag(wfc_cplx), kind=c_double)

      call librpa_set_wfc_c(this%ptr_c_handle, ispin, ik, nstates_local, nbasis_local, wfc_real, wfc_imag)
      deallocate(wfc_real, wfc_imag)
   end subroutine librpa_set_wfc

   !> @brief Set the wave-function atomic basis
   !>
   !> @param natoms   number of atoms
   !> @param nbs_aux  number of wave-function basis on each atom
   subroutine librpa_set_ao_basis_wfc(this, natoms, nbs_wfc)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      integer, intent(in) :: natoms
      integer, intent(in) :: nbs_wfc(natoms)

      integer(c_int) :: natoms_c
      integer(c_size_t), allocatable :: nbs_c(:)
      integer :: iatom

      natoms_c = int(natoms, kind=c_int)
      allocate(nbs_c(natoms))

      do iatom = 1, natoms
         nbs_c(iatom) = int(nbs_wfc(iatom), kind=c_size_t)
      end do

      call librpa_set_ao_basis_wfc_c(this%ptr_c_handle, natoms_c, nbs_c)
      deallocate(nbs_c)
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

      integer(c_int) :: natoms_c
      integer(c_size_t), allocatable :: nbs_c(:)
      integer :: iatom

      natoms_c = int(natoms, kind=c_int)
      allocate(nbs_c(natoms))

      do iatom = 1, natoms
         nbs_c(iatom) = int(nbs_aux(iatom), kind=c_size_t)
      end do

      call librpa_set_ao_basis_aux_c(this%ptr_c_handle, natoms_c, nbs_c)
      deallocate(nbs_c)
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
      real(8), dimension(3, 3), intent(in) :: latt, recplatt

      real(c_double), allocatable :: latt_c(:), recplatt_c(:)
      integer :: ir, ic, idx

      allocate(latt_c(9), recplatt_c(9))
      idx = 0
      do ic = 1, 3
         do ir = 1, 3
            idx = idx + 1
            latt_c(idx) = real(latt(ir, ic), kind=c_double)
            recplatt_c(idx) = real(recplatt(ir, ic), kind=c_double)
         end do
      end do
      call librpa_set_latvec_and_G_c(this%ptr_c_handle, latt_c, recplatt_c)
      deallocate(latt_c, recplatt_c)
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
      real(8), dimension(3, natoms), intent(in) :: posi_cart

      call librpa_set_atoms_c(this%ptr_c_handle, int(natoms, kind=c_int), &
                              int(types, kind=c_int), real(posi_cart, kind=c_double))
   end subroutine librpa_set_atoms

   subroutine librpa_set_kgrids_kvec(this, nk1, nk2, nk3, kvecs)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      integer, intent(in) :: nk1, nk2, nk3
      real(8), intent(in) :: kvecs(3, nk1*nk2*nk3)

      call librpa_set_kgrids_kvec_c(this%ptr_c_handle, &
         int(nk1, kind=c_int), int(nk2, kind=c_int), int(nk3, kind=c_int), real(kvecs, kind=c_double))
   end subroutine librpa_set_kgrids_kvec

   !> @brief Set the mapping from full k-point list to the irreducbile sector
   !>
   !> @param nkpts     number of k-points in the full Brillouin zone
   !> @param map_ibzk  mapping to the k-point in the irreducible sector
   !>
   !>
   !> Example: four-k-point case where the first two and last points are in the irreducbile sector,
   !>          and the third point is mapped to the second, then map_ibzk should be (1, 2, 2, 4)
   subroutine librpa_set_ibz_mapping(this, nkpts, map_ibzk)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      integer, intent(in) :: nkpts
      integer, dimension(nkpts), intent(in) :: map_ibzk

      integer :: ik
      integer(c_int), allocatable :: map_ibzk_c(:)

      allocate(map_ibzk_c(nkpts))
      do ik = 1, nkpts
         map_ibzk_c(ik) = int(map_ibzk(ik), kind=c_int) - 1
      end do
      call librpa_set_ibz_mapping_c(this%ptr_c_handle, int(nkpts, kind=c_int), map_ibzk_c)
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
      real(8), intent(in) :: coeff(naux_i, nao_j, nao_i)

      ! Check if the routing is a valid LIBRPA_ROUTING_ parameter
      select case (routing)
         case (LIBRPA_ROUTING_AUTO)
         case (LIBRPA_ROUTING_RTAU)
         case (LIBRPA_ROUTING_ATOMPAIR)
         case (LIBRPA_ROUTING_LIBRI)
         case default
            write(*,*) "Invalid routing parameter:", routing
            stop
      end select

      call librpa_set_lri_coeff_c(this%ptr_c_handle, &
         int(routing, kind=c_int), int(i_atom-1, kind=c_int), int(j_atom-1, kind=c_int), &
         int(nao_i, kind=c_int), int(nao_j, kind=c_int), int(naux_i, kind=c_int), &
         int(r, kind=c_int), real(coeff, kind=c_double) &
         )
   end subroutine librpa_set_lri_coeff

   subroutine librpa_set_aux_bare_coulomb_k_atom_pair &
         (this, ik, i_atom, j_atom, naux_i, naux_j, vq, vq_threshold)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      integer, intent(in) :: ik, i_atom, j_atom, naux_i, naux_j
      complex(8), intent(in) :: vq(naux_i, naux_j)
      real(8), intent(in) :: vq_threshold

      call librpa_set_aux_bare_coulomb_k_atom_pair_c(this%ptr_c_handle, &
         int(ik-1, kind=c_int), int(i_atom-1, kind=c_int), int(j_atom-1, kind=c_int), &
         int(naux_i, kind=c_int), int(naux_j, kind=c_int), &
         transpose(real(vq, kind=c_double)), &
         transpose(real(aimag(vq), kind=c_double)), &
         real(vq_threshold, kind=c_double))
   end subroutine librpa_set_aux_bare_coulomb_k_atom_pair

   subroutine librpa_set_aux_cut_coulomb_k_atom_pair &
         (this, ik, i_atom, j_atom, naux_i, naux_j, vq, vq_threshold)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      integer, intent(in) :: ik, i_atom, j_atom, naux_i, naux_j
      complex(8), intent(in) :: vq(naux_i, naux_j)
      real(8), intent(in) :: vq_threshold

      call librpa_set_aux_cut_coulomb_k_atom_pair_c(this%ptr_c_handle, &
         int(ik-1, kind=c_int), int(i_atom-1, kind=c_int), int(j_atom-1, kind=c_int), &
         int(naux_i, kind=c_int), int(naux_j, kind=c_int), &
         transpose(real(vq, kind=c_double)), &
         transpose(real(aimag(vq), kind=c_double)), &
         real(vq_threshold, kind=c_double))
   end subroutine librpa_set_aux_cut_coulomb_k_atom_pair

   subroutine librpa_set_aux_bare_coulomb_k_2d_block &
         (this, ik, mu_begin, mu_end, nu_begin, nu_end, vq)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      integer, intent(in) :: ik, mu_begin, mu_end, nu_begin, nu_end
      complex(8), intent(in) :: vq(mu_end-mu_begin+1, nu_end-nu_begin+1)

      call librpa_set_aux_bare_coulomb_k_2d_block_c(this%ptr_c_handle, &
         int(ik-1, kind=c_int), int(mu_begin-1, kind=c_int), int(mu_end-1, kind=c_int), &
         int(nu_begin-1, kind=c_int), int(nu_end-1, kind=c_int), &
         transpose(real(vq, kind=c_double)), &
         transpose(real(aimag(vq), kind=c_double)))
   end subroutine librpa_set_aux_bare_coulomb_k_2d_block

   subroutine librpa_set_aux_cut_coulomb_k_2d_block &
         (this, ik, mu_begin, mu_end, nu_begin, nu_end, vq)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      integer, intent(in) :: ik, mu_begin, mu_end, nu_begin, nu_end
      complex(8), intent(in) :: vq(mu_end-mu_begin+1, nu_end-nu_begin+1)

      call librpa_set_aux_cut_coulomb_k_2d_block_c(this%ptr_c_handle, &
         int(ik-1, kind=c_int), int(mu_begin-1, kind=c_int), int(mu_end-1, kind=c_int), &
         int(nu_begin-1, kind=c_int), int(nu_end-1, kind=c_int), &
         transpose(real(vq, kind=c_double)), &
         transpose(real(aimag(vq), kind=c_double)))
   end subroutine librpa_set_aux_cut_coulomb_k_2d_block

   real(8) function librpa_get_rpa_correlation_energy(this, opts, nkpts_ibz, contrib_ibzk)
      class(LibrpaHandler), intent(inout) :: this
      type(LibrpaOptions), intent(inout) :: opts
      integer, intent(in) :: nkpts_ibz
      complex(8), dimension(nkpts_ibz), intent(inout) :: contrib_ibzk

      real(8), allocatable :: contrib_ibzk_re(:), contrib_ibzk_im(:)

      call sync_opts(opts, SYNC_OPTS_F2C)
      allocate(contrib_ibzk_re(nkpts_ibz), contrib_ibzk_im(nkpts_ibz))

      librpa_get_rpa_correlation_energy = &
         librpa_get_rpa_correlation_energy_c( &
            this%ptr_c_handle, opts%opts_c, nkpts_ibz, contrib_ibzk_re, contrib_ibzk_im)
      contrib_ibzk = contrib_ibzk_re + contrib_ibzk_im * CIMAG
      deallocate(contrib_ibzk_re, contrib_ibzk_im)
   end function librpa_get_rpa_correlation_energy

end module librpa_f03
