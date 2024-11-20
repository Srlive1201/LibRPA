module librpa

   use iso_c_binding, only: c_char, c_int, c_double
   implicit none

   private
   public :: LibRPAParams
   public :: initialize_librpa_environment
   public :: finalize_librpa_environment
   public :: get_rpa_correlation_energy
   public :: set_dimension
   public :: set_wg_ekb_efermi
   public :: set_ao_basis_wfc
   public :: set_latvec_and_G
   public :: set_kgrids_kvec_tot
   public :: set_ibz2bz_index_and_weight
   public :: set_ao_basis_aux
   public :: set_aux_bare_coulomb_k_2D_block
   public :: set_aux_cut_coulomb_k_atom_pair
   public :: set_aux_bare_coulomb_k_atom_pair
   public :: set_aux_cut_coulomb_k_2D_block
   public :: set_librpa_params
   public :: get_default_librpa_params
   public :: run_librpa_main

   !!! Interface declaration

   !> @brief User interface for control parameters in LibRPA
   type :: LibRPAParams
      character(len=100) :: output_file
      character(len=100) :: output_dir
      character(len=20)  :: parallel_routing
      character(len=20)  :: tfgrids_type
      character(len=20)  :: DFT_software

      integer :: nfreq

      logical :: debug
      logical :: use_scalapack_ecrpa

      real*8 :: gf_R_threshold
      real*8 :: cs_threshold
      real*8 :: vq_threshold
      real*8 :: sqrt_coulomb_threshold
      real*8 :: libri_chi0_threshold_C
      real*8 :: libri_chi0_threshold_G
      real*8 :: libri_exx_threshold_CSM
      real*8 :: libri_exx_threshold_C
      real*8 :: libri_exx_threshold_D
      real*8 :: libri_exx_threshold_V
      real*8 :: libri_gw_threshold_C
      real*8 :: libri_gw_threshold_G
      real*8 :: libri_gw_threshold_W
   end type LibRPAParams

   !> @brief C interface to the LibRPAParams struct
   type, bind(c) :: LibRPAParams_c
      character(kind=c_char, len=1) :: output_file(100)
      character(kind=c_char, len=1) :: output_dir(100)
      character(kind=c_char, len=1) :: parallel_routing(20)
      character(kind=c_char, len=1) :: tfgrids_type(20)
      character(kind=c_char, len=1) :: DFT_software(20)

      integer(c_int) :: nfreq

      integer(c_int) :: debug
      integer(c_int) :: use_scalapack_ecrpa

      real(c_double) :: gf_R_threshold
      real(c_double) :: cs_threshold
      real(c_double) :: vq_threshold
      real(c_double) :: sqrt_coulomb_threshold
      real(c_double) :: libri_chi0_threshold_C
      real(c_double) :: libri_chi0_threshold_G
      real(c_double) :: libri_exx_threshold_CSM
      real(c_double) :: libri_exx_threshold_C
      real(c_double) :: libri_exx_threshold_D
      real(c_double) :: libri_exx_threshold_V
      real(c_double) :: libri_gw_threshold_C
      real(c_double) :: libri_gw_threshold_G
      real(c_double) :: libri_gw_threshold_W
   end type LibRPAParams_c

   interface
      subroutine initialize_librpa_environment( &
         comm_in, is_fortran_comm, redirect_stdout, output_filename) bind(c, name="initialize_librpa_environment")
         use, intrinsic :: iso_c_binding, only: c_int, c_char
         integer(c_int), value :: comm_in, is_fortran_comm, redirect_stdout
         character(kind=c_char), dimension(*), intent(in) :: output_filename
      end subroutine
   end interface

   interface
      subroutine finalize_librpa_environment() bind(c, name="finalize_librpa_environment")
      end subroutine
   end interface

   interface
      subroutine set_dimension(nspins, nkpts, nstates, nbasis, natoms) bind(c, name="set_dimension")
         use, intrinsic :: iso_c_binding
         integer(c_int), value :: nspins, nkpts, nstates, nbasis, natoms
      end subroutine
   end interface

   interface
      subroutine set_wg_ekb_efermi(nspins, nkpts, nstates, wg, ekb, efermi) bind(c, name="set_wg_ekb_efermi")
         use, intrinsic :: iso_c_binding
         integer(c_int), value :: nspins, nkpts, nstates
         real(c_double), dimension(*), intent(inout) :: wg
         real(c_double), dimension(*), intent(inout) :: ekb
         real(c_double), value :: efermi
      end subroutine
   end interface

   interface
      subroutine set_ao_basis_wfc(nspins, nkpts, wfc_real, wfc_imag) bind(c, name="set_ao_basis_wfc")
         use, intrinsic :: iso_c_binding
         integer(c_int), value :: nspins, nkpts
         real(c_double), dimension(*), intent(inout) :: wfc_real
         real(c_double), dimension(*), intent(inout) :: wfc_imag
      end subroutine
   end interface

   interface
      subroutine set_latvec_and_G(lat_mat, G_mat) bind(c, name="set_latvec_and_G")
         use, intrinsic :: iso_c_binding
         real(c_double), dimension(*), intent(inout) :: lat_mat
         real(c_double), dimension(*), intent(inout) :: G_mat
      end subroutine
   end interface

   interface
      subroutine set_kgrids_kvec_tot(nkx, nky, nkz, kvecs) bind(c, name="set_kgrids_kvec_tot")
         use, intrinsic :: iso_c_binding
         integer(c_int), value :: nkx, nky, nkz
         real(c_double), dimension(*), intent(inout) :: kvecs
      end subroutine
   end interface

   interface
      subroutine set_ibz2bz_index_and_weight(nk_irk, ibz2bz_index, wk_irk) bind(c, name="set_ibz2bz_index_and_weight")
         use, intrinsic :: iso_c_binding
         integer(c_int), value :: nk_irk
         integer(c_int), dimension(*), intent(inout) :: ibz2bz_index
         real(c_double), dimension(*), intent(inout) :: wk_irk
      end subroutine
   end interface

   interface
      subroutine set_ao_basis_aux(I, J, nbasis_i, nbasis_j, naux_mu, R, Cs_in, insert_index_only) &
            bind(c, name="set_ao_basis_aux")
         use, intrinsic :: iso_c_binding
         integer(c_int), value :: I, J, nbasis_i, nbasis_j, naux_mu
         integer(c_int), dimension(3), intent(inout) :: R
         real(c_double), dimension(*), intent(inout) :: Cs_in
         integer(c_int), value :: insert_index_only
      end subroutine
   end interface

   interface
      subroutine set_aux_bare_coulomb_k_atom_pair(ik, I, J, naux_mu, naux_nu, Vq_real_in, Vq_imag_in) &
            bind(c, name="set_aux_bare_coulomb_k_atom_pair")
         use, intrinsic :: iso_c_binding
         integer(c_int), value :: I, J, naux_mu, naux_nu, ik
         real(c_double), dimension(*), intent(inout) :: Vq_real_in
         real(c_double), dimension(*), intent(inout) :: Vq_imag_in
      end subroutine
   end interface

   interface
      subroutine set_aux_cut_coulomb_k_atom_pair(ik, I, J, naux_mu, naux_nu, Vq_real_in, Vq_imag_in) &
            bind(c, name="set_aux_cut_coulomb_k_atom_pair")
         use, intrinsic :: iso_c_binding
         integer(c_int), value :: I, J, naux_mu, naux_nu, ik
         real(c_double), dimension(*), intent(inout) :: Vq_real_in
         real(c_double), dimension(*), intent(inout) :: Vq_imag_in
      end subroutine
   end interface

   interface
      subroutine set_aux_bare_coulomb_k_2D_block(ik, max_naux, mu_begin, mu_end, nu_begin, nu_end, Vq_real_in, Vq_imag_in) &
            bind(c, name="set_aux_bare_coulomb_k_2D_block")
         use, intrinsic :: iso_c_binding
         integer(c_int), value :: ik, max_naux, mu_begin, mu_end, nu_begin, nu_end
         real(c_double), dimension(*), intent(inout) :: Vq_real_in
         real(c_double), dimension(*), intent(inout) :: Vq_imag_in
      end subroutine
   end interface

   interface
      subroutine set_aux_cut_coulomb_k_2D_block(ik, max_naux, mu_begin, mu_end, nu_begin, nu_end, Vq_real_in, Vq_imag_in) &
            bind(c, name="set_aux_cut_coulomb_k_2D_block")
         use, intrinsic :: iso_c_binding
         integer(c_int), value :: ik, max_naux, mu_begin, mu_end, nu_begin, nu_end
         real(c_double), dimension(*), intent(inout) :: Vq_real_in
         real(c_double), dimension(*), intent(inout) :: Vq_imag_in
      end subroutine
   end interface

   interface
      subroutine run_librpa_main() bind(c, name="run_librpa_main")
      end subroutine
   end interface

   interface
      subroutine set_librpa_params_c(params_c) bind(c, name="set_librpa_params")
         use, intrinsic :: iso_c_binding
         import :: LibRPAParams_c
         type(LibRPAParams_c), intent(in) :: params_c
      end subroutine
   end interface

   interface
      subroutine get_default_librpa_params_c(params_c) bind(c, name="get_default_librpa_params")
         use, intrinsic :: iso_c_binding
         import :: LibRPAParams_c
         type(LibRPAParams_c), intent(inout) :: params_c
      end subroutine
   end interface

   interface
      subroutine get_rpa_correlation_energy(rpa_corr, rpa_corr_irk_contrib) bind(c, name="get_rpa_correlation_energy")
         import :: c_double
         real(c_double), dimension(2), intent(out) :: rpa_corr
         real(c_double), dimension(*), intent(out) :: rpa_corr_irk_contrib
      end subroutine
   end interface

contains

   !> @brief Copy a Fortran character varaible to C char array
   !!
   !> adapted from https://fortranwiki.org/fortran/show/c_interface_module
   subroutine f_c_string_chars(f_string, c_string, c_string_len, trim_f)
      use iso_c_binding, only: c_null_char
      implicit none

      character(len=*), intent(in) :: f_string
      character(len=1, kind=c_char), dimension(*), intent(out) :: c_string
      ! Max string length, INCLUDING THE TERMINAL NUL
      integer, intent(in), optional :: c_string_len
      logical, intent(in), optional :: trim_f

      integer :: i, strlen

      if (present(trim_f) .and. trim_f) then
         strlen = len(trim(f_string))
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
         c_bi = 1
      else
         c_bi = 0
      endif

   end subroutine f_c_bool


   !> @brief Convert C integer as boolean to Fortran logical
   subroutine c_f_bool(c_bi, f_logical)
      implicit none

      integer(kind=c_int), intent(in) :: c_bi
      logical, intent(out) :: f_logical

      f_logical = (c_bi .ne. 0)

   end subroutine c_f_bool


   !> @brief helper function to convert data between Fortran and C derived type
   subroutine communicate_librpa_params(params, params_c, f2c)
      implicit none

      type(LibRPAParams), intent(inout)   :: params
      type(LibRPAParams_c), intent(inout) :: params_c
      logical, intent(in) :: f2c

      if (f2c) then
         ! Copy Fortran parameters to C
         call f_c_string_chars(params%output_file      , params_c%output_file     , trim_f=.true.)
         call f_c_string_chars(params%output_dir       , params_c%output_dir      , trim_f=.true.)
         call f_c_string_chars(params%parallel_routing , params_c%parallel_routing, trim_f=.true.)
         call f_c_string_chars(params%tfgrids_type     , params_c%tfgrids_type    , trim_f=.true.)
         call f_c_string_chars(params%DFT_software     , params_c%DFT_software    , trim_f=.true.)

         params_c%nfreq = params%nfreq

         call f_c_bool(params%debug,               params_c%debug)
         call f_c_bool(params%use_scalapack_ecrpa, params_c%use_scalapack_ecrpa)

         params_c%gf_R_threshold           = params%gf_R_threshold
         params_c%cs_threshold             = params%cs_threshold
         params_c%vq_threshold             = params%vq_threshold
         params_c%sqrt_coulomb_threshold   = params%sqrt_coulomb_threshold
         params_c%libri_chi0_threshold_C   = params%libri_chi0_threshold_C
         params_c%libri_chi0_threshold_G   = params%libri_chi0_threshold_G
         params_c%libri_exx_threshold_CSM  = params%libri_exx_threshold_CSM
         params_c%libri_exx_threshold_C    = params%libri_exx_threshold_C
         params_c%libri_exx_threshold_D    = params%libri_exx_threshold_D
         params_c%libri_exx_threshold_V    = params%libri_exx_threshold_V
         params_c%libri_gw_threshold_C     = params%libri_gw_threshold_C
         params_c%libri_gw_threshold_G     = params%libri_gw_threshold_G
         params_c%libri_gw_threshold_W     = params%libri_gw_threshold_W
      else
         ! Copy C parameters to Fortran
         call c_f_string_chars(params_c%output_file      , params%output_file)
         call c_f_string_chars(params_c%output_dir       , params%output_dir)
         call c_f_string_chars(params_c%parallel_routing , params%parallel_routing)
         call c_f_string_chars(params_c%tfgrids_type     , params%tfgrids_type)
         call c_f_string_chars(params_c%DFT_software     , params%DFT_software)

         params%nfreq = params_c%nfreq

         call c_f_bool(params_c%debug,               params%debug)
         call c_f_bool(params_c%use_scalapack_ecrpa, params%use_scalapack_ecrpa)

         params%gf_R_threshold           = params_c%gf_R_threshold
         params%cs_threshold             = params_c%cs_threshold
         params%vq_threshold             = params_c%vq_threshold
         params%sqrt_coulomb_threshold   = params_c%sqrt_coulomb_threshold
         params%libri_chi0_threshold_C   = params_c%libri_chi0_threshold_C
         params%libri_chi0_threshold_G   = params_c%libri_chi0_threshold_G
         params%libri_exx_threshold_CSM  = params_c%libri_exx_threshold_CSM
         params%libri_exx_threshold_C    = params_c%libri_exx_threshold_C
         params%libri_exx_threshold_D    = params_c%libri_exx_threshold_D
         params%libri_exx_threshold_V    = params_c%libri_exx_threshold_V
         params%libri_gw_threshold_C     = params_c%libri_gw_threshold_C
         params%libri_gw_threshold_G     = params_c%libri_gw_threshold_G
         params%libri_gw_threshold_W     = params_c%libri_gw_threshold_W
      endif

   end subroutine communicate_librpa_params


   subroutine set_librpa_params(params)

      implicit none
      type(LibRPAParams), intent(inout) :: params

      type(LibRPAParams_c) :: params_c

      call communicate_librpa_params(params, params_c, .true.)
      call set_librpa_params_c(params_c)

   end subroutine set_librpa_params


   subroutine get_default_librpa_params(params)

      implicit none
      type(LibRPAParams), intent(inout) :: params

      type(LibRPAParams_c) :: params_c

      call get_default_librpa_params_c(params_c)
      call communicate_librpa_params(params, params_c, .false.)

   end subroutine get_default_librpa_params


end module
