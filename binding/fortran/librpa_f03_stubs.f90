module librpa_f03
   !=======================================================================
   ! Stub Fortran module for using LibRPA
   !=======================================================================
   implicit none

   private

   !=======================================================================
   ! Public types, constants, and functions
   !=======================================================================
   public :: LibrpaOptions
   public :: LibrpaHandler

   integer, parameter, public :: dp = 8

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

   ! High-level Fortran wrapper
   type :: LibrpaOptions
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

   ! High-level Fortran wrapper
   type :: LibrpaHandler
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
         ! Compute
         procedure :: get_imaginary_frequency_grids => librpa_get_imaginary_frequency_grids
         procedure :: get_rpa_correlation_energy => librpa_get_rpa_correlation_energy
         procedure :: build_exx => librpa_build_exx
         procedure :: get_exx_pot_kgrid => librpa_get_exx_pot_kgrid
         procedure :: build_g0w0_sigma => librpa_build_g0w0_sigma
         procedure :: get_g0w0_qpe_kgrid => librpa_get_g0w0_qpe_kgrid
   end type LibrpaHandler

contains

   ! Can be customized by actual hosting code
   subroutine error_on_call(func)
      implicit none
      character(len=*), intent(in) :: func
      character(len=200) :: info_str

      write(info_str,'(1X,A,A,A)') '* You have called a librpa_f03_stub routine: ', &
         trim(func), '. Make sure you have linked LibRPA library.'
      write(*,'(A)') info_str
      stop
   end subroutine error_on_call

   !=======================================================================
   ! Usually no need change things below
   !=======================================================================

   subroutine librpa_init_options(opts)
      implicit none
      class(LibrpaOptions), intent(inout) :: opts
      call error_on_call("librpa_init_options")
   end subroutine librpa_init_options

   subroutine librpa_set_output_dir(opts, output_dir)
      implicit none
      class(LibrpaOptions), intent(inout) :: opts
      character(len=*), intent(in) :: output_dir
      call error_on_call("librpa_set_output_dir")
   end subroutine librpa_set_output_dir

   subroutine librpa_init_global(sw_redirect, redirect_path, sw_process)
      implicit none
      logical, intent(in), optional :: sw_redirect, sw_process
      character(len=*), intent(in), optional :: redirect_path
      call error_on_call("librpa_init_global")
   end subroutine librpa_init_global

   subroutine librpa_finalize_global()
      implicit none
      call error_on_call("librpa_finalize_global")
   end subroutine librpa_finalize_global

   integer function librpa_get_major_version() result(v)
      v = -1
      call error_on_call("librpa_get_major_version")
   end function librpa_get_major_version

   integer function librpa_get_minor_version() result(v)
      v = -1
      call error_on_call("librpa_get_minor_version")
   end function librpa_get_minor_version

   integer function librpa_get_patch_version() result(v)
      v = -1
      call error_on_call("librpa_get_patch_version")
   end function librpa_get_patch_version

   subroutine librpa_test()
      implicit none
      call error_on_call("librpa_test")
   end subroutine librpa_test

   subroutine librpa_create_handler(this, comm)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      integer, intent(in) :: comm
      call error_on_call("librpa_create_handler")
   end subroutine librpa_create_handler

   subroutine librpa_destroy_handler(this)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      call error_on_call("librpa_destroy_handler")
   end subroutine librpa_destroy_handler

   ! Input functions
   subroutine librpa_set_scf_dimension(this, nspins, nkpts, nstates, nbasis, st_istate, nstates_local, st_ibasis, nbasis_local)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      integer, intent(in) :: nspins, nkpts, nstates, nbasis
      integer, intent(in) :: st_istate, nstates_local, st_ibasis, nbasis_local
      call error_on_call("librpa_set_scf_dimension")
   end subroutine librpa_set_scf_dimension

   subroutine librpa_set_wg_ekb_efermi(this, nspins, nkpts, nstates, wg, ekb, efermi)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      integer, intent(in) :: nspins, nkpts, nstates
      real(dp), intent(in) :: wg(nstates, nkpts, nspins)
      real(dp), intent(in) :: ekb(nstates, nkpts, nspins)
      real(dp), intent(in) :: efermi
      call error_on_call("librpa_set_wg_ekb_efermi")
   end subroutine librpa_set_wg_ekb_efermi

   subroutine librpa_set_wfc(this, ispin, ik, nstates_local, nbasis_local, wfc_cplx)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      integer, intent(in) :: ispin, ik, nstates_local, nbasis_local
      complex(dp), intent(in) :: wfc_cplx(nbasis_local, nstates_local)
      call error_on_call("librpa_set_wfc")
   end subroutine librpa_set_wfc

   subroutine librpa_set_ao_basis_wfc(this, natoms, nbs_wfc)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      integer, intent(in) :: natoms
      integer, intent(in) :: nbs_wfc(natoms)
      call error_on_call("librpa_set_ao_basis_wfc")
   end subroutine librpa_set_ao_basis_wfc

   subroutine librpa_set_ao_basis_aux(this, natoms, nbs_aux)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      integer, intent(in) :: natoms
      integer, intent(in) :: nbs_aux(natoms)
      call error_on_call("librpa_set_ao_basis_aux")
   end subroutine librpa_set_ao_basis_aux

   subroutine librpa_set_latvec_and_G(this, latt, recplatt)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      real(dp), dimension(3, 3), intent(in) :: latt, recplatt
      call error_on_call("librpa_set_latvec_and_G")
   end subroutine librpa_set_latvec_and_G

   subroutine librpa_set_atoms(this, natoms, types, posi_cart)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      integer, intent(in) :: natoms
      integer, dimension(natoms), intent(in) :: types
      real(dp), dimension(3, natoms), intent(in) :: posi_cart
      call error_on_call("librpa_set_atoms")
   end subroutine librpa_set_atoms

   subroutine librpa_set_kgrids_kvec(this, nk1, nk2, nk3, kvecs)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      integer, intent(in) :: nk1, nk2, nk3
      real(dp), intent(in) :: kvecs(3, nk1*nk2*nk3)
      call error_on_call("librpa_set_kgrids_kvec")
   end subroutine librpa_set_kgrids_kvec

   subroutine librpa_set_ibz_mapping(this, nkpts, map_ibzk)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      integer, intent(in) :: nkpts
      integer, dimension(nkpts), intent(in) :: map_ibzk
      call error_on_call("librpa_set_ibz_mapping")
   end subroutine librpa_set_ibz_mapping

   subroutine librpa_set_lri_coeff(this, routing, i_atom, j_atom, nao_i, nao_j, naux_i, r, coeff)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      integer, intent(in) :: routing, i_atom, j_atom, nao_i, nao_j, naux_i
      integer, dimension(3), intent(in) :: r
      real(dp), intent(in) :: coeff(naux_i, nao_j, nao_i)
      call error_on_call("librpa_set_lri_coeff")
   end subroutine librpa_set_lri_coeff

   subroutine librpa_set_aux_bare_coulomb_k_atom_pair &
         (this, ik, i_atom, j_atom, naux_i, naux_j, vq, vq_threshold)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      integer, intent(in) :: ik, i_atom, j_atom, naux_i, naux_j
      complex(dp), intent(in) :: vq(naux_i, naux_j)
      real(dp), intent(in) :: vq_threshold
      call error_on_call("librpa_set_aux_bare_coulomb_k_atom_pair")
   end subroutine librpa_set_aux_bare_coulomb_k_atom_pair

   subroutine librpa_set_aux_cut_coulomb_k_atom_pair &
         (this, ik, i_atom, j_atom, naux_i, naux_j, vq, vq_threshold)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      integer, intent(in) :: ik, i_atom, j_atom, naux_i, naux_j
      complex(dp), intent(in) :: vq(naux_i, naux_j)
      real(dp), intent(in) :: vq_threshold
      call error_on_call("librpa_set_aux_cut_coulomb_k_atom_pair")
   end subroutine librpa_set_aux_cut_coulomb_k_atom_pair

   subroutine librpa_set_aux_bare_coulomb_k_2d_block &
         (this, ik, mu_begin, mu_end, nu_begin, nu_end, vq)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      integer, intent(in) :: ik, mu_begin, mu_end, nu_begin, nu_end
      complex(dp), intent(in) :: vq(mu_end-mu_begin+1, nu_end-nu_begin+1)
      call error_on_call("librpa_set_aux_bare_coulomb_k_2d_block")
   end subroutine librpa_set_aux_bare_coulomb_k_2d_block

   subroutine librpa_set_aux_cut_coulomb_k_2d_block &
         (this, ik, mu_begin, mu_end, nu_begin, nu_end, vq)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      integer, intent(in) :: ik, mu_begin, mu_end, nu_begin, nu_end
      complex(dp), intent(in) :: vq(mu_end-mu_begin+1, nu_end-nu_begin+1)
      call error_on_call("librpa_set_aux_cut_coulomb_k_2d_block")
   end subroutine librpa_set_aux_cut_coulomb_k_2d_block

   subroutine librpa_set_dielect_func_imagfreq(this, nfreq, omegas_imag, dielect_func)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      integer, intent(in) :: nfreq
      real(dp), dimension(nfreq), intent(in) :: omegas_imag
      real(dp), dimension(nfreq), intent(in) :: dielect_func
      call error_on_call("librpa_set_dielect_func_imagfreq")
   end subroutine librpa_set_dielect_func_imagfreq

   ! Compute functions
   subroutine librpa_get_imaginary_frequency_grids(this, opts, omegas, weights)
      class(LibrpaHandler), intent(inout) :: this
      type(LibrpaOptions), intent(inout) :: opts
      real(dp), allocatable, intent(inout) :: omegas(:), weights(:)
      call error_on_call("librpa_get_imaginary_frequency_grids")
   end subroutine librpa_get_imaginary_frequency_grids

   real(dp) function librpa_get_rpa_correlation_energy(this, opts, nkpts_ibz, contrib_ibzk) result(e)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      type(LibrpaOptions), intent(inout) :: opts
      integer, intent(in) :: nkpts_ibz
      complex(dp), dimension(nkpts_ibz), intent(inout) :: contrib_ibzk

      e = 0.0d0
      call error_on_call("librpa_get_rpa_correlation_energy")
   end function librpa_get_rpa_correlation_energy

   subroutine librpa_build_exx(this, opts)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      type(LibrpaOptions), intent(inout) :: opts
      call error_on_call("librpa_build_exx")
   end subroutine librpa_build_exx

   subroutine librpa_get_exx_pot_kgrid(this, opts, n_spins, n_kpoints_local, iks_local, &
                                       i_state_low, i_state_high, vexx)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      type(LibrpaOptions), intent(inout) :: opts
      integer, dimension(n_kpoints_local), intent(in) :: iks_local
      integer, intent(in) :: n_spins, n_kpoints_local, i_state_low, i_state_high
      real(dp), dimension(i_state_high - i_state_low + 1, n_kpoints_local, n_spins), intent(inout) :: vexx
      call error_on_call("librpa_get_exx_pot_kgrid")
   end subroutine librpa_get_exx_pot_kgrid

   subroutine librpa_build_g0w0_sigma(this, opts)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      type(LibrpaOptions), intent(inout) :: opts
      call error_on_call("librpa_build_g0w0_sigma")
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
      call error_on_call("librpa_get_g0w0_qpe_kgrid")
   end subroutine librpa_get_g0w0_qpe_kgrid

end module librpa_f03
