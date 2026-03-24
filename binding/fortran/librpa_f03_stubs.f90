!> @file librpa_f03.f90
!> @brief Fortran 2003 binding for LibRPA
!>
!> This module provides Fortran interfaces to LibRPA functionality for
!> performing RPA correlation energy, exact exchange, and G0W0 calculations.
!>
!> ## Usage
!>
!> Typical workflow:
!> @code{.f90}
!> use librpa_f03
!> implicit none
!>
!> type(LibrpaOptions) :: opts
!> type(LibrpaHandler) :: h
!>
!> ! Initialize LibRPA environment
!> call librpa_init_global()
!>
!> ! Initialize options
!> call opts%init()
!> opts%output_level = LIBRPA_VERBOSE_INFO
!>
!> ! Create handler
!> call h%init(MPI_COMM_WORLD)
!>
!> ! Set input data
!> call h%set_scf_dimension(nspins, nkpts, nstates, nbasis)
!> ! ... set more input data ...
!>
!> ! Perform calculation
!> Ec = h%get_rpa_correlation_energy(h, opts, nkpts_ibz, contrib_ibzk(:))
!>
!> ! Clean up
!> call h%free()
!> call librpa_finalize_global()
!> @endcode

!> @brief Fortran 2003 module for LibRPA API
module librpa_f03
   implicit none

   private

   !=======================================================================
   ! Public types, constants, and functions
   !=======================================================================
   public :: LibrpaOptions
   public :: LibrpaHandler

   !> @brief Double precision kind for user data.
   !>
   integer, parameter, public :: dp = 8

   !> @brief Maximum length for string parameters.
   integer, parameter, public :: LIBRPA_MAX_STRLEN = 200

   !> @name Verbosity levels
   !> @brief Controls the amount of output during computation.
   !> @{
   integer, parameter, public :: LIBRPA_VERBOSE_DEBUG = 4      !< Debug output
   integer, parameter, public :: LIBRPA_VERBOSE_WARN = 3       !< Warnings and above
   integer, parameter, public :: LIBRPA_VERBOSE_INFO = 2       !< Informational messages
   integer, parameter, public :: LIBRPA_VERBOSE_CRITICAL = 1   !< Critical errors only
   integer, parameter, public :: LIBRPA_VERBOSE_SILENT = 0     !< No output
   !> @}

   !> @brief Undefined or unset value for integer parameters.
   integer, parameter :: LIBRPA_UNSET = -101

   !> @brief Automatic selection value. LibRPA will choose appropriate setting.
   integer, parameter :: LIBRPA_AUTO = -51

   !> @name Parallel routing strategies
   !> @brief Specifies how computation is distributed across MPI processes.
   !> @{
   integer, parameter, public :: LIBRPA_ROUTING_UNSET = LIBRPA_UNSET  !< Unset
   integer, parameter, public :: LIBRPA_ROUTING_AUTO = LIBRPA_AUTO    !< Auto-select
   integer, parameter, public :: LIBRPA_ROUTING_RTAU = 0              !< Real-space tau decomposition
   integer, parameter, public :: LIBRPA_ROUTING_ATOMPAIR = 1          !< Atom-pair parallelization
   integer, parameter, public :: LIBRPA_ROUTING_LIBRI = 2             !< Use LibRI for RI basis
   !> @}

   !> @name Time-frequency grid types
   !> @brief Different grid types for numerical integration.
   !> @{
   integer, parameter, public :: LIBRPA_TFGRID_UNSET = LIBRPA_UNSET     !< Unset
   integer, parameter, public :: LIBRPA_TFGRID_GL = 0                   !< Gauss-Legendre
   integer, parameter, public :: LIBRPA_TFGRID_GCI= 1                   !< Gauss-Chebyshev type I
   integer, parameter, public :: LIBRPA_TFGRID_GCII = 2                 !< Gauss-Chebyshev type II
   integer, parameter, public :: LIBRPA_TFGRID_MINIMAX = 3              !< Minimax grid
   integer, parameter, public :: LIBRPA_TFGRID_EVENSPACED = 4           !< Evenly spaced
   integer, parameter, public :: LIBRPA_TFGRID_EVENSPACED_TF = 5        !< Evenly spaced in time-frequency
   !> @}

   public :: librpa_init_global
   public :: librpa_finalize_global
   public :: librpa_get_major_version
   public :: librpa_get_minor_version
   public :: librpa_get_patch_version
   public :: librpa_test
   public :: librpa_print_profile

   !> @brief High-level Fortran wrapper for runtime options.
   !>
   !> This type provides a Fortran-friendly interface to LibRPA options.
   !> Initialize with init() method or call librpa_init_options() before use.
   !>
   !> @note The data layout must match the C struct. Do not add or remove members.
   type :: LibrpaOptions
      character(len=LIBRPA_MAX_STRLEN) :: output_dir
      integer :: parallel_routing
      integer :: output_level
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

   !> @brief High-level Fortran wrapper for LibRPA handler.
   !>
   !> This type encapsulates the LibRPA handler and provides member procedures
   !> for setting input data and performing calculations.
   !>
   !> Usage:
   !> @code{.f90}
   !> type(LibrpaHandler) :: h
   !> call h%create(MPI_COMM_WORLD)
   !> call h%set_scf_dimension(nspins, nkpts, nstates, nbasis)
   !> ! ... set more input ...
   !> Ec = h%get_rpa_correlation_energy(opts)
   !> call h%destroy()
   !> @endcode
   type :: LibrpaHandler
      contains
         ! Initialization and destruction
         procedure :: init => librpa_create_handler
         procedure :: free => librpa_destroy_handler
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
         procedure :: get_g0w0_sigc_kgrid => librpa_get_g0w0_sigc_kgrid
         procedure :: get_g0w0_sigc_band_k => librpa_get_g0w0_sigc_band_k
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

   !> @brief Initialize runtime options to default values.
   !>
   !> Sets all options to their default settings. Must be called before
   !> modifying options and passing to computation functions.
   !>
   !> @param[in,out] opts Options structure to initialize.
   subroutine librpa_init_options(opts)
      implicit none
      class(LibrpaOptions), intent(inout) :: opts
      call error_on_call("librpa_init_options")
   end subroutine librpa_init_options

   !> @brief Set the output directory for LibRPA results.
   !>
   !> @param[in,out] opts       Options structure.
   !> @param[in]     output_dir Path to directory for output files.
   subroutine librpa_set_output_dir(opts, output_dir)
      implicit none
      class(LibrpaOptions), intent(inout) :: opts
      character(len=*), intent(in) :: output_dir
      call error_on_call("librpa_set_output_dir")
   end subroutine librpa_set_output_dir

   !> @brief Initialize the global computing environment of LibRPA
   !>
   !> It should be called after MPI initialization and before other LibRPA functions.
   !>
   !> @param  sw_redirect    Switch of redirecting standard output (default false)
   !> @param  redirect_path  Path of redirected output, only used when `sw_redirect` is true
   !> @param  sw_process     Switch of writing per-process output (default true)
   subroutine librpa_init_global(sw_redirect, redirect_path, sw_process)
      implicit none
      logical, intent(in), optional :: sw_redirect, sw_process
      character(len=*), intent(in), optional :: redirect_path
      call error_on_call("librpa_init_global")
   end subroutine librpa_init_global

   !> @brief Release all internal data and finalize the global computing environment of LibRPA
   !>
   !> It should be called after all LibRPA operations are finished.
   subroutine librpa_finalize_global()
      implicit none
      call error_on_call("librpa_finalize_global")
   end subroutine librpa_finalize_global

   !> @brief Get major version number.
   !> @return Major version (X in X.Y.Z).
   integer function librpa_get_major_version() result(v)
      implicit none
      v = -1
      call error_on_call("librpa_get_major_version")
   end function librpa_get_major_version

   !> @brief Get minor version number.
   !> @return Minor version (Y in X.Y.Z).
   integer function librpa_get_minor_version() result(v)
      implicit none
      v = -1
      call error_on_call("librpa_get_minor_version")
   end function librpa_get_minor_version

   !> @brief Get patch version number.
   !> @return Patch version (Z in X.Y.Z).
   integer function librpa_get_patch_version() result(v)
      implicit none
      v = -1
      call error_on_call("librpa_get_patch_version")
   end function librpa_get_patch_version

   !> @brief Run internal self-tests.
   !>
   !> Performs basic sanity checks on LibRPA functionality.
   !> Useful for debugging issues.
   subroutine librpa_test()
      implicit none
      call error_on_call("librpa_test")
   end subroutine librpa_test

   !> @brief Print profiling information.
   !>
   !> Outputs timing and memory usage statistics.
   subroutine librpa_print_profile()
      implicit none
      call error_on_call("librpa_print_profile")
   end subroutine librpa_print_profile

   !> @brief Create a new LibRPA handler instance.
   !>
   !> Allocates and initializes a new LibRPA handler associated with the given
   !> MPI communicator.
   !>
   !> @param[in,out] this  Handler to create.
   !> @param[in]     comm  MPI communicator (e.g., MPI_COMM_WORLD).
   subroutine librpa_create_handler(this, comm)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      integer, intent(in) :: comm
      call error_on_call("librpa_create_handler")
   end subroutine librpa_create_handler

   !> @brief Destroy a LibRPA handler instance.
   !>
   !> Frees all internal resources associated with the handler.
   !>
   !> @param[in,out] this  Handler to destroy.
   subroutine librpa_destroy_handler(this)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      call error_on_call("librpa_destroy_handler")
   end subroutine librpa_destroy_handler

   !> @brief Set SCF wavefunction dimension.
   !>
   !> @param[in,out] this     Handler.
   !> @param[in]     nspins   Number of spin channels.
   !> @param[in]     nkpts    Number of k-points.
   !> @param[in]     nstates  Number of electronic states.
   !> @param[in]     nbasis   Number of basis functions.
   subroutine librpa_set_scf_dimension(this, nspins, nkpts, nstates, nbasis)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      integer, intent(in) :: nspins, nkpts, nstates, nbasis
      call error_on_call("librpa_set_scf_dimension")
   end subroutine librpa_set_scf_dimension

   !> @brief Set occupation numbers, eigenvalues, and Fermi level.
   !>
   !> @param[in,out] this     Handler.
   !> @param[in]     nspins   Number of spin channels.
   !> @param[in]     nkpts    Number of k-points.
   !> @param[in]     nstates  Number of electronic states.
   !> @param[in]     wg       Occupation numbers (nstates x nkpts x nspins).
   !> @param[in]     ekb      Eigenvalues (nstates x nkpts x nspins).
   !> @param[in]     efermi   Fermi level.
   subroutine librpa_set_wg_ekb_efermi(this, nspins, nkpts, nstates, wg, ekb, efermi)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      integer, intent(in) :: nspins, nkpts, nstates
      real(dp), intent(in) :: wg(nstates, nkpts, nspins)
      real(dp), intent(in) :: ekb(nstates, nkpts, nspins)
      real(dp), intent(in) :: efermi
      call error_on_call("librpa_set_wg_ekb_efermi")
   end subroutine librpa_set_wg_ekb_efermi

   !> @brief Set the wave-function expansion coefficients
   !>
   !> @param ispin          spin index (starting from 1) of the wave function
   !> @param ik             (global) k-point index (starting from 1) of the wave function
   !> @param nstates_local  local dimenstion (number of states) of the parsed wave-function
   !> @param nbasis_local   local dimenstion (number of basis functions) of the parsed wave-function
   !> @param wfc_cplx       Complex-valued wave function to parse
   subroutine librpa_set_wfc(this, ispin, ik, nstates_local, nbasis_local, wfc_cplx)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      integer, intent(in) :: ispin, ik, nstates_local, nbasis_local
      complex(dp), intent(in), target :: wfc_cplx(nbasis_local, nstates_local)
      call error_on_call("librpa_set_wfc")
   end subroutine librpa_set_wfc

   !> @brief Set the wave-function atomic basis
   !>
   !> @param natoms   number of atoms
   !> @param nbs_wfc  number of wave-function basis on each atom
   subroutine librpa_set_ao_basis_wfc(this, natoms, nbs_wfc)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      integer, intent(in) :: natoms
      integer, intent(in) :: nbs_wfc(natoms)
      call error_on_call("librpa_set_ao_basis_wfc")
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
      call error_on_call("librpa_set_ao_basis_aux")
   end subroutine librpa_set_ao_basis_aux

   !> @brief Set the direct and reciprocal lattice vectors
   !>
   !> Each column is a lattice/reciprocal lattice vector.
   !>
   !> @param latt     lattice vectors (in Bohr)
   !> @param recplatt reciprocal lattice vectors (in Bohr^-1)
   !>
   subroutine librpa_set_latvec_and_G(this, latt, recplatt)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      real(dp), dimension(3, 3), intent(in) :: latt, recplatt
      call error_on_call("librpa_set_latvec_and_G")
   end subroutine librpa_set_latvec_and_G

   !> @brief Set types and coordinates of the atoms in the model
   !>
   !> @param natoms     number of atoms
   !> @param types      species type of each atom
   !> @param pos_cart   Cartesian coordinates of each atom
   !>
   subroutine librpa_set_atoms(this, natoms, types, posi_cart)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      integer, intent(in) :: natoms
      integer, dimension(natoms), intent(in) :: types
      real(dp), dimension(3, natoms), intent(in) :: posi_cart
      call error_on_call("librpa_set_atoms")
   end subroutine librpa_set_atoms

   !> @brief Set k-point grid vectors.
   !>
   !> @param[in,out] this  Handler.
   !> @param[in]     nk1    Number of k-points along direction 1.
   !> @param[in]     nk2    Number of k-points along direction 2.
   !> @param[in]     nk3    Number of k-points along direction 3.
   !> @param[in]     kvecs  K-point vectors (3 x nk1*nk2*nk3, Cartesian).
   !>
   subroutine librpa_set_kgrids_kvec(this, nk1, nk2, nk3, kvecs)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      integer, intent(in) :: nk1, nk2, nk3
      real(dp), intent(in) :: kvecs(3, nk1*nk2*nk3)
      call error_on_call("librpa_set_kgrids_kvec")
   end subroutine librpa_set_kgrids_kvec

   !> @brief Set the mapping from full k-point list to the irreducbile sector
   !>
   !> Example: four-k-point case where the first two and last points are in the irreducbile sector,
   !>          and the third point is mapped to the second, then map_ibzk should be (1, 2, 2, 4)
   !>
   !> @param nkpts     number of k-points in the full Brillouin zone
   !> @param map_ibzk  mapping to the k-point in the irreducible sector
   !>
   subroutine librpa_set_ibz_mapping(this, nkpts, map_ibzk)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      integer, intent(in) :: nkpts
      integer, dimension(nkpts), intent(in) :: map_ibzk
      call error_on_call("librpa_set_ibz_mapping")
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
      call error_on_call("librpa_set_lri_coeff")
   end subroutine librpa_set_lri_coeff

   !> @brief Set bare Coulomb matrix elements (atom-pair format).
   !>
   !> @param[in,out] this         Handler.
   !> @param[in]     ik           K-point index (1-based).
   !> @param[in]     i_atom       Atom I index (1-based).
   !> @param[in]     j_atom       Atom J index (1-based).
   !> @param[in]     naux_i       Number of aux functions for i.
   !> @param[in]     naux_j       Number of aux functions for j.
   !> @param[in]     vq           Coulomb matrix (naux_i x naux_j, complex).
   !> @param[in]     vq_threshold  Threshold for screening.
   subroutine librpa_set_aux_bare_coulomb_k_atom_pair &
         (this, ik, i_atom, j_atom, naux_i, naux_j, vq, vq_threshold)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      integer, intent(in) :: ik, i_atom, j_atom, naux_i, naux_j
      complex(dp), intent(in) :: vq(naux_i, naux_j)
      real(dp), intent(in) :: vq_threshold
      call error_on_call("librpa_set_aux_bare_coulomb_k_atom_pair")
   end subroutine librpa_set_aux_bare_coulomb_k_atom_pair

   !> @brief Set truncated Coulomb matrix elements (atom-pair format).
   !>
   !> @param[in,out] this         Handler.
   !> @param[in]     ik           K-point index (1-based).
   !> @param[in]     i_atom       Atom I index (1-based).
   !> @param[in]     j_atom       Atom J index (1-based).
   !> @param[in]     naux_i       Number of aux functions for i.
   !> @param[in]     naux_j       Number of aux functions for j.
   !> @param[in]     vq           Coulomb matrix (naux_i x naux_j, complex).
   !> @param[in]     vq_threshold  Threshold for screening.
   subroutine librpa_set_aux_cut_coulomb_k_atom_pair &
         (this, ik, i_atom, j_atom, naux_i, naux_j, vq, vq_threshold)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      integer, intent(in) :: ik, i_atom, j_atom, naux_i, naux_j
      complex(dp), intent(in) :: vq(naux_i, naux_j)
      real(dp), intent(in) :: vq_threshold
      call error_on_call("librpa_set_aux_cut_coulomb_k_atom_pair")
   end subroutine librpa_set_aux_cut_coulomb_k_atom_pair

   !> @brief Set bare Coulomb matrix elements (2D block format).
   !>
   !> @param[in,out] this        Handler.
   !> @param[in]     ik          K-point index (1-based).
   !> @param[in]     mu_begin    Starting mu index (1-based).
   !> @param[in]     mu_end      Ending mu index (inclusive).
   !> @param[in]     nu_begin    Starting nu index (1-based).
   !> @param[in]     nu_end      Ending nu index (inclusive).
   !> @param[in]     vq          Coulomb matrix (complex).
   subroutine librpa_set_aux_bare_coulomb_k_2d_block &
         (this, ik, mu_begin, mu_end, nu_begin, nu_end, vq)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      integer, intent(in) :: ik, mu_begin, mu_end, nu_begin, nu_end
      complex(dp), intent(in) :: vq(mu_end-mu_begin+1, nu_end-nu_begin+1)
      call error_on_call("librpa_set_aux_bare_coulomb_k_2d_block")
   end subroutine librpa_set_aux_bare_coulomb_k_2d_block

   !> @brief Set truncated Coulomb matrix elements (2D block format).
   !>
   !> @param[in,out] this        Handler.
   !> @param[in]     ik          K-point index (1-based).
   !> @param[in]     mu_begin    Starting mu index (1-based).
   !> @param[in]     mu_end      Ending mu index (inclusive).
   !> @param[in]     nu_begin    Starting nu index (1-based).
   !> @param[in]     nu_end      Ending nu index (inclusive).
   !> @param[in]     vq          Coulomb matrix (complex).
   subroutine librpa_set_aux_cut_coulomb_k_2d_block &
         (this, ik, mu_begin, mu_end, nu_begin, nu_end, vq)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      integer, intent(in) :: ik, mu_begin, mu_end, nu_begin, nu_end
      complex(dp), intent(in) :: vq(mu_end-mu_begin+1, nu_end-nu_begin+1)
      call error_on_call("librpa_set_aux_cut_coulomb_k_2d_block")
   end subroutine librpa_set_aux_cut_coulomb_k_2d_block

   !> @brief Set dielectric function on imaginary frequency axis.
   !> @param[in,out] this            Handler.
   !> @param[in]     nfreq           Number of frequency points.
   !> @param[in]     omegas_imag     Imaginary frequency values.
   !> @param[in]     dielect_func    Dielectric function values.
   subroutine librpa_set_dielect_func_imagfreq(this, nfreq, omegas_imag, dielect_func)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      integer, intent(in) :: nfreq
      real(dp), dimension(nfreq), intent(in) :: omegas_imag
      real(dp), dimension(nfreq), intent(in) :: dielect_func
      call error_on_call("librpa_set_dielect_func_imagfreq")
   end subroutine librpa_set_dielect_func_imagfreq

   !> @brief Set k-points for band structure calculations.
   !>
   !> @param[in,out] this        Handler.
   !> @param[in]     nkpts_band  Number of band k-points.
   !> @param[in]     kfrac_band  Band k-point coordinates (3 x nkpts_band, fractional).
   subroutine librpa_set_band_kvec(this, nkpts_band, kfrac_band)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      integer, intent(in) :: nkpts_band
      real(dp), intent(in) :: kfrac_band(3, nkpts_band)
      call error_on_call("librpa_set_band_kvec")
   end subroutine librpa_set_band_kvec

   !> @brief Set occupation numbers and eigenvalues for band k-points.
   !>
   !> @param[in,out] this         Handler.
   !> @param[in]     nspins       Number of spin channels.
   !> @param[in]     nkpts_band   Number of band k-points.
   !> @param[in]     nstates      Number of states.
   !> @param[in]     occ          Occupation numbers (nstates x nkpts_band x nspins).
   !> @param[in]     eig          Eigenvalues (nstates x nkpts_band x nspins).
   subroutine librpa_set_band_occ_eigval(this, nspins, nkpts_band, nstates, occ, eig)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      integer, intent(in) :: nspins, nkpts_band, nstates
      real(dp), intent(in) :: occ(nstates, nkpts_band, nspins)
      real(dp), intent(in) :: eig(nstates, nkpts_band, nspins)
      call error_on_call("librpa_set_band_occ_eigval")
   end subroutine librpa_set_band_occ_eigval

   !> @brief Set the wave-function expansion coefficients for band calculation
   !>
   !> @param ispin          spin index (starting from 1) of the wave function
   !> @param ik_band        (global) k-point index (starting from 1) of the wave function
   !> @param nstates_local  local dimenstion (number of states) of the parsed wave-function
   !> @param nbasis_local   local dimenstion (number of basis functions) of the parsed wave-function
   !> @param wfc_cplx       Complex-valued wave function to parse
   subroutine librpa_set_wfc_band(this, ispin, ik_band, nstates_local, nbasis_local, wfc_cplx)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      integer, intent(in) :: ispin, ik_band, nstates_local, nbasis_local
      complex(dp), intent(in), target :: wfc_cplx(nbasis_local, nstates_local)
      call error_on_call("librpa_set_wfc_band")
   end subroutine librpa_set_wfc_band

   !> @brief Reset band structure data.
   !> @param[in,out] this  Handler.
   subroutine librpa_reset_band_data(this)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      call error_on_call("librpa_reset_band_data")
   end subroutine librpa_reset_band_data

   !> @brief Construct and return frequency grids.
   !>
   !> @param[in,out] this    Handler.
   !> @param[in,out] opts    Runtime options.
   !> @param[out]    omegas  Frequency values.
   !>
   !> @param[out]    weights Quadrature weights.
   subroutine librpa_get_imaginary_frequency_grids(this, opts, omegas, weights)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      type(LibrpaOptions), intent(inout) :: opts
      real(dp), allocatable, intent(inout) :: omegas(:), weights(:)
      call error_on_call("librpa_get_imaginary_frequency_grids")
   end subroutine librpa_get_imaginary_frequency_grids

   !> @brief Compute RPA correlation energy.
   !>
   !> @param[in,out] this          Handler.
   !> @param[in,out] opts          Runtime options.
   !> @param[in]     nkpts_ibz     Number of irreducible k-points.
   !> @param[out]    contrib_ibzk  Complex correlation contribution per k-point.
   !>
   !> @return Total RPA correlation energy.
   real(dp) function librpa_get_rpa_correlation_energy(this, opts, nkpts_ibz, contrib_ibzk) result(e)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      type(LibrpaOptions), intent(inout) :: opts
      integer, intent(in) :: nkpts_ibz
      complex(dp), dimension(nkpts_ibz), intent(inout) :: contrib_ibzk
      e = 0.0d0
      call error_on_call("librpa_get_rpa_correlation_energy")
   end function librpa_get_rpa_correlation_energy

   !> @brief Build exact-exchange matrix in real space.
   !> @param[in,out] this  Handler.
   !> @param[in,out] opts  Runtime options.
   subroutine librpa_build_exx(this, opts)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      type(LibrpaOptions), intent(inout) :: opts
      call error_on_call("librpa_build_exx")
   end subroutine librpa_build_exx

   !> @brief Get exact-exchange potential for k-grid states.
   !> @param[in,out] this          Handler.
   !> @param[in,out] opts          Runtime options.
   !> @param[in]     n_spins       Number of spin channels.
   !> @param[in]     n_kpts_this   Number of k-points on this process.
   !> @param[in]     iks_this      List of k-point indices (1-based).
   !> @param[in]     i_state_low   First state index (1-based, inclusive).
   !> @param[in]     i_state_high  Last state index (1-based, inclusive).
   !> @param[out]    vexx          Exact-exchange potentials.
   subroutine librpa_get_exx_pot_kgrid(this, opts, n_spins, n_kpts_this, iks_this, &
                                       i_state_low, i_state_high, vexx)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      type(LibrpaOptions), intent(inout) :: opts
      integer, contiguous, dimension(:), intent(in) :: iks_this
      integer, intent(in) :: n_spins, n_kpts_this, i_state_low, i_state_high
      real(dp), dimension(i_state_high - i_state_low + 1, n_kpts_this, n_spins), intent(inout) :: vexx
      call error_on_call("librpa_get_exx_pot_kgrid")
   end subroutine librpa_get_exx_pot_kgrid

   !> @brief Get exact-exchange potential for band k-points.
   !> @param[in,out] this              Handler.
   !> @param[in,out] opts              Runtime options.
   !> @param[in]     n_spins           Number of spin channels.
   !> @param[in]     n_kpts_band_this  Number of band k-points on this process.
   !> @param[in]     iks_band_this     List of band k-point indices (1-based).
   !> @param[in]     i_state_low       First state index (1-based, inclusive).
   !> @param[in]     i_state_high      Last state index (1-based, inclusive).
   !> @param[out]    vexx_band         Exact-exchange potentials for band k-points.
   subroutine librpa_get_exx_pot_band_k(this, opts, n_spins, n_kpts_band_this, iks_band_this, &
                                        i_state_low, i_state_high, vexx_band)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      type(LibrpaOptions), intent(inout) :: opts
      integer, contiguous, dimension(:), intent(in) :: iks_band_this
      integer, intent(in) :: n_spins, n_kpts_band_this, i_state_low, i_state_high
      real(dp), dimension(i_state_high - i_state_low + 1, n_kpts_band_this, n_spins), intent(inout) :: vexx_band
      call error_on_call("librpa_get_exx_pot_band_k")
   end subroutine librpa_get_exx_pot_band_k

   !> @brief Build G0W0 self-energy matrix in real space.
   !> @param[in,out] this  Handler.
   !> @param[in,out] opts  Runtime options.
   subroutine librpa_build_g0w0_sigma(this, opts)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      type(LibrpaOptions), intent(inout) :: opts
      call error_on_call("librpa_build_g0w0_sigma")
   end subroutine librpa_build_g0w0_sigma

   !> @brief Get G0W0 correlation self-energy for k-grid states.
   !> @param[in,out] this           Handler.
   !> @param[in,out] opts           Runtime options.
   !> @param[in]     n_spins        Number of spin channels.
   !> @param[in]     n_kpts_this    Number of k-points on this process.
   !> @param[in]     iks_this       List of k-point indices (1-based).
   !> @param[in]     i_state_low    First state index (1-based, inclusive).
   !> @param[in]     i_state_high   Last state index (1-based, inclusive).
   !> @param[in]     vxc            XC potential for selected states.
   !> @param[in]     vexx           Exact-exchange potential for selected states.
   !> @param[out]    sigc           Correlation self-energy (complex).
   subroutine librpa_get_g0w0_sigc_kgrid(this, opts, n_spins, n_kpts_this, iks_this, &
                                         i_state_low, i_state_high, vxc, vexx, sigc)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      type(LibrpaOptions), intent(inout) :: opts
      integer, contiguous, dimension(:), intent(in) :: iks_this
      integer, intent(in) :: n_spins, n_kpts_this, i_state_low, i_state_high
      real(dp), dimension(i_state_high - i_state_low + 1, n_kpts_this, n_spins), intent(in) :: vxc
      real(dp), dimension(i_state_high - i_state_low + 1, n_kpts_this, n_spins), intent(in) :: vexx
      complex(dp), dimension(i_state_high - i_state_low + 1, n_kpts_this, n_spins), intent(inout) :: sigc
      call error_on_call("librpa_get_g0w0_sigc_kgrid")
   end subroutine librpa_get_g0w0_sigc_kgrid

   !> @brief Get G0W0 correlation self-energy for band k-points.
   !> @param[in,out] this              Handler.
   !> @param[in,out] opts              Runtime options.
   !> @param[in]     n_spins           Number of spin channels.
   !> @param[in]     n_kpts_band_this  Number of band k-points on this process.
   !> @param[in]     iks_band_this     List of band k-point indices (1-based).
   !> @param[in]     i_state_low       First state index (1-based, inclusive).
   !> @param[in]     i_state_high      Last state index (1-based, inclusive).
   !> @param[in]     vxc_band          XC potential for band states.
   !> @param[in]     vexx_band         Exact-exchange potential for band states.
   !> @param[out]    sigc_band        Correlation self-energy for band (complex).
   subroutine librpa_get_g0w0_sigc_band_k(this, opts, n_spins, n_kpts_band_this, iks_band_this, &
                                          i_state_low, i_state_high, vxc_band, vexx_band, sigc_band)
      implicit none
      class(LibrpaHandler), intent(inout) :: this
      type(LibrpaOptions), intent(inout) :: opts
      integer, contiguous, dimension(:), intent(in) :: iks_band_this
      integer, intent(in) :: n_spins, n_kpts_band_this, i_state_low, i_state_high
      real(dp), dimension(i_state_high - i_state_low + 1, n_kpts_band_this, n_spins), intent(in) :: vxc_band
      real(dp), dimension(i_state_high - i_state_low + 1, n_kpts_band_this, n_spins), intent(in) :: vexx_band
      complex(dp), dimension(i_state_high - i_state_low + 1, n_kpts_band_this, n_spins), intent(inout) :: sigc_band
      call error_on_call("librpa_get_g0w0_sigc_band_k")
   end subroutine librpa_get_g0w0_sigc_band_k

end module librpa_f03
