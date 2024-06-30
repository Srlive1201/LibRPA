module librpa

   use iso_c_binding
   implicit none

   public

   type, bind(c) :: LibRPAParams
      character(kind=c_char, len=1) :: task(20)
      ! character(kind=c_char, len=1)  :: parallel_routing(20)
      ! character(kind=c_char, len=1)  :: tfgrids_type(20)
      character(kind=c_char, len=1) :: output_file(100)
      character(kind=c_char, len=1) :: output_dir(100)

      integer(c_int) :: nfreq
      ! integer(c_int) :: debug
      ! integer(c_int) :: use_scalapack_ecrpa

      ! real(c_double) :: gf_R_threshold
      ! real(c_double) :: cs_threshold
      ! real(c_double) :: vq_threshold
      ! real(c_double) :: sqrt_coulomb_threshold
      ! real(c_double) :: libri_chi0_threshold_C
      ! real(c_double) :: libri_chi0_threshold_G
      ! real(c_double) :: libri_exx_threshold_CSM
      ! real(c_double) :: libri_exx_threshold_C
      ! real(c_double) :: libri_exx_threshold_D
      ! real(c_double) :: libri_exx_threshold_V
   end type LibRPAParams

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
      subroutine set_librpa_params(params) bind(c, name="set_librpa_params")
         use, intrinsic :: iso_c_binding
         import :: LibRPAParams
         type(LibRPAParams), intent(in) :: params
      end subroutine
   end interface

   interface
      subroutine get_default_librpa_params(params) bind(c, name="get_default_librpa_params")
         use, intrinsic :: iso_c_binding
         import :: LibRPAParams
         type(LibRPAParams), intent(inout) :: params
      end subroutine
   end interface

contains

end module
