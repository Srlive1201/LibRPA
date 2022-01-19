!****s* FHI-aims/post_scf_correlation_treatment_p0
!  NAME
!   post_scf_correlation_treatment_p0
!  SYNOPSIS

      subroutine post_scf_correlation_treatment_p0()

!  PURPOSE
!  This subroutine evaluate the necessary matrix elements (3-center overlap,
!  coulomb matrix) used for accurate correlation energy calculations beyond DFT
!  (e.g. MP2, RPA, RPA+, etc).
!
!  USES

      use dimensions
      use species_data
      use runtime_choices
      use pbc_lists, only : k_point_list
      use prodbas, only : OVLP_TYPE_COULOMB
      use hartree_fock
      use hartree_fock_p0
      use physics
      use gw_para
      use my_lvl_triples, only : my_cleanup_lvl_triples, my_initialize_lvl_triples
      use lvl_triples, only : cleanup_lvl_triples, initialize_lvl_triples
      use tight_binding_auxmat
      use timing
      use mpi_tasks
      use sbt_overlap_aims
      use calculate_fock_matrix_p0, only : cleanup_fock_matrix_calculations, &
          init_fock_matrix_calculations, evaluate_exchange_matr_realspace_p0
      use crpa_blacs
      use lvl_tricoeff
      use KS_optical_properties, only : KS_dielectric_calculation_imagfreq
      use pbc_lists, only: k_phase, n_cells, inv_irk_point_mapping
      use scalapack_wrapper,     only : mxld, mxcol, eigenvec, eigenvec_complex
      implicit none

!  ARGUMENTS

!  INPUT
!    none
!  OUTPUT
!    none
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

!  local variables
      character*20 filename
      logical :: need_o3fn
      real*8 :: time_prodbas_add, clock_time_prodbas_add
      real*8 :: rpa_tot_en
      real*8 :: pt2_tot_en
      real*8 :: en_se, en_rse
      real*8 :: exx_ene, d_exx_ene(3) !SVL dummy variables, not actually used anywhere but needed for hf_realspace call 

      integer, allocatable :: n_homo_k(:,:)
      integer :: n_low_state_polar
      integer :: n_k_points_min
      integer :: bb_row_s, bb_col_s

!   GW self-energy
      complex*16, allocatable :: gw_selfenergy(:,:,:,:)
      complex*16, allocatable :: gw_selfe_band(:,:,:,:)
!   exact-exhcange term
      real*8, allocatable :: exact_x_kspace(:,:,:)
      real*8, allocatable :: xc_kspace(:,:,:)
      real*8, allocatable :: qp_energy(:,:,:)
!      complex*16, allocatable :: coulomb_gamma(:,:)
      complex*16, allocatable :: coulomb_gamma_blacs(:,:)

!   the XC part of the DFT hamiltonian in realspace
      real*8, allocatable :: xc_realspace(:,:)
      real*8, allocatable :: x_realspace(:,:)
      real*8, allocatable :: c_realspace(:,:)
! Counter
      integer i_spin
      integer i_state
      integer i_k_point, i_k_point_local
      integer i_band, i_cell
      integer i_k_point_index
      integer i_prodbas_1, i_prodbas_2
      integer i_bb_1, i_bb_2, i_bb_s_1, i_bb_s_2
      integer i_irk_point,i_irk_point_local

      integer :: info
      character*150 :: info_str
      character(*), parameter :: &
          func='post_scf_correlation_treatment_p0'
      ! for parallel threads
      integer :: omp_get_max_threads, omp_get_thread_num, nrThreads, myThreadId, i_thread ,mode
      integer :: ks_i,ks_j,ks_k,ks_l,i_k
      real*8, allocatable :: lvl_tricoeff_bravais(:,:,:,:)
      real*8 :: k_vec(3)
      character(len=100) eFile
      real*8 :: nt_d_nk
      nrThreads   = 1
      myThreadId = 1
!$    nrThreads = omp_get_max_threads()
!$    myThreadId = omp_get_thread_num() + 1

      call perfinit('evpost')
      call perfon('pscf')
      call get_timestamps(time_prodbas_add, clock_time_prodbas_add)
!      call perfon('ini1')
!  n_irk_points, irk_point_list, irk_point_mapping, inv_irk_point_mapping,
!  irk_point_included (all defined in pbc_lists.f90) will be determined here
      call determine_irreducible_k_grid ()

      if(real_eigenvectors) then
        allocate(KS_eigenvector_irk(n_basis,n_states,n_spin,n_irk_points_task),stat=info) 
        call check_allocation(info, 'KS_eigenvector_irk', func)
        allocate(KS_eigenvector_complex_irk(1,1,1,1),stat=info) 
        call check_allocation(info, 'KS_eigenvector_complex_irk', func)
      else
        allocate(KS_eigenvector_irk(1,1,1,1),stat=info) 
        call check_allocation(info, 'KS_eigenvector_irk', func)
        allocate(KS_eigenvector_complex_irk(n_basis,n_states,n_spin,n_irk_points_task),stat=info) 
        call check_allocation(info, 'KS_eigenvector_complex_irk', func)
      endif


      call distribute_irreducible_eigenvectors &
           ( KS_eigenvector, KS_eigenvector_complex, &
             KS_eigenvector_irk, KS_eigenvector_complex_irk )

      ! --- Ensure ovlp_3fn / coeff_3fn_ten&coulomb_matr_lvl

      ! for non-hartree-fock self-consistent calculations, we need to
      ! construct the product (auxiliary) basis functions, and evaluate
      ! the 3-center overlap and coulomb matrix elements here.
!call perfoff

!       if(use_gw) then
!           if(.not. allocated(dielec_func_imagfreq)) then
!              allocate(dielec_func_imagfreq(n_full_freq),stat=info)
!              call check_allocation(info, 'dielec_func_imagfreq            ')
!           endif
!           if(use_scalapack) then
!              call KS_dielectric_calculation_imagfreq( mxld, mxcol, n_spin, n_states, KS_eigenvalue, eigenvec, &
!                     eigenvec_complex, occ_numbers, chemical_potential, partition_tab, l_shell_max )
!
!           else
!              call KS_dielectric_calculation_imagfreq( n_basis, n_states, n_spin, n_states, KS_eigenvalue, KS_eigenvector, &
!                     KS_eigenvector_complex, occ_numbers, chemical_potential, partition_tab, l_shell_max )
!           endif
!           stop
!       endif

      if (.not.use_hf_kspace) then
!call perfon('ini2')
        if (.not.allocated(n_homo)) then
         allocate (n_homo(n_spin))
          n_homo(:)=0
        endif

        call initialize_prodbas()


!call perfoff
!call perfon('ini3')
        if(.not.use_hf_realspace)then
           call allocate_hartree_fock_p0()
           if(.not.use_hf_kspace_with_rpa)then
              call cleanup_fock_matrix_calculations
              call init_fock_matrix_calculations(.false., .false.)
              call evaluate_exchange_matr_realspace_p0(KS_eigenvector,KS_eigenvector_complex,occ_numbers, &
                   exx_ene, d_exx_ene, .false., .false.)
              call cleanup_fock_matrix_calculations
           endif
        endif
!call perfoff
call perfon('ini4')

        if(myid.eq.0) then
           write(use_unit,*)
           write(use_unit,'(2X,A)')"------------------------------------------------------------------"
           write(use_unit,*)
           write(use_unit,'(2X,A)') "Post-SCF correlation calculation starts ..."
           write(use_unit,*)
        endif
        call init_crpa_blacs_distribution(n_full_freq)
if (myid.eq.0) then
   print*,'reso: bb,b,st',n_basbas,n_basis,n_states
   print*,'reso: mnbb,nc,nf',max_n_basis_sp,n_cells,n_full_freq
end if
        if (myid_col.eq.0) then
           allocate(lvl_tricoeff_mod_r(lbb_row:ubb_row,max_n_basis_sp,n_states,n_spin,n_ks_points_task),stat=info)
        else
           allocate(lvl_tricoeff_mod_r(lbb_row:lbb_row,1,1,1,1),stat=info)
        end if
        call check_allocation(info, 'lvl_tricoeff_mod_r', func)


        allocate(coulomb_matr_blacs(lbb_row:ubb_row,lbb_col:ubb_col, n_irkq_points_task),stat=info)
        call check_allocation(info, 'coulomb_matr_blacs                 ')

        allocate(coulomb_cut_blacs(lbb_row:ubb_row,lbb_col:ubb_col, n_irkq_points_task),stat=info)
        call check_allocation(info, 'coulomb_matr_blacs                 ')


        allocate(coulomb_gamma_blacs(lbb_row:ubb_row,lbb_col:ubb_col),stat=info)
        call check_allocation(info, 'coulomb_gamma_blacs                 ')

        n_k_points_min=max(n_k_points_task,1)

        allocate(kq_point_list(n_k_points,n_k_points),stat=info) 
        call check_allocation(info, 'kq_point_list', func)
        allocate(kpq_point_list(n_k_points,n_k_points),stat=info) 
        call check_allocation(info, 'kpq_point_list', func)

        if (use_threadsafe_gwinit) then
          call my_initialize_lvl_triples(OVLP_TYPE_COULOMB)
        else
          call initialize_lvl_triples(OVLP_TYPE_COULOMB)
        endif
        if(myid.eq.0) then
           write(use_unit,*)
           write(use_unit,'(2X,A)') "Finished with initialization of triples (Coulomb norm) ..."
        endif

        call initialize_tb_auxmat(1, OVLP_TYPE_COULOMB)
        if(myid.eq.0) then
           write(use_unit,*)
           write(use_unit,'(2X,A)') "Finished with initialization of Coulomb interaction matrix ..."
           write(use_unit,*)
        endif
!shirong add 2021-01-19       
        open(1, file='band_out')
        write(1,*)n_k_points
        write(1,*)n_spin
        write(1,*)n_states
        write(1,*)n_basis
        write(1,*)chemical_potential
        do i_k_point = 1,n_k_points,1
          do i_spin = 1, n_spin, 1
            write(1,*)i_k_point,i_spin
            do i_state = 1, n_states, 1
                write(1,'(2X,I5,6X,F8.5,5X,F17.9,4X,F15.5)') &
                     i_state, occ_numbers(i_state, i_spin, i_k_point), &
                     (KS_eigenvalue(i_state, i_spin, i_k_point)), &
                     (KS_eigenvalue(i_state, i_spin, i_k_point)*hartree)
             end do
          end do
        end do
        close(1)
        ! print *, ' INFO in post_scf ',n_basis,n_states, n_spin,n_k_points_task,n_k_points
        ! print *, ' KS_eigenvector  size ',size(KS_eigenvector)
        ! print *, ' KS_eigenvector  dim1 ',size(KS_eigenvector,dim=1)
        ! print *, ' KS_eigenvector  dim2 ',size(KS_eigenvector,dim=2)
        ! print *, ' KS_eigenvector  dim3 ',size(KS_eigenvector,dim=3)
        ! print *, ' KS_eigenvector  dim4 ',size(KS_eigenvector,dim=4)
        i_k=0
        
        write(eFile,*)myid
        eFile='KS_eigenvector_'//trim(adjustl(eFile))//'.txt'
        
       ! write(2,*)'KS_eigenvector dim',real_eigenvectors,n_basis,n_states, n_spin,n_k_points_task
        nt_d_nk=n_tasks*1.d0/n_k_points
        open(2,file=eFile)
        if(use_scalapack) then
        do i_k_point=1,n_k_points
          if(myid == CEILING((i_k_point-1)*nt_d_nk)) then
            write(2,*)i_k_point
            if(real_eigenvectors) then
              do ks_i=1, n_basis,1
                do ks_j=1, n_states,1
                  do ks_k=1,n_spin,1
                    write(2,'(2F24.15)')KS_eigenvector(ks_i,ks_j,ks_k,1),0.d0
                  end do
                end do
              end do
            else
              do ks_i=1, n_basis,1
                do ks_j=1, n_states,1
                  do ks_k=1,n_spin,1
                    write(2,'(2F24.15)')KS_eigenvector_complex(ks_i,ks_j,ks_k,1)
                  end do
                end do
              end do
            end if
          end if
        end do
      else
        do i_k_point = 1, n_k_points
          if(myid == MOD(i_k_point, n_tasks)  .and. myid <= n_k_points) then
            i_k=i_k+1
            write(2,*)i_k_point
            if(real_eigenvectors) then
              do ks_i=1, n_basis,1
                do ks_j=1, n_states,1
                  do ks_k=1,n_spin,1
                    write(2,'(2F24.15)')KS_eigenvector(ks_i,ks_j,ks_k,i_k),0.d0
                  end do
                end do
              end do
            else
              do ks_i=1, n_basis,1
                do ks_j=1, n_states,1
                  do ks_k=1,n_spin,1
                    write(2,'(2F24.15)')KS_eigenvector_complex(ks_i,ks_j,ks_k,i_k)
                  end do
                end do
              end do
            end if
          end if
        end do 
      end if
      close(2)

        ! open(3,file='eigenvalue')
        ! write(3,*)'KS_eigenvalue dim', n_states,n_spin,n_k_points,real_eigenvectors
        ! do ks_i=1,n_states,1
        !   do ks_j=1,n_spin,1
        !     do ks_k=1,n_band_kpoints,1
        !       write(3,'(F15.10)')KS_eigenvalue(ks_i,ks_j,ks_k)
        !     end do
        !   end do
        ! end do
        ! close(3)

        open(4,file='stru_out')
        write(4,*)lattice_vector(:,1)
        write(4,*)lattice_vector(:,2)
        write(4,*)lattice_vector(:,3)
        write(4,*)recip_lattice_vector(:,1)
        write(4,*)recip_lattice_vector(:,2)
        write(4,*)recip_lattice_vector(:,3)
        write(4,*)n_k_points_xyz(:)
        do i_k_point = 1, n_k_points, 1
          k_vec = matmul(recip_lattice_vector,k_point_list(i_k_point,:))
          write(4,*)k_vec
        end do


        
       ! print*, 'post_scf n_cells  ',n_cells,'      n_task  ',n_tasks
        if(mod(n_cells, n_tasks).eq.0) then
          n_cells_task = n_cells/n_tasks
        else
          n_cells_task = n_cells/n_tasks+1
        endif 
       ! print *, 'post_scf  n_cell_tasks   ', n_cells_task
        ! if(.not. allocated(lvl_tricoeff_bravais)) then 
        !   allocate(lvl_tricoeff_bravais(n_basis,n_basis,n_basbas,n_cells_task),stat=info)
        !   call check_allocation(info,'lvl_tricoeff_bravais',func)
        ! endif
        call get_lvl_tricoeff_bravais(n_cells_task)

        ! if(allocated(lvl_tricoeff_bravais)) then 
        !   deallocate(lvl_tricoeff_bravais)
        ! endif
!        if(n_k_points.gt.1) then
           !ACHTUNG: this section can be removed by removing lvl_tricoef_recip* from evaluate_exchange_matr_kspace_p0!!
!           allocate(lvl_tricoeff_recip1(max_n_basbas_sp,n_basis,n_basis,n_k_points_min),stat=info) 
!           call check_allocation(info, 'lvl_tricoeff_recip1', func)
!           allocate(lvl_tricoeff_recip2(max_n_basbas_sp,n_basis,n_basis),stat=info) 
!           call check_allocation(info, 'lvl_tricoeff_recip2', func)
!           call get_lvl_tricoeff_recip(n_cells_task,lvl_tricoeff_recip1,lvl_tricoeff_recip2)
!        end if


        ! compute the number of real-space unit cells locally. parallelization in comm_col
        ! the entry n_cells+1 contains the contribution from the second atom
        n_cells_task = 0
        do i_cell = 1, n_cells+1
           if(myid_col.eq.mod(i_cell,n_tasks_col) .and.(myid_col .le. n_cells+1)) then
              n_cells_task = n_cells_task + 1
           endif
        enddo
        
        !compute the lvl_tricoeff_cell array only once
        call gw_init_lvl_tricoeff_recip(n_cells,n_cells_task,n_k_points, n_k_points_task, n_ks_points_task,&
             k_phase)
        if(myid.eq.0) then
           write(use_unit,*)
           write(use_unit,'(2X,A)') "Finished with initialization of the periodic LVL triple coefficients ..."
           write(use_unit,*)
        endif

        if (use_threadsafe_gwinit) then
           call my_cleanup_lvl_triples()
        else
           call cleanup_lvl_triples()
        endif

        call gw_get_lvl_tricoeff_recip(n_cells, n_k_points, n_k_points_task, n_ks_points_task,&
             k_phase, KS_eigenvector, KS_eigenvector_complex, lvl_tricoeff_mod_r)
        if(myid.eq.0) then
           write(use_unit,*)
           write(use_unit,'(2X,A)') "Finishing computation of the periodic LVL triple coefficients ..."
           write(use_unit,*)
        endif
        
        !for GW, the lvl_tricoeff_cell array is still needed
        if(.not.use_periodic_gw) call gw_cleanup_lvl_tricoeff_recip

        call deallocate_tb_auxmat()
        if(use_hse .and. hse_omega_hf /= 0.d0 .and. (.not. use_gw_and_hse) &
            .and. (.not. use_dftpt2_and_hse)) then
           call initialize_tb_auxmat(1, OVLP_TYPE_HSE)
        else
           call initialize_tb_auxmat(1, OVLP_TYPE_CUT)
           !call initialize_tb_auxmat(1, OVLP_TYPE_CUT_ANALYTIC)
!          call initialize_periodic_tb_auxmat(1, 1.d0)
        endif

        call determine_k_minus_q_list(kq_point_list,kpq_point_list)

!        call get_coulomb_matr_recip(coulomb_matr_recip,1)
!call perfon('gcmb')
        call get_coulomb_matr_blacs(coulomb_matr_blacs,1)
        if(use_periodic_gw) then
          coulomb_cut_blacs = coulomb_matr_blacs
        endif
!call perfoff
        if((.not.use_hf_realspace).and.use_hf_kspace_with_rpa)then

            print*,'ACHTUNG: noch nicht an coulomb_matr_blacs angepasst'
            stop
           call evaluate_exchange_matr_kspace_p0 &
                (KS_eigenvalue,KS_eigenvector,KS_eigenvector_complex,occ_numbers)
        endif
!<ADDED CODE>
        call deallocate_tb_auxmat()
!</ADDED CODE>


       if((.not. gamma_cut_coulomb) .or. use_periodic_gw) then
         do i_irk_point_local = 1, n_irkq_points_task
           i_irk_point=n_tasks_irkq*(i_irk_point_local-1) + myid_irkq + 1
           i_k_point = inv_irk_point_mapping(i_irk_point)
           if(all(abs(k_point_list(i_k_point,:)).lt.1.e-10)) then
              coulomb_gamma_blacs(:,:)= coulomb_matr_blacs(:,:,i_irk_point_local)        
           endif
         enddo

         call deallocate_tb_auxmat()
         call initialize_periodic_tb_auxmat(1, 1.d0)

         call get_coulomb_matr_blacs(coulomb_matr_blacs,1)

!         if(use_gw .and. use_gw_gamma_corr) then
!            call get_coulomb_coeff_blacs(coulomb_cut_blacs,1)
!         endif

! test, full Coulomb operator
         if(.not.use_periodic_gw) then
           do i_irk_point_local = 1, n_irkq_points_task

             i_irk_point=n_tasks_irkq*(i_irk_point_local-1) + myid_irkq + 1
             i_k_point = inv_irk_point_mapping(i_irk_point)
             if(all(abs(k_point_list(i_k_point,:)).lt.1.e-10)) then
                coulomb_matr_blacs(:,:,i_irk_point_local) = coulomb_gamma_blacs(:,:)        
             endif
! end of do i_irk_point_local
           enddo
         endif

! end of if(.not. gamma_cut_coulomb) 
       endif
       call perfoff
! end of if(.not.use_hf_realspace)
   endif


      call deallocate_tb_auxmat()
!      call initialize_tb_auxmat(1, OVLP_TYPE_COULOMB)

      call allocate_gw()

!     print * ,'freq_grid_type:' ,freq_grid_type
      if(freq_grid_type.eq.0) then
         call tf_ini(n_freq,n_full_freq, omegamax,omegamax, &
              omega_grid,omega_full_grid, &
              womega,womega_full,.true.)
      elseif(freq_grid_type.eq.1) then
         n_freq=n_full_freq
         call tf_ini_trans(n_freq,n_full_freq, omegamax,omegamax, &
              omega_grid,omega_full_grid, &
              womega,womega_full,.true.)
      endif
      call get_timestamps(time_rse_corr, clock_time_rse_corr)

      if(use_rpa_ene) then
call perfon('pscrp')
!        if (use_threadsafe_gwinit) then
!          if(.not.use_hf_kspace) call my_cleanup_lvl_triples()
!        else
!          if(.not.use_hf_kspace) call cleanup_lvl_triples()
!        endif

        call evaluate_single_excitation_correction_p0 &
           ( n_high_state, &
             occ_numbers, &
             hartree_potential,rho,rho_gradient, &
             kinetic_density, &
             partition_tab, l_shell_max, &
             hamiltonian, &
             KS_eigenvalue,KS_eigenvector, &
             KS_eigenvector_complex, &
             en_se, en_rse &
            )
        call get_times(time_rse_corr, clock_time_rse_corr)

         call rpa_calculation_p0(rpa_tot_en)
         if(myid.eq.0) then
            write(use_unit,'(2X,A,f21.8,A,f21.8,A)') &
             "  RPA+SE total energy          : ",  rpa_tot_en+en_se, &
             " Ha,", (rpa_tot_en+en_se)*hartree, " eV"
             write(use_unit,'(2X,A,f21.8,A,f21.8,A)') &
             "  RPA+rSE total energy         : ",  rpa_tot_en+en_rse, &
             " Ha,", (rpa_tot_en+en_rse)*hartree, " eV"
             write(use_unit,*)
         endif
call perfoff
       elseif(use_ci) then  ! add by igor.
           call ci_calculation_slab()

       elseif(use_mp2 .or. use_dftpt2) then  ! add by igor.
           call pt2_calculation(pt2_tot_en)

       elseif(use_periodic_gw) then
           allocate(n_homo_k(n_spin, n_k_points),stat=info)
           call check_allocation(info,'n_homo_k',func) 
!  determine the maximal HOMO for each spin channel
           do i_spin = 1, n_spin, 1
             do i_k_point = 1, n_k_points, 1
               do i_state = 1, n_states, 1
                 if(occ_numbers(i_state,i_spin,i_k_point).gt.1.e-6) then
                    n_homo_k(i_spin, i_k_point) = i_state
                 endif
               enddo
               if(n_homo(i_spin) .lt. n_homo_k(i_spin,i_k_point)) &
                  n_homo(i_spin) = n_homo_k(i_spin,i_k_point)
             enddo
           enddo
! determine lower and upper limits of the states to be treated
           if(n_low_state .ge. int(n_electrons/2)) then
              n_low_state = 1
           endif
           if (n_high_state .le. int(n_electrons/2)+1) then
              n_high_state = int(n_electrons/2) + 10
           endif
           if (n_high_state.gt.n_states) then
              n_high_state=n_states
           elseif (n_high_state.lt.max(n_homo(1),n_homo(n_spin))) then
              n_high_state = max(n_homo(1),n_homo(n_spin))
           endif

           if (flag_frozen_core_postSCF) then ! count the frozen core states
              call count_frozen_core_states(n_low_state_polar)
           else
              n_low_state_polar = 1
           endif
           n_low_state = max(n_low_state, n_low_state_polar)


           if(.not.allocated(gw_selfenergy)) then
             allocate(gw_selfenergy(n_freq,n_low_state:n_high_state,n_spin,n_irk_points_task),stat=info)
             call check_allocation(info,'gw_self_energy',func)
           endif
          
           if(out_band) then
              n_band_kpoints=sum(n_points_in_band(:))
           else
              n_band_kpoints=1
           endif

           n_band_kpoints_task = 0
           i_k_point_index = 0
!           do i_band = 1, n_plot_band, 1
             do i_k_point = 1, n_band_kpoints, 1
               i_k_point_index = i_k_point_index + 1
               if(myid ==  MOD(i_k_point_index, n_tasks) .and. myid <= n_band_kpoints )then
                  n_band_kpoints_task = n_band_kpoints_task + 1
               endif
             enddo
!           end do

           if(out_band) then

!  hamiltonian is needed for generating KS eigenvectors on new k point mesh.
             call integrate_real_hamiltonian_matrix_p2 &
                  ( hartree_potential, rho, rho_gradient, kinetic_density, partition_tab, &
                    l_shell_max, en_xc, en_pot_xc, hamiltonian, en_vdw, en_pot_vdw )

             allocate(gw_selfe_band(n_freq,n_low_state:n_high_state,n_spin,n_band_kpoints_task),stat=info)
             call check_allocation(info,'gw_selfe_band',func)
             call initialize_tb_auxmat(1, OVLP_TYPE_COULOMB)
           endif

! dielec_func_imagfreq defined in module gw_para
           if(.not. allocated(dielec_func_imagfreq)) then
              allocate(dielec_func_imagfreq(n_full_freq),stat=info)
              call check_allocation(info, 'dielec_func_imagfreq            ')
           endif

! dielec_func_imagfreq is calculated here.
           if(use_scalapack) then
              call KS_dielectric_calculation_imagfreq( mxld, mxcol, n_spin, n_states, KS_eigenvalue, eigenvec, &
                     eigenvec_complex, occ_numbers, chemical_potential, partition_tab, l_shell_max )

           else
               call KS_dielectric_calculation_imagfreq( n_basis, n_states, n_spin, n_states, KS_eigenvalue, KS_eigenvector, &
                      KS_eigenvector_complex, occ_numbers, chemical_potential, partition_tab, l_shell_max )
           endif
           call evaluate_periodic_gw_selfenergy &
             ( n_low_state_polar, n_low_state, n_high_state, &
                occ_numbers, n_freq, n_full_freq, &
                omega_grid, omega_full_grid, womega_full, &
                chemical_potential_spin, dielec_func_imagfreq, &
                KS_eigenvalue, KS_eigenvector, KS_eigenvector_complex, &
                KS_eigenvector_irk, KS_eigenvector_complex_irk, &
                gw_selfenergy, gw_selfe_band, out_self_energy &
              )

           if (allocated (dielec_func_imagfreq)) then
             deallocate (dielec_func_imagfreq)
           endif


!            if (use_threadsafe_gwinit) then
!              call my_cleanup_lvl_triples()
!            else
!              call cleanup_lvl_triples()
!            endif

           call deallocate_tb_auxmat()

           if(out_gw_regular_kgrid) then

              if(.not.allocated(sigma_par_p0)) then
                allocate(sigma_par_p0(n_max_par,n_low_state:n_high_state,n_spin,n_irk_points_task))
              endif

              call analy_continue_self_energy_p0 &
              ( anacon_type, n_irk_points, n_irk_points_task, & 
                n_max_par, n_low_state,n_high_state,n_freq, &
                omega_grid, gw_selfenergy, &
                sigma_par_p0)
            endif

           if(.not.allocated(xc_realspace)) then
            allocate(xc_realspace(n_hamiltonian_matrix_size,n_spin))
           endif
           if(.not.allocated(x_realspace)) then
            allocate(x_realspace(n_hamiltonian_matrix_size,n_spin))
           endif
           if(.not.allocated(c_realspace)) then
            allocate(c_realspace(n_hamiltonian_matrix_size,n_spin))
           endif
           call integrate_xc_realspace_p2 &
           (hartree_potential, rho, rho_gradient,  &
            kinetic_density, &
            partition_tab, l_shell_max, en_xc, en_pot_xc, &
            xc_realspace, &
            x_realspace, &
            c_realspace &
           )

           if(out_gw_regular_kgrid) then

             if(.not.allocated(exact_x_kspace)) then
               allocate(exact_x_kspace(n_low_state:n_high_state,n_spin,n_irk_points_task))
             endif
             if(.not.allocated(xc_kspace)) then
               allocate(xc_kspace(n_low_state:n_high_state,n_spin,n_irk_points_task))
             endif

             call evaluate_exx_matr_kspace &
             (n_irk_points, n_irk_points_task, n_low_state, n_high_state, &
              KS_eigenvector_irk, KS_eigenvector_complex_irk, &
              exact_x_kspace &
              )

             call evaluate_xc_matr_kspace &
             (n_irk_points, n_irk_points_task, n_low_state, n_high_state, &
              KS_eigenvector_irk, KS_eigenvector_complex_irk, &
              xc_realspace, xc_kspace &
              )

              if(.not.allocated(qp_energy)) then
                allocate(qp_energy(n_low_state:n_high_state,n_spin,n_irk_points_task))
              endif
              call quasi_particle_energy_p0 &
              (anacon_type, n_max_par, &
               n_low_state, n_high_state, &
               n_freq, omega_grid, &
               sigma_par_p0, occ_numbers, KS_eigenvalue, &
               chemical_potential_spin, &
               exact_x_kspace, xc_kspace, &
               qp_energy )
! end of  if(out_gw_regular_kgrid)
           endif

           if (out_band) then

             call get_gw_band_struct_info &
               (xc_realspace, gw_selfe_band)

          endif
          call gw_cleanup_lvl_tricoeff_recip
           if(allocated(xc_realspace)) then
             deallocate(xc_realspace)
           endif
           if(allocated(x_realspace)) then
             deallocate(x_realspace)
           endif
           if(allocated(c_realspace)) then
             deallocate(c_realspace)
           endif
           if(allocated(exact_x_kspace)) then
             deallocate(exact_x_kspace)
           endif
           if(allocated(xc_kspace)) then
             deallocate(xc_kspace)
           endif
           if(allocated(gw_selfenergy)) then
             deallocate(gw_selfenergy)
           endif
           if(allocated(gw_selfe_band)) then
             deallocate(gw_selfe_band)
           endif
           if(allocated(sigma_par_p0)) then
             deallocate(sigma_par_p0)
           endif
           if(allocated(qp_energy)) then
             deallocate(qp_energy)
           endif
           if(allocated(kq_point_list)) then
             deallocate(kq_point_list)
           endif
           if(allocated(kpq_point_list)) then
             deallocate(kpq_point_list)
           endif
           if(allocated(KS_eigenvector_irk)) then
             deallocate(KS_eigenvector_irk)
           endif
           if(allocated(KS_eigenvector_complex_irk)) then
             deallocate(KS_eigenvector_complex_irk)
           endif
           if(allocated(n_states_k)) then
             deallocate(n_states_k)
           endif
 
       endif

       call deallocate_gw()

       deallocate(lvl_tricoeff_mod_r,coulomb_matr_blacs,coulomb_cut_blacs,coulomb_gamma_blacs)

       call perfoff
       call perfout('pscf')
      end subroutine post_scf_correlation_treatment_p0
!***************
