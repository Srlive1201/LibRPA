!****s* FHI-aims/evaluate_crpa_energy_kspace
!  NAME
!   evaluate_crpa_energy_kspace
!  SYNOPSIS

      subroutine evaluate_crpa_energy_kspace &
           (n_low_state, occ_numbers, n_full_freq, &
            omega_full, womega_full, &
            KS_eigenvalue, KS_eigenvector, &
            KS_eigenvector_complex, &
            rpa_c_energy &
           )

!  PURPOSE
!  Subroutine evaluate_crpa_energy_kspace evaluates the correlation
!  energy at the RPA level using the adiabatic connection fluctuation
!  dissipation theorem.
!
!  E_RPA = 1/2pi \int dw { ln(det(1-v_times_polar)) + tr(v_times_polar) }

! USES
      use dimensions
      use prodbas
      use pbc_lists
      use hartree_fock
      use constants
      use mpi_tasks
      use synchronize_mpi
      use timing
      use runtime_choices
      use crpa_blacs
      use evaluate_polarisability_kspace_mod
      implicit none

! ARGUMENTS 

      integer :: n_full_freq
      integer :: n_low_state
      integer :: n_high_state

      real*8  :: occ_numbers(n_states,n_spin,n_k_points)
      real*8  :: omega_full(n_full_freq)
      real*8  :: womega_full(n_full_freq)
      real*8  :: KS_eigenvalue(n_states,n_spin,n_k_points)
      real*8  :: KS_eigenvector(n_basis,n_states,n_spin,n_k_points_task)
      complex*16  :: KS_eigenvector_complex(n_basis,n_states,n_spin,n_k_points_task)

!     output
      real*8  :: rpa_c_energy
      complex*16, external:: pzlatra
      real*8, external:: pdlatra

! INPUTS
! o  n_full_freq -- integer number,
!            the number of frequency points for the screened Coulomb interaction W
! o  n_low_state  -- integer number,
!            the lowest KS/HF eigenstate taken into account in the polarisability calcua            ltions
! o  n_high_state -- integer number,
!            the highest KS/HF eigenstate. In the present case, n_high_state >= n_homo
!            should be fine. 
! o  n_electrons -- real number
!            the total number of electrons in the system
! o  occ_numbers -- real 2-dimentianal array of length (n_states, n_spin)
!            the occupation number of the electrons for each eigenstate and each spin
! o  omega_full(n_freq) -- real array
!            the Gauss-Legendre frequency grid for the screened Coulomb interaction
! o  womega_full(n_freq) -- real array
!            the weigth of the Gauss-Legendre frequency grid for the screened Coulomb 
!            in teraction
! o  chemical_potential -- real number, the chemical potential of the system
! o  KS_eigenvalue -- real array,
!            the eigenvalues of the single-particle calculation. For DFT calculation,
!            this is the KS eigenvalue, but for HF calculation, this is then the HF
!            eigenvalue
! o  KS_eigenvector -- real array,
!            the eigenvector of the single-particle calculation
! o  KS_eigenvector_complex -- complex array,
!            the complex eigenvector of the single-particle calculation,
!            used when "real_eigenvectors == .false."
!           
!
! OUTPUT
! o  rpa_c_energy -- real number, the calculated RPA correlation energy
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

      complex*16  det_v_times_polar
      complex*16  trace_v_times_polar
      complex*16  trace_v(n_irk_points), trace_polar(n_irk_points)
      real*8  rpa_c_integrand
      real*8, dimension(:), allocatable :: rpa_c_kgrid
      
!    local timing
      real*8  temp_time_rpa
      real*8  temp_clock_time_rpa

!     auxiliary matrices for Level 3 Blas matrix multiplications
!     n_first : the first orbital which is NOT fully occupied
      integer :: n_first(n_spin)
      integer :: info

      complex*16, dimension(:,:,:), allocatable :: polar_kspace_complex
      real*8,dimension(:,:,:), allocatable:: polar_kspace_real
      real*8,dimension(:,:), allocatable:: coulomb_matr_blacs_real
      complex*16, dimension(:,:), allocatable :: v_times_polar ,v_times_polar_t, v_times_polar_2
      real*8, dimension(:,:), allocatable :: rv_times_polar, rv_times_polar_t, rv_times_polar_2

      character*50  filename

!     counters

      integer :: i_state
      integer :: i_freq
      integer :: i_spin
      integer :: i_index
      integer :: i_k_point
      integer :: i_k_point_local
      integer :: i_irk_point
      integer :: i_irk_point_local
      integer :: i_prodbas_1
      integer :: max_irk_points_task
      integer:: n_freqs_block, i_fblock, lbf, ubf, n_blocks
      integer:: i_rc,j_rc,mu,nu,two_nbas
!      comm
      real*8  temp_crpa_grid
      real*8 :: k_vec(3)

      integer :: mpierr

!     begin work

      call perfon('ev_rpa')
      if(myid.eq.0) then
        write(use_unit,*)
        write(use_unit,*)"-----------------------------------------------------------------------"
        write(use_unit,'(2X,A)') &
              "Start to calculate the periodic RPA correlation energy  ... "
      endif

      if(flag_KS_eigenfunc_conjg) then
         KS_eigenvector_complex = conjg(KS_eigenvector_complex)
      endif

!     determine the highest occupied orbital level
!     such complication occurs when some of the orbitals are
!     either not fully occupied or not fully empty
      n_first(:) = 1
      do i_k_point = 1, n_k_points, 1
        do i_spin = 1, n_spin
          do i_state = 1, n_states
           if (abs(occ_numbers(i_state,i_spin,i_k_point)-dble(2/n_spin)) &
                           .lt.1.d-8) then
             n_first(i_spin)= i_state + 1
           endif
          enddo
          if(n_first(i_spin) .gt. n_states) then
           n_first(i_spin) = n_states
          endif
        enddo
      enddo

! work array

      allocate(rpa_c_kgrid(n_irk_points),stat=i_index)

      time_rpa_corr = 0.d0
      clock_time_rpa_corr = 0.d0
      time_polar = 0.d0
      clock_time_polar = 0.d0

      rpa_c_energy = 0.d0 
      rpa_c_kgrid = 0.d0 

      n_low_state=max(1,n_low_state)
      if (flag_frozen_core_postSCF) then ! count the frozen core states
          call count_frozen_core_states(n_low_state)
      endif
      if (myid .eq. 0) then
          write(use_unit,'(2X,A,I12)') &
              'The first valence state in the frozen-core algorithm :', &
              n_low_state
      endif

      !because of the global mpi_fence operations, all mpi tasks have to run all iteration so of the irkq loop
      if (mod(n_irk_points,n_tasks_irkq).eq.0) then    
         max_irk_points_task=n_irk_points/n_tasks_irkq
      else
         max_irk_points_task=n_irk_points/n_tasks_irkq+1
      end if

      !split computation into several blocks if memory consumption of polar_kspace becomes too big
      if(n_k_points.eq.1) then
         !polar_kspace is double precision (8 bytes per element)
         n_freqs_block=min(80000000/(bb_bl_row*bb_bl_col),n_full_freq)

         !take into account coulomb_matr_blacs
         if ((n_freqs_block+2)*bb_bl_row*bb_bl_col.gt.80000000) then
            n_freqs_block=n_freqs_block-2
         end if
      else
         !polar_kspace is double complex (16 bytes per element)
         n_freqs_block=min(40000000/(bb_bl_row*bb_bl_col),n_full_freq)

         !take into account coulomb_matr_blacs
         if ((n_freqs_block+1)*bb_bl_row*bb_bl_col.gt.40000000) then
            n_freqs_block=n_freqs_block-1
         end if
      end if
      n_blocks=n_full_freq/n_freqs_block
      if(n_blocks*n_freqs_block.lt.n_full_freq) n_blocks=n_blocks+1
      if(myid.eq.0) print*,'using n_freqs_block=',n_freqs_block, ', n_blocks=',n_blocks

      if (real_eigenvectors) then
         allocate(polar_kspace_real(lbb_row:ubb_row, lbb_col:ubb_col, n_freqs_block),stat=i_index)
      else
         allocate(polar_kspace_real(0, lbb_col:ubb_col, n_freqs_block),stat=i_index)
      end if
      call check_allocation(i_index, 'polar_kspace_real                    ')

      if(n_k_points.gt.1) then
         allocate(polar_kspace_complex(lbb_row:ubb_row, lbb_col:ubb_col, n_freqs_block),stat=i_index)
      else
         allocate(polar_kspace_complex(1,1,1),stat=i_index)
      end if
      call check_allocation(i_index, 'polar_kspace_complex                    ')
      open(11,file="real_coulomb.txt")
      write(11,*)max_irk_points_task
      do i_irk_point_local = 1, max_irk_points_task
         i_irk_point=n_tasks_irkq*(i_irk_point_local-1) + myid_irkq + 1
         do i_fblock=1,n_blocks
!         do i_fblock=1,1
            lbf=(i_fblock-1)*n_freqs_block+1
            ubf=min(i_fblock*n_freqs_block,n_full_freq)

            call get_timestamps(temp_time_rpa, temp_clock_time_rpa)
            call evaluate_polarisability_kspace_list &
                 (n_low_state, i_irk_point_local,ubf-lbf+1, omega_full(lbf:ubf), KS_eigenvalue, KS_eigenvector, &
                 KS_eigenvector_complex, occ_numbers, polar_kspace_real(:,:,1:ubf-lbf+1), polar_kspace_complex(:,:,1:ubf-lbf+1))

            call get_timestamps(rtime, clock_rtime)
            time_polar = time_polar + rtime - temp_time_rpa
            clock_time_polar = clock_time_polar + clock_rtime - temp_clock_time_rpa
            
            
            if(i_irk_point_local.gt.n_irkq_points_task) cycle
            !print*,"IN real_coulomb k: ",i_irk_point_local,i_irk_point
           ! k_vec = matmul(recip_lattice_vector,k_point_list(i_irk_point,:))
           ! print*,k_vec
         !this seems to hold, but i don't know why and use real_eigenvectors instead
!         i_k_point=inv_irk_point_mapping(i_irk_point)
!         if (i_k_point.eq.kq_point_list(1,i_k_point)) then
            write(11,*)size(coulomb_matr_blacs(:,:,i_irk_point_local)),size(coulomb_matr_blacs(:,:,i_irk_point_local),dim=1),size(coulomb_matr_blacs(:,:,i_irk_point_local),dim=2)
            write(11,*)inv_irk_point_mapping(i_irk_point_local),irk_weight(i_irk_point)
               !print*, i_irk_point_local
            do i_rc=1,size(coulomb_matr_blacs(:,:,i_irk_point_local),dim=1)
               do j_rc = 1,size(coulomb_matr_blacs(:,:,i_irk_point_local),dim=2)
                  write(11,'(2F20.15)')coulomb_matr_blacs(i_rc,j_rc,i_irk_point_local)
                  !print*, coulomb_matr_blacs_real(i_rc,j_rc)
               end do
            end do
            if(real_eigenvectors) then 
               !print*,"   go real_eigenvectors"
               allocate(coulomb_matr_blacs_real(lbb_row:ubb_row, lbb_col:ubb_col),stat=i_index)
               allocate(rv_times_polar(lbb_row:ubb_row, lbb_col:ubb_col),stat=i_index)
               allocate(rv_times_polar_t(lbb_row:ubb_row, lbb_col:ubb_col),stat=i_index)

               call check_allocation(i_index, 'v_times_polar                   ')
               coulomb_matr_blacs_real=real(coulomb_matr_blacs(:,:,i_irk_point_local))

               
               
               
               do i_freq = lbf, ubf
                  if (n_k_points.gt.1) polar_kspace_real(:,:,i_freq-lbf+1)=real(polar_kspace_complex(:,:,i_freq-lbf+1))

                  if(i_freq == 1) then
                     trace_polar(i_irk_point)=pdlatra(n_basbas,polar_kspace_real,1,1,bb2desc)
                     trace_v(i_irk_point)=pdlatra(n_basbas,coulomb_matr_blacs_real,1,1,bb2desc)
                  endif
                  !            call perfon('rdge')
                  !  Multiply \chi_0 with bare coulomb matrix v.
                  call pdgemm('N', 'N', n_basbas, n_basbas, n_basbas, 1.d0, &
                       coulomb_matr_blacs_real(lbb_row,lbb_col), 1, 1, bb2desc, &
                       polar_kspace_real(lbb_row,lbb_col,i_freq-lbf+1), 1, 1, bb2desc, 0.d0, &
                       rv_times_polar(lbb_row, lbb_col), 1, 1, bb2desc)
                  
                  !            call perfoff
                  trace_v_times_polar=pdlatra(n_basbas,rv_times_polar,1,1,bb2desc)
                  
                  do i_prodbas_1 = lbb_row, ubb_row
                     if((i_prodbas_1.ge.lbb_col).and.(i_prodbas_1.le.ubb_col)) then
                        rv_times_polar(i_prodbas_1, i_prodbas_1) = rv_times_polar(i_prodbas_1, i_prodbas_1) - 1.
                     end if
                  enddo
                  rv_times_polar(:,:) = - rv_times_polar(:,:)
                  
                  !symmetrize for Cholesky factorization
                  call pdtran( n_basbas, n_basbas, &
                       1.d0, rv_times_polar, 1, 1, bb2desc, &
                       0.d0, rv_times_polar_t, 1, 1, bb2desc )
                  call pdgemm('N', 'N', n_basbas, n_basbas, n_basbas, 1.d0, &
                       rv_times_polar, 1, 1, bb2desc, &
                       rv_times_polar_t, 1, 1, bb2desc, 0.d0, &
                       polar_kspace_real(lbb_row,lbb_col,i_freq-lbf+1), 1, 1, bb2desc)
                  
                  !factorize
                  call  pdpotrf ('L',n_basbas,polar_kspace_real(lbb_row,lbb_col,i_freq-lbf+1),1,1,bb2desc,info)
                  
                  det_v_times_polar = (1.d0,0.d0)
                  do i_prodbas_1 = lbb_row, ubb_row
                     if((i_prodbas_1.ge.lbb_col).and.(i_prodbas_1.le.ubb_col)) then
                        det_v_times_polar = det_v_times_polar * polar_kspace_real(i_prodbas_1,i_prodbas_1,i_freq-lbf+1)
                     end if
                  enddo
                  if (irkblacs_member) call mpi_allreduce(MPI_IN_PLACE,det_v_times_polar,1,MPI_DOUBLE_COMPLEX,&
                       MPI_PROD,comm_blacs,mpierr)                  
                  
                  if ((real(det_v_times_polar).lt.0).and.(myid_bl.eq.0)) write (use_unit,'(2X,A)') &
                       'WARNING: det(1-chi_0*V) in evaluate_crpa_energy_kspace is negative! Please contact the developers.'
                  
                  rpa_c_integrand = log (abs(real(det_v_times_polar))) + &
                       real(trace_v_times_polar)
                  
                  rpa_c_energy = rpa_c_energy + &
                       rpa_c_integrand * womega_full(i_freq) * &
                       irk_weight(i_irk_point)
                  
                  rpa_c_kgrid(i_irk_point) = rpa_c_kgrid(i_irk_point) + &
                       rpa_c_integrand * womega_full(i_freq) 
                  ! end of loop over i_freq
                  
               enddo
               !temp_crpa_grid=rpa_c_kgrid(i_irk_point)*irk_weight(i_irk_point)/2.d0/pi
               !print*,"cRPA_kgrid",i_irk_point,rpa_c_kgrid(i_irk_point),temp_crpa_grid,irk_weight(i_irk_point)
               deallocate(coulomb_matr_blacs_real,rv_times_polar,rv_times_polar_t)
            else
               !print*,"   go complex_eigenvectors"
               allocate(v_times_polar(lbb_row:ubb_row, lbb_col:ubb_col),stat=i_index)
               allocate(v_times_polar_t(lbb_row:ubb_row, lbb_col:ubb_col),stat=i_index)
               allocate(v_times_polar_2(lbb_row:ubb_row, lbb_col:ubb_col),stat=i_index)
               call check_allocation(i_index, 'v_times_polar                   ')

              ! open(12,file='pi_mat.txt',position='append')
              ! write(12,*)i_irk_point_local,irk_weight(i_irk_point)
               do i_freq = lbf, ubf
                  if(i_freq == 1) then
                     trace_polar(i_irk_point)=pzlatra(n_basbas,polar_kspace_complex,1,1,bb2desc)
                     trace_v(i_irk_point)=pzlatra(n_basbas,coulomb_matr_blacs(lbb_row,lbb_col,i_irk_point_local),1,1,bb2desc)
                  endif
                  
                  !  Multiply \chi_0 with bare coulomb matrix v.
                  
                  call pzgemm('N', 'N', n_basbas, n_basbas, n_basbas, (1.d0,0.d0), &
                       coulomb_matr_blacs(lbb_row,lbb_col,i_irk_point_local), 1, 1, bb2desc, &
                       polar_kspace_complex(lbb_row,lbb_col,i_freq-lbf+1), 1, 1, bb2desc, (0.d0, 0.d0), &
                       v_times_polar(lbb_row, lbb_col), 1, 1, bb2desc)
                  
                  ! if(i_freq==6) then
                  !    write(12,*)"coulomb_mat"
                  !    write(12,200) (((coulomb_matr_blacs(mu,nu,i_irk_point_local)),nu=1,n_basbas),mu=1,n_basbas)
                  !    write(12,*)"chi0_mat"
                  !    write(12,200)(((polar_kspace_complex(mu,nu,i_freq-lbf+1)),nu=1,n_basbas),mu=1,n_basbas)
                  !    write(12,*)"pi_mat"
                  !    write(12,200)(((v_times_polar(mu,nu)),nu=1,n_basbas),mu=1,n_basbas)
                  ! end if
                  !200 FORMAT (1X, 52F10.6)
                     trace_v_times_polar=pzlatra(n_basbas,v_times_polar,1,1,bb2desc)
                  
                  do i_prodbas_1 = lbb_row, ubb_row
                     if((i_prodbas_1.ge.lbb_col).and.(i_prodbas_1.le.ubb_col)) then
                        v_times_polar(i_prodbas_1, i_prodbas_1) = v_times_polar(i_prodbas_1, i_prodbas_1) - (1.d0,0.d0)
                     end if
                  enddo
                  v_times_polar(:,:) = - v_times_polar(:,:)
                  
                  !symmetrize to get self-adjoint matrix for Cholesky factorization
                  call pztranc( n_basbas, n_basbas, &
                       (1.d0,0.d0), v_times_polar, 1, 1, bb2desc, &
                       (0.d0,0.d0), v_times_polar_t, 1, 1, bb2desc )
                  call pzgemm('N', 'N', n_basbas, n_basbas, n_basbas, (1.d0,0.d0), &
                       v_times_polar, 1, 1, bb2desc, &
                       v_times_polar_t, 1, 1, bb2desc, (0.d0, 0.d0), &
                       v_times_polar_2, 1, 1, bb2desc)
                  
                  !Cholesky factorization               
                  call  pzpotrf ('L',n_basbas,v_times_polar_2,1,1,bb2desc,info)
                  
                  det_v_times_polar = (1.d0,0.d0)
                  do i_prodbas_1 = lbb_row, ubb_row
                     if((i_prodbas_1.ge.lbb_col).and.(i_prodbas_1.le.ubb_col)) then
                        det_v_times_polar = det_v_times_polar * v_times_polar_2(i_prodbas_1,i_prodbas_1)
                     end if
                  enddo
                  
                  if (irkblacs_member) call mpi_allreduce(MPI_IN_PLACE,det_v_times_polar,1,MPI_DOUBLE_COMPLEX,&
                       MPI_PROD,comm_blacs,mpierr)
                  
                  if ((real(det_v_times_polar).lt.0).and.(myid_bl.eq.0)) write (use_unit,'(2X,A)') &
                       'WARNING: det(1-chi_0*V) in evaluate_crpa_energy_kspace is negative! Please contact the developers.'
                  
                  
                  rpa_c_integrand = log (abs(real(det_v_times_polar))) + &
                       real(trace_v_times_polar)
                  
                  rpa_c_energy = rpa_c_energy + &
                       rpa_c_integrand * womega_full(i_freq) * &
                       irk_weight(i_irk_point)
                  
                  rpa_c_kgrid(i_irk_point) = rpa_c_kgrid(i_irk_point) + &
                       rpa_c_integrand * womega_full(i_freq) 
                  ! end of loop over i_freq
               enddo
               !close(12)
               !temp_crpa_grid=rpa_c_kgrid(i_irk_point)*irk_weight(i_irk_point)/2.d0/pi
               !print*,"cRPA_kgrid",i_irk_point,rpa_c_kgrid(i_irk_point),temp_crpa_grid,irk_weight(i_irk_point)
               deallocate (v_times_polar,v_times_polar_t,v_times_polar_2)
            end if
            
            call get_timestamps(temp_time_rpa, temp_clock_time_rpa )
            time_rpa_corr = time_rpa_corr +  temp_time_rpa - rtime
            clock_time_rpa_corr = clock_time_rpa_corr + temp_clock_time_rpa - clock_rtime
            ! end of loop over i_irk_point_local
         enddo
      end do
      close(11)
      call sync_timing(time_polar)
      call sync_timing(time_rpa_corr)

      if (irkblacs_member) then
         call sync_vector(rpa_c_kgrid,n_irk_points,comm_irkq)
         call sync_vector_complex(trace_v,n_irk_points,comm_irkq)
         call sync_vector_complex(trace_polar,n_irk_points,comm_irkq)
      end if
      if(myid.eq.0) then
        do i_k_point = 1, n_k_points, 1
            if(.not. irk_point_included(i_k_point) ) cycle

            i_irk_point = irk_point_mapping(i_k_point)
            write(use_unit,'(2I6,4f18.8)') i_k_point, i_irk_point, k_point_list(i_k_point,:), &
                  rpa_c_kgrid(i_irk_point)/2.d0
            write(use_unit,'(4f18.8)') trace_v(i_irk_point), trace_polar(i_irk_point)
        enddo
      endif

!      call sync_real_number(rpa_c_energy)     
      if (irkblacs_member) call mpi_allreduce(MPI_IN_PLACE,rpa_c_energy,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm_irkq,mpierr)
      rpa_c_energy=rpa_c_energy/2.d0/pi
!      rpa_c_energy= - rpa_c_energy/2.d0

! Delta term correction
! This deals with the case when there are frational occupations and 
! the intra-orbtial excitations should be taken into account. The contribution
! is a delta function at zero frequency and a finite contribution when integrated
! out over frequency axis.

      if(myid.eq.0) then
        write(use_unit,*)
        write(use_unit,*)"----------------------------------------------------", &
                  "-------------------------"
        write(use_unit,'(2X,A,2X,f19.8,2X,A,f19.8,2X,A)') &
            " RPA correlation energy :", rpa_c_energy, "Ha,", &
             rpa_c_energy*hartree, "eV"
!        write(use_unit,*)"----------------------------------------------------", &
!                 "-------------------------"
        write(use_unit,*)
      endif

      deallocate (polar_kspace_real)      
      deallocate (polar_kspace_complex)

      call perfoff

      return

      end subroutine evaluate_crpa_energy_kspace
