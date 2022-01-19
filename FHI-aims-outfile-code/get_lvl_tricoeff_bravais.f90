  !******
  !----------------------------------------------------------------------------
  !****s* FHI-aims/get_lvl_tricoeff_bravais
  !  NAME
  !    get_lvl_tricoeff_bravais
  !  SYNOPSIS

  subroutine get_lvl_tricoeff_bravais(n_cells_task)

    !  PURPOSE
    !
    !    Compute the LVL triple expansion coefficients in real space (i.e., for a set
    !    of Bravais vectors which separate the two unit cells where the basis pairs
    !    live.
    !
    !  USES

    use dimensions
    use prodbas
    use sbt_overlap_tb
    use bravais
    use pbc_lists
    use geometry
    use basis
    use mpi_tasks, only: myid, n_tasks, check_allocation
    use localorb_io, only: localorb_info
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: n_cells_task
   ! real*8, intent(OUT) :: lvl_tricoeff_bravais(n_basis,n_basis,n_basbas,n_cells_task)
     ! n_cells_task: the number of unit cells (in the Born-von Karmen supercell) per task

    !  INPUTS
    !    o n_cells_task :: the number of unit cells (within the Born-von Karmen supercell) per task
    !  OUTPUTS
    !    o lvl_tricoeff_bravais: LVL triple coefficents in real space
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    integer :: i_atom_1, i_atom_2, i_atom_aux 
    integer :: i_species_1, i_species_2, i_species_aux, i_species_other
    integer :: basis_off_1, basis_off_2, n_sp_basis_1, n_sp_basis_2
    integer :: basis_off_own, basis_off_other, n_sp_basis_own, n_sp_basis_other
    integer :: bboff_1, n_spbb_1, bboff_2, n_spbb_2
    integer :: i_sp_basis_1, i_sp_basis_2
    integer :: i_cell_1, i_cell_2, i_cell_3
    integer :: i_cell, i_cell_local, n_cell_local
    integer :: i_basis_1,i_basis_2,i_basis_3,i_basis_4,i_prodbas_1,i_prodbas_2
    real*8 :: Dvec(3), Cvec(3)
    real*8 :: Rvecs(3,n_cells_task), dummy(1,1,1,1,1,1)
    real*8, allocatable :: coeff_3fn(:,:,:,:,:)
    integer :: info
    character*150 :: info_str
    character(*), parameter :: func = 'get_lvl_tricoeff_bravais'

    integer :: i_c,j_c,k_c
    real*8 :: dR
    character(len=100) cFile

    write(info_str,'(2X,A)') "Computing the triple LVL expansion coefficents in real space ..."
    call localorb_info(info_str)
    if(.not. allocated(coeff_3fn)) then
      allocate(coeff_3fn(max_n_basis_sp,max_n_basis_sp,max_n_basbas_sp,2, n_cells_task), stat=info)
      call check_allocation(info, 'coeff_3fn', func)
    endif
    
    !lvl_tricoeff_bravais(:,:,:,:) = 0.d0
    !print *, 'LVL  n_cell_tasks  ', n_cells_task
    write(cFile,*)myid
    cFile='Cs_data_'//trim(adjustl(cFile))//'.txt'
    open(1,file=cFile)
    write(1,*) n_atoms,n_cells
    do i_atom_1 = 1, n_atoms, 1
       i_species_1 = species(i_atom_1)
       basis_off_1 = atom2basis_off(i_atom_1)
       n_sp_basis_1 = sp2n_basis_sp(i_species_1)
       bboff_1 = atom2basbas_off(i_atom_1)
       n_spbb_1 = sp2n_basbas_sp(i_species_1)

       do i_atom_2 = 1, n_atoms, 1 
          i_species_2 = species(i_atom_2)
          basis_off_2 = atom2basis_off(i_atom_2)
          n_sp_basis_2 = sp2n_basis_sp(i_species_2)
          bboff_2 = atom2basbas_off(i_atom_2)
          n_spbb_2 = sp2n_basbas_sp(i_species_2)

          Dvec = coords(:, i_atom_2) - coords(:, i_atom_1)

          i_cell_local=0
          do i_cell = 1, n_cells, 1

             if(myid .ne. mod(i_cell, n_tasks)) cycle
             i_cell_local = (i_cell-1)/n_tasks + 1 

             i_cell_1 = cell_index(i_cell, 1)
             i_cell_2 = cell_index(i_cell, 2)
             i_cell_3 = cell_index(i_cell, 3)
             ! distance between two unit cells
             Cvec = matmul(lattice_vector, (/i_cell_1, i_cell_2, i_cell_3/))

             Rvecs(:, i_cell_local) = Dvec + Cvec
          end do
          n_cell_local = i_cell_local   ! Last value
          
          if(n_cell_local.gt.0) then
            call get_pairwise_coeff_3fn(i_species_1, i_species_2, n_cell_local, Rvecs, coeff_3fn, dummy, .false.)
          endif

          do i_cell = 1, n_cells, 1
             if(myid .ne. mod(i_cell, n_tasks)) cycle
             i_cell_local = (i_cell-1)/n_tasks + 1

             i_cell_1 = cell_index(i_cell, 1)
             i_cell_2 = cell_index(i_cell, 2)
             i_cell_3 = cell_index(i_cell, 3)
             !print *, 'LVL  ', i_atom_1,'   ',i_atom_2,'   ',i_cell,'   ',i_cell_local
             dR=sqrt(Rvecs(1,i_cell_local)**2+Rvecs(2,i_cell_local)**2+Rvecs(3,i_cell_local)**2)
             if(dR .gt. 16.d0) cycle
           ! write(1,*)i_atom_1,i_atom_2,i_cell_1,i_cell_2,i_cell_3, n_sp_basis_1,n_sp_basis_2,n_spbb_1,Rvecs(1,i_cell_local),Rvecs(2,i_cell_local),Rvecs(3,i_cell_local),dR
            write(1,*)i_atom_1,i_atom_2,i_cell_1,i_cell_2,i_cell_3, n_sp_basis_1,n_sp_basis_2,n_spbb_1
            if(myid.eq.mod(1,n_tasks)) then
              if((i_cell == 1) .and. (i_atom_1 .eq. i_atom_2))  then
                do i_c=1, n_sp_basis_1,1
                  do j_c=1, n_sp_basis_2,1
                    do k_c=1,n_spbb_1,1
                      write(1,'(F20.15)')coeff_3fn(i_c,j_c,k_c,1,i_cell_local)/2.0
                    end do
                  end do
                end do
              else
                do i_c=1, n_sp_basis_1,1
                  do j_c=1, n_sp_basis_2,1
                    do k_c=1,n_spbb_1,1
                      write(1,'(F20.15)') coeff_3fn(i_c,j_c,k_c,1,i_cell_local)
                    end do
                  end do
                end do
              end if
            else
              if((i_cell == 1) .and. (i_atom_1 .eq. i_atom_2)) then
                do i_c=1, n_sp_basis_1,1
                  do j_c=1, n_sp_basis_2,1
                    do k_c=1,n_spbb_1,1
                      write(1,'(F20.15)')coeff_3fn(i_c,j_c,k_c,1,i_cell_local)/2.0
                    end do
                  end do
                end do
              else
                do i_c=1, n_sp_basis_1,1
                  do j_c=1, n_sp_basis_2,1
                    do k_c=1,n_spbb_1,1
                      write(1,'(F20.15)')coeff_3fn(i_c,j_c,k_c,1,i_cell_local)
                    end do
                  end do
                end do
              end if
            end if
              
            !  lvl_tricoeff_bravais(basis_off_1+1:basis_off_1+n_sp_basis_1, &
            !  & basis_off_2+1:basis_off_2+n_sp_basis_2, bboff_1+1:bboff_1+n_spbb_1, i_cell_local) &
            !  & = coeff_3fn(1:n_sp_basis_1,1:n_sp_basis_2,1:n_spbb_1, 1, i_cell_local)
           
            ! if((i_cell == 1) .and. (i_atom_2 .ne. i_atom_1)) then
            !  lvl_tricoeff_bravais(basis_off_1+1:basis_off_1+n_sp_basis_1, &
            !  & basis_off_2+1:basis_off_2+n_sp_basis_2, bboff_2+1:bboff_2+n_spbb_2, i_cell_local) &
            !  & = coeff_3fn(1:n_sp_basis_1,1:n_sp_basis_2,1:n_spbb_2, 2, i_cell_local)
            ! endif
! end loop over i_cell
          enddo

! end loop over i_atom_2
       enddo
! end loop over i_atom_1
    enddo
    close(1)
! The on-site coefficients (which are stored in the first task) have to be halved
!    do i_cell = 1, n_cells, 1
!        if(myid .ne. mod(i_cell, n_tasks)) cycle
!        i_cell_local = (i_cell-1)/n_tasks + 1
!        write(use_unit,*) "cell_index:", cell_index(i_cell,1), cell_index(i_cell,2), cell_index(i_cell,3)
!        do i_basis_1 = 1, n_basis
!          do i_basis_2 = 1, n_basis
!            do i_prodbas_1 = 1, n_basbas, 1
!              write(use_unit,'(3I4,f16.8)') i_basis_1, i_basis_2, i_prodbas_1, &
!                  lvl_tricoeff_bravais(i_basis_1,i_basis_2,i_prodbas_1,i_cell_local)
!            enddo
!          enddo
!        enddo
!    enddo
    !print *, 'myid:  ',myid,'   mod:  ', mod(1,n_tasks)
    ! if(myid.eq.mod(1,n_tasks)) then
    !  lvl_tricoeff_bravais(:,:,:,1)=lvl_tricoeff_bravais(:,:,:,1)/2.d0
    ! endif
    ! open(2,file="Cs_full.txt")
    ! write(2,*)size(lvl_tricoeff_bravais),size(lvl_tricoeff_bravais,dim=1),size(lvl_tricoeff_bravais,dim=2),size(lvl_tricoeff_bravais,dim=3),size(lvl_tricoeff_bravais,dim=4)
    
    ! do j_c=1,size(lvl_tricoeff_bravais,dim=2)
    !   do k_c=1,size(lvl_tricoeff_bravais,dim=3)
    !     do i_c=1,size(lvl_tricoeff_bravais,dim=1)
    !       write(2,*)lvl_tricoeff_bravais(i_c,j_c,k_c,1)
    !     end do
    !   end do
    ! end do
    ! close(2)
    if(allocated(coeff_3fn)) then
      deallocate(coeff_3fn)
    endif
  end subroutine get_lvl_tricoeff_bravais
