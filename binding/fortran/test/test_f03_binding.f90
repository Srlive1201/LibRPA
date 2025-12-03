program test_f03_binding

   use mpi
   use librpa_f03
   implicit none

   ! Example dimensions
   integer, parameter :: nspins = 1
   integer, parameter :: nkpts = 2
   integer, parameter :: nstates = 3
   integer, parameter :: nbasis = nstates

   integer :: ierr, ispin, ik, istate
   integer :: provided
   type(LibrpaHandler) :: h
   type(LibrpaOptions) :: opts
   real*8, allocatable :: wg(:,:,:), ekb(:,:,:)
   real*8 :: efermi

   call MPI_Init_thread(MPI_THREAD_MULTIPLE, provided, ierr);
   if (MPI_THREAD_MULTIPLE .ne. provided) then
      write(*,*) "Warning: MPI_Init_thread provide ", provided, " != required", MPI_THREAD_MULTIPLE
   endif
   call librpa_init_global(.true., "output.txt", .false.)

   write(*,*) "Host code output: initialize default options"
   call opts%init()
   write(*,*) "Default frequency numbers:", opts%nfreq

   call h%create(MPI_COMM_WORLD)

   ! Play with the handler
   ! Initialize dimensions
   call h%set_scf_dimension(nspins, nkpts, nstates, nbasis)
   ! Parse the energy levels and occupation numbers
   allocate(wg(nstates, nkpts, nspins), ekb(nstates, nkpts, nspins))
   ! only the first state is occupied
   wg(:, :, :) = 0.0d0
   wg(1, :, :) = 2.0d0
   efermi = 0.0d0
   do ispin = 1, nspins
      do ik = 1, nkpts
         do istate = 1, nstates
            ekb(istate, ik, ispin) = efermi + sign(istate + ik - 0.5d0, istate - 1.5d0)
         end do
      end do
   end do
   print*, wg
   print*, ekb
   call h%set_wg_ekb_efermi(nspins, nkpts, nstates, wg, ekb, efermi)

   call h%destroy()

   deallocate(wg, ekb)

   call librpa_finalize_global
   call MPI_Finalize(ierr)

end program
