program test_f03_binding_call_internal_test

   use mpi
   use librpa_f03
   implicit none

   character(len=*), parameter :: prog = "test_f03_binding_call_internal_test"
   integer :: ierr

   call initialize()

   call librpa_test()

   call mpi_barrier(MPI_COMM_WORLD, ierr)

   call finalize(.true.)

contains

   subroutine initialize
      implicit none
      integer :: provided
      call MPI_Init_thread(MPI_THREAD_MULTIPLE, provided, ierr)
      if (MPI_THREAD_MULTIPLE .ne. provided) then
         write(*,*) "Warning: MPI_Init_thread provide ", provided, " != required", MPI_THREAD_MULTIPLE
      endif
      call librpa_init_global(.false., "output.txt", .true.)
   end subroutine initialize

   subroutine finalize(success)
      logical, intent(in) :: success

      if (success) then
         write(*,*) prog, " succeeded"
      else
         write(*,*) prog, " failed"
      end if

      call librpa_finalize_global
      call MPI_Finalize(ierr)
   end subroutine finalize

end program

