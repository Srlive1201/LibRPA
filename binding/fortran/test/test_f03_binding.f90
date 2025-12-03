program test_f03_binding

   use mpi
   use librpa_f03
   implicit none

   integer :: ierr
   integer :: provided
   type(LibrpaHandler) :: h
   type(LibrpaOptions) :: opts

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
   call h%destroy()


   call librpa_finalize_global
   call MPI_Finalize(ierr)

end program
