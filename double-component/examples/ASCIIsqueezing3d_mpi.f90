program ASCIIsqueezing3d_mpi

   use mpi
   implicit none
   integer fh, ierr, i
   double precision en
   character(len=*), parameter :: file = "Energy3D_M=5_ua=100p4_ub=100p4_uab=95p0.bin"

   call MPI_INIT(ierr)

   call MPI_FILE_OPEN(MPI_COMM_WORLD, file, MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ierr)

   do i = 0, 127
      call MPI_FILE_READ(fh, en, 1, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)
      print *, en
   enddo

   call MPI_FILE_CLOSE(fh, ierr)

   call MPI_FINALIZE(ierr)

end program
