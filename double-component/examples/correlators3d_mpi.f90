program correlator3d_mpi

use mpi
use solvers
use measures
use parameters, only:dp

implicit none

! for MPI processes
integer ierr, comsize, rank, Num
integer fh
integer i, num_units
integer(kind=MPI_OFFSET_KIND) disp

character(len=200) :: file
integer L, nmax
real( kind = dp ) ja, jb, ua, ub, uab, meaBase, mebBase
real( kind = dp ) dmean, jmax, dj
real( kind = dp ) :: corr2(3), corr3(4)
complex( kind = dp ), allocatable :: f(:,:,:,:,:)

! namelists
namelist /BHparams/ L, nmax, ua, ub, uab
namelist /TMMparams/ dmean, jmax, dj, file

open(138, file='TMM_inp.nml')
read(138, BHparams)
read(138, TMMparams)
close(138)

allocate( f(0:nmax,0:nmax,1:L,1:L,1:L) )

! redefine interaction energies
ub  = ub/ua
uab = uab/ua
ua  = ua/ua

! mean number of particles
meaBase = 0.5_dp * real( L**3, dp )
mebBase = 0.5_dp * real( L**3, dp )

call MPI_INIT( ierr ) 
call MPI_COMM_SIZE( MPI_COMM_WORLD, comsize, ierr )
call MPI_COMM_RANK( MPI_COMM_WORLD, rank, ierr )

! Open the file
call MPI_FILE_OPEN( MPI_COMM_WORLD, trim(file), MPI_MODE_WRONLY + MPI_MODE_CREATE, &
                    MPI_INFO_NULL, fh, ierr ) 

! Total number of points
Num = int(jmax/dj)

call process_split( Num, comsize, num_units, disp, rank )
ja = jmax - real(disp, dp) * dj
jb = ja
disp = disp * 8

call MPI_FILE_SET_VIEW(fh, disp, MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION,&
                        "native", MPI_INFO_NULL, ierr)

do i = 1, num_units

   call MPI_file_write( fh, ja, 1, MPI_double_precision, &
                        MPI_status_ignore, ierr )

   call InitUniformNC( f, meaBase, mebBase )
   call GroundStateNC( f, ja, jb, ua, ub, uab, meaBase, mebBase )
   corr2 = G2correlator( f, 1, 1, 1 )
   corr3 = G3correlator( f, 1, 1, 1 )

   call MPI_file_write( fh, corr2, 3, MPI_double_precision, &
                        MPI_status_ignore, ierr )
   call MPI_file_write( fh, corr3, 4, MPI_double_precision, &
                        MPI_status_ignore, ierr )

   ja = ja - dj
   jb = jb - dj
   disp = disp + 8

enddo

call MPI_file_close( fh, ierr )
call MPI_finalize( ierr )

end program

! ---------------------------------------
! Determine spliting of the N dimensional
! data into M processes
! ---------------------------------------
subroutine process_split( N, M, num_units, index, rank )

   use mpi
   implicit none
   integer, intent(in) :: N, M, rank
   integer, intent(out) :: num_units 
   integer(kind=MPI_OFFSET_KIND) :: index

   integer avg, rem
   ! average number of units per process
   avg = N/M
   ! remaining number 
   rem = mod(N, M)

   if( rank < rem ) then
      num_units = avg + 1
      index = (avg + 1) * rank
   else
      num_units = avg
      index = (avg + 1) * rem + avg * (rank - rem) 
   endif

end subroutine
