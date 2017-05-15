program two_mode_model3d_mpi

use mpi
use solvers
use measures
use parameters, only:dp

implicit none

! for MPI processes
integer ierr, comsize, rank
integer fh
integer(kind=MPI_OFFSET_KIND) disp

character(len=100) :: file
integer L, nmax
real( kind = dp ) ja, jb, ua, ub, uab, mea, meb, meaBase, mebBase
real( kind = dp ) dmean, jmax, dj
real( kind = dp ) :: en
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

! Always use 9 processes
if( comsize < 9 ) then
   print *, "Total number of processes in smaller than 9"
   stop
endif

select case (rank)

   case( 0 )
      mea = meaBase
      meb = mebBase
      disp = 0_MPI_OFFSET_KIND
   case( 1 )
      mea = meaBase- dmean
      meb = mebBase
      disp = rank + 1
   case( 2 )
      mea = meaBase + dmean
      meb = mebBase
      disp = rank + 1
   case( 3 )
      mea = meaBase
      meb = mebBase - dmean
      disp = rank + 1
   case( 4 )
      mea = meaBase
      meb = mebBase + dmean
      disp = rank + 1
   case( 5 ) 
      mea = meaBase + dmean
      meb = mebBase - dmean
      disp = rank + 1
   case( 6 )
      mea = meaBase + dmean
      meb = mebBase + dmean
      disp = rank + 1
   case( 7 )
      mea = meaBase - dmean
      meb = mebBase - dmean
      disp = rank + 1
   case( 8 )
      mea = meaBase - dmean
      meb = mebBase + dmean
      disp = rank + 1

end select

call MPI_FILE_OPEN( MPI_COMM_WORLD, trim(file), MPI_MODE_WRONLY + MPI_MODE_CREATE, &
                    MPI_INFO_NULL, fh, ierr ) 

! ----------------------------------------------------------------
ja = jmax
jb = jmax
do while( ja > 0.0_dp )

   call InitUniformNC( f, mea, meb )
   call GroundStateNC( f, ja, jb, ua, ub, uab, mea, meb )
   en = TotEnergy( f, ja, jb, ua, ub, uab )

   call MPI_FILE_SET_VIEW(fh, disp, MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION,&
                           "native", MPI_INFO_NULL, ierr)

   if( rank == 0) then
      call MPI_FILE_WRITE(fh, ja, 1, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)
   endif

   call MPI_FILE_WRITE(fh, en, 1, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)

   ja = ja - dj
   jb = jb - dj
   disp = disp + 10

enddo

call MPI_FILE_CLOSE(fh, ierr)

call MPI_FINALIZE( ierr )

end program
