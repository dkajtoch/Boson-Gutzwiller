program two_mode_model3d_mpi

use mpi
use solvers
use measures
use parameters, only:dp

implicit none

! for MPI processes
integer ierr, comsize, rank, Num
integer fh
integer mean_id, mean_size, mean_rank, num_units, i, MPI_COMM_MEAN
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

! Open the file
call MPI_FILE_OPEN( MPI_COMM_WORLD, trim(file), MPI_MODE_WRONLY + MPI_MODE_CREATE, &
                    MPI_INFO_NULL, fh, ierr ) 

! -------------------------------------------------
! If less than 9 pocesses:
! 
!  Split only hoping loop among processes
! -------------------------------------------------
if( comsize < 9 ) then
   ! Total number of points
   Num = int(jmax/dj)

   call process_split( Num, comsize, num_units, disp, rank )
   ja = jmax - real(disp, dp) * dj
   jb = ja
   disp = disp * 10

   call MPI_FILE_SET_VIEW(fh, disp, MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION,&
                           "native", MPI_INFO_NULL, ierr)

   do i = 1, num_units

      call MPI_file_write( fh, ja, 1, MPI_double_precision, &
                           MPI_status_ignore, ierr )
      ! 0) Na, Nb
      mea = meaBase
      meb = mebBase
      call InitUniformNC( f, mea, meb )
      call GroundStateNC( f, ja, jb, ua, ub, uab, mea, meb )
      en = TotEnergy( f, ja, jb, ua, ub, uab )

      call MPI_file_write( fh, en, 1, MPI_double_precision, &
                           MPI_status_ignore, ierr )

      ! 1) Na - dN, Nb 
      mea = meaBase - dmean
      meb = mebBase
      call InitUniformNC( f, mea, meb )
      call GroundStateNC( f, ja, jb, ua, ub, uab, mea, meb )
      en = TotEnergy( f, ja, jb, ua, ub, uab )

      call MPI_file_write( fh, en, 1, MPI_double_precision, &
                           MPI_status_ignore, ierr )

      ! 2) Na + dN, Nb
      mea = meaBase + dmean
      meb = mebBase
      call InitUniformNC( f, mea, meb )
      call GroundStateNC( f, ja, jb, ua, ub, uab, mea, meb )
      en = TotEnergy( f, ja, jb, ua, ub, uab )

      call MPI_file_write( fh, en, 1, MPI_double_precision, &
                           MPI_status_ignore, ierr )

      ! 3) Na, Nb - dN
      mea = meaBase
      meb = mebBase - dmean
      call InitUniformNC( f, mea, meb )
      call GroundStateNC( f, ja, jb, ua, ub, uab, mea, meb )
      en = TotEnergy( f, ja, jb, ua, ub, uab )

      call MPI_file_write( fh, en, 1, MPI_double_precision, &
                           MPI_status_ignore, ierr )

      ! 4) Na, Nb + dN
      mea = meaBase
      meb = mebBase + dmean
      call InitUniformNC( f, mea, meb )
      call GroundStateNC( f, ja, jb, ua, ub, uab, mea, meb )
      en = TotEnergy( f, ja, jb, ua, ub, uab )

      call MPI_file_write( fh, en, 1, MPI_double_precision, &
                           MPI_status_ignore, ierr )

      ! 5) Na + dN, Nb - dN
      mea = meaBase + dmean
      meb = mebBase - dmean
      call InitUniformNC( f, mea, meb )
      call GroundStateNC( f, ja, jb, ua, ub, uab, mea, meb )
      en = TotEnergy( f, ja, jb, ua, ub, uab )

      call MPI_file_write( fh, en, 1, MPI_double_precision, &
                           MPI_status_ignore, ierr )

      ! 6) Na + dN, Nb + dN
      mea = meaBase + dmean
      meb = mebBase + dmean
      call InitUniformNC( f, mea, meb )
      call GroundStateNC( f, ja, jb, ua, ub, uab, mea, meb )
      en = TotEnergy( f, ja, jb, ua, ub, uab )

      call MPI_file_write( fh, en, 1, MPI_double_precision, &
                           MPI_status_ignore, ierr )

      ! 7) Na - dN, Nb - dN
      mea = meaBase - dmean
      meb = mebBase - dmean
      call InitUniformNC( f, mea, meb )
      call GroundStateNC( f, ja, jb, ua, ub, uab, mea, meb )
      en = TotEnergy( f, ja, jb, ua, ub, uab )

      call MPI_file_write( fh, en, 1, MPI_double_precision, &
                           MPI_status_ignore, ierr )

      ! 8) Na - dN, Nb + dN
      mea = meaBase - dmean
      meb = mebBase + dmean
      call InitUniformNC( f, mea, meb )
      call GroundStateNC( f, ja, jb, ua, ub, uab, mea, meb )
      en = TotEnergy( f, ja, jb, ua, ub, uab )

      call MPI_file_write( fh, en, 1, MPI_double_precision, &
                           MPI_status_ignore, ierr )

      ja = ja - dj
      jb = jb - dj

   enddo

! -------------------------------------------------
! More than 9 pocesses 
! -------------------------------------------------
else
   ! create group of processes
   mean_id = mod(rank, 9)
   call MPI_comm_split( MPI_COMM_WORLD, mean_id, rank, MPI_COMM_MEAN, ierr )

   ! Total number of processes per each energy computation
   call MPI_comm_size( MPI_COMM_MEAN, mean_size, ierr )
   call MPI_comm_rank( MPI_COMM_MEAN, mean_rank, ierr )

   ! Total number of points
   Num = int(jmax/dj)

   call process_split( Num, mean_size, num_units, disp, mean_rank )
   ja = jmax - real(disp, dp) * dj
   jb = ja
   disp = disp * 10

   if( mean_id /= 0 ) then
      disp = disp + (mean_id + 1)
   endif

   select case (mean_id)
      case( 0 )
         mea = meaBase
         meb = mebBase
      case( 1 )
         mea = meaBase- dmean
         meb = mebBase
      case( 2 )
         mea = meaBase + dmean
         meb = mebBase
      case( 3 )
         mea = meaBase
         meb = mebBase - dmean
      case( 4 )
         mea = meaBase
         meb = mebBase + dmean
      case( 5 ) 
         mea = meaBase + dmean
         meb = mebBase - dmean
      case( 6 )
         mea = meaBase + dmean
         meb = mebBase + dmean
      case( 7 )
         mea = meaBase - dmean
         meb = mebBase - dmean
      case( 8 )
         mea = meaBase - dmean
         meb = mebBase + dmean
   end select

   do i = 1, num_units

      ! redefine the file view
      call MPI_FILE_SET_VIEW(fh, disp, MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION,&
                              "native", MPI_INFO_NULL, ierr)
      if( mean_id == 0 ) then
         call MPI_file_write( fh, ja, 1, MPI_double_precision, &
                              MPI_status_ignore, ierr )
      endif

      call InitUniformNC( f, mea, meb )
      call GroundStateNC( f, ja, jb, ua, ub, uab, mea, meb )
      en = TotEnergy( f, ja, jb, ua, ub, uab )

      call MPI_file_write( fh, en, 1, MPI_double_precision, &
                           MPI_status_ignore, ierr )

      ja = ja - dj
      jb = jb - dj
      disp = disp + 10

   enddo

   call MPI_comm_free( MPI_COMM_MEAN, ierr )

endif

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
