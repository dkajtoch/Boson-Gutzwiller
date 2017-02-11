program test02_c1d

use mpi
use solvers
use measures
use squeezing
use parameters, only:dp

implicit none

! for the MPI processes
integer ierr, comsize, rank, i
integer status(MPI_STATUS_SIZE) 

integer, parameter :: nmax = 8
integer M, Num, rem, avg, num_units
real( kind = dp ) ja, jb, ua, ub, uab, en, mea, meb
real( kind = dp ) dmean, jmax, dj
real( kind = dp ) sq
complex( kind = dp ), allocatable :: f(:,:,:)
complex( kind = dp ), allocatable :: master_f(:,:,:)
namelist /BHparams/ M, ja, jb, ua, ub, uab

open(138, file='BH_input.nml')
read(138, BHparams)
close(138)

allocate( f(0:nmax,0:nmax,1:M) )

! redefine on-site interaction
uab = uab/ua
ub  = ub/ua
ua  = ua/ua

! total number of states
Num = M + 1

call MPI_INIT( ierr ) 
call MPI_COMM_SIZE( MPI_COMM_WORLD, comsize, ierr )
call MPI_COMM_RANK( MPI_COMM_WORLD, rank, ierr )

! Master process which collects the data
if( rank == 0 ) then
   print *, "Hello from master process!"
   print *, "Process rank: ", rank
   print *, "Total number of processes: ", comsize

   ! allocate memory
   allocate( master_f(0:nmax,0:nmax,1:M) )

   ! flush array
   f(:,:,:) = (0.0_dp,0.0_dp)

   ! receive data from slave processes 
   do i=1, Num
      call MPI_RECV( master_f, size(master_f), MPI_DOUBLE_COMPLEX, MPI_ANY_SOURCE, &
                     MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)

      f(:,:,:) = f(:,:,:) + master_f(:,:,:)
   enddo

   ! normalize and calculate squeezing
   call normalize(f)

   sq = SpinSqueezing(f)

   print *, "Not normalized squeezing: ", sq
   print *, "Normalized squeezing: ", sq*M

else

   ! decrease total number of processes by 1
   comsize = comsize - 1

   ! average number of units per process
   avg = Num/comsize
   ! remaining number 
   rem = mod(Num,comsize)

   ! split squeezing into all processes
   if( rank-1 < rem ) then
      num_units = avg + 1

      ! intial mean values
      mea = real( num_units*(rank-1), dp )
      meb = real(M, dp) - mea

   else
      num_units = avg

      mea = real( (avg+1)*rem + num_units*(rank-1-rem), dp )
      meb = real(M, dp) - mea

   endif

   do i = 1, num_units

      call InitUniformNC( f, mea, meb )
      call GroundStateNC( f, ja, jb, ua, ub, uab, mea, meb )

      mea = mea + dmean
      meb = meb - dmean

      call MPI_SEND( f, size(f), MPI_DOUBLE_COMPLEX, 0, 2017, MPI_COMM_WORLD, ierr)

   enddo

endif

call MPI_FINALIZE( ierr )

end program
