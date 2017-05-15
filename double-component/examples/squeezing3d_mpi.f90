program squeezing3d_mpi

use mpi
use solvers
use measures
use parameters, only:dp

implicit none

! =============================================================
! L     - box size in one direction
! nmax  - Fock state basis cutoff
! Num   - total number of states: N + 1
! ja    - hoping for component A
! jb    - hoping for component B
! dj    - step size for hoping
! jmax  - maximal hoping (starting value for 'do' loop)
! ua    - scattering length a_11
! ub    - scattering length a_22
! uab   - scattering length a_12
! TotEn - Total Gutzwiller energy
! dmean - step size for the mean
! ordA  - on-site order parameter of component A
! ordB  - on-site order parameter of component B
! flucA - on-site atom number fluctuations in component A
! flucB - on-site atom number fluctuations in component B
! corr2 - on-site g^2 correlator
! corr3 - on-site g^3 correlator
! mea   - total mean number of atoms in component A
! meb   - total mean number of atoms in component B
! =============================================================

integer L, nmax, Num
real( kind = dp ) ja, jb, dj, jmax, ua, ub, uab 
real( kind = dp ) TotEn, dmean, flucA, flucB, mea, meb
complex( kind = dp ) ordA, ordB
real( kind = dp ) corr2(3), corr3(4)

! =============================================================
! f  - Gutzwiller vector 
! =============================================================

complex( kind = dp ), allocatable :: f(:,:,:,:,:)

! =============================================================
! MPI communication variables:
!
!    ierr      - error feedback
!    comsize   - total number of processes in MPI_COMM_WORLD
!    rank      - label of the proces
!    avg       - average number of units per proces
!    rem       - remaining number of units
!    num_units - number of units per proces
!    disp      - displacement from the begining of a file
!    fh...     - file handlers
!    k         - counter for the 'do' loop
! =============================================================

integer ierr, comsize, rank, avg, mean_int, num_units, rem, k, i
integer fhTotalEnergy, fhOrderA, fhOrderB, fhFlucA, fhFlucB, fhCorr2, fhCorr3
integer(kind=MPI_OFFSET_KIND) disp

! =============================================================
! Filenames:
! 
!    fileTotalEnergy - Gutzwiller total energy
!    fileOrderA      - order parameter A
!    fileOrderB      - order parameter B
!    fileFlucA       - fluctuations in A
!    fileFlucB       - fluctuations in B
!    fileCorr2       - correlator g^2
!    fileCorr3       - correlator g^3
! =============================================================

character(len=100) fileTotalEnergy
character(len=100) fileOrderA
character(len=100) fileOrderB
character(len=100) fileFlucA
character(len=100) fileFlucB
character(len=100) fileCorr2
character(len=100) fileCorr3

! =============================================================
! Namelists
! =============================================================

namelist /BHparams/ L, nmax, ua, ub, uab
namelist /SQparams/ dmean, jmax, dj
namelist /FLparams/ fileTotalEnergy, fileOrderA, fileOrderB, &
                    fileFlucA, fileFlucB, fileCorr2, fileCorr3

open(138, file='BH_inp.nml')
read(138, BHparams)
read(138, SQparams)
read(138, FLparams)
close(138)

! allocate memory
allocate( f(0:nmax,0:nmax,1:L,1:L,1:L) )

! redefine interaction energies
ub  = ub/ua
uab = uab/ua
ua  = ua/ua

! total number of states
Num = L**3 + 1

call MPI_INIT( ierr ) 
call MPI_COMM_SIZE( MPI_COMM_WORLD, comsize, ierr )
call MPI_COMM_RANK( MPI_COMM_WORLD, rank, ierr )

! average number of units per process
avg = Num/comsize
! remaining number 
rem = mod(Num,comsize)

! split squeezing into all processes
if( rank < rem ) then
   num_units = avg + 1

   ! intial mean values
   mean_int = num_units*rank
   mea = real( mean_int, dp )
   meb = real(L**3, dp) - mea

else
   num_units = avg

   mean_int = (avg+1)*rem + num_units*(rank-rem)
   mea = real( mean_int, dp )
   meb = real(L**3, dp) - mea

endif

! Open All the files
call MPI_FILE_OPEN( MPI_COMM_WORLD, trim(fileTotalEnergy), MPI_MODE_WRONLY + MPI_MODE_CREATE, &
                    MPI_INFO_NULL, fhTotalEnergy, ierr )
call MPI_FILE_OPEN( MPI_COMM_WORLD, trim(fileOrderA), MPI_MODE_WRONLY + MPI_MODE_CREATE, &
                    MPI_INFO_NULL, fhOrderA, ierr )
call MPI_FILE_OPEN( MPI_COMM_WORLD, trim(fileOrderB), MPI_MODE_WRONLY + MPI_MODE_CREATE, &
                    MPI_INFO_NULL, fhOrderB, ierr )
call MPI_FILE_OPEN( MPI_COMM_WORLD, trim(fileFlucA), MPI_MODE_WRONLY + MPI_MODE_CREATE, &
                    MPI_INFO_NULL, fhFlucA, ierr )
call MPI_FILE_OPEN( MPI_COMM_WORLD, trim(fileFlucB), MPI_MODE_WRONLY + MPI_MODE_CREATE, &
                    MPI_INFO_NULL, fhFlucB, ierr )
call MPI_FILE_OPEN( MPI_COMM_WORLD, trim(fileCorr2), MPI_MODE_WRONLY + MPI_MODE_CREATE, &
                    MPI_INFO_NULL, fhCorr2, ierr )
call MPI_FILE_OPEN( MPI_COMM_WORLD, trim(fileCorr3), MPI_MODE_WRONLY + MPI_MODE_CREATE, &
                    MPI_INFO_NULL, fhCorr3, ierr )

do i = 1, num_units

   ja = jmax
   jb = jmax

   k = 0 
   do while( ja > 0.0_dp )

      call InitUniformNC( f, mea, meb )
      call GroundStateNC( f, ja, jb, ua, ub, uab, mea, meb )

      ! ==============
      ! Total energy
      ! ==============
      TotEn = TotEnergy( f, ja, jb, ua, ub, uab )
      ! change the file view
      if( mean_int == 0 ) then
         disp = k * (Num + 1) * SIZEOF(TotEn)
         call MPI_FILE_SET_VIEW( fhTotalEnergy, disp, MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION, &
                                 "native", MPI_INFO_NULL, ierr )
         call MPI_FILE_WRITE( fhTotalEnergy, ja, 1, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr )
         call MPI_FILE_WRITE( fhTotalEnergy, TotEn, 1, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr )
      else
         disp = (k * (Num + 1) + mean_int + 1) * SIZEOF(TotEn)
         call MPI_FILE_SET_VIEW( fhTotalEnergy, disp, MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION, &
                                 "native", MPI_INFO_NULL, ierr )
         call MPI_FILE_WRITE( fhTotalEnergy, TotEn, 1, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr )
      endif

      ! ==================
      ! Order parameter A
      ! ==================
      ordA  = OrderA(f, 1, 1, 1)
      ! change the file view
      if( mean_int == 0 ) then
         disp = k * (2*Num + 1) * SIZEOF(TotEn)
         call MPI_FILE_SET_VIEW( fhOrderA, disp, MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION, &
                                 "native", MPI_INFO_NULL, ierr )
         call MPI_FILE_WRITE( fhOrderA, ja, 1, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr )
         call MPI_FILE_WRITE( fhOrderA, real(ordA), 1, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr )
         call MPI_FILE_WRITE( fhOrderA, aimag(ordA), 1, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr )
      else
         disp = (k * (2*Num + 1) + 2*mean_int + 1) * SIZEOF(TotEn)
         call MPI_FILE_SET_VIEW( fhOrderA, disp, MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION, &
                                 "native", MPI_INFO_NULL, ierr )
         call MPI_FILE_WRITE( fhOrderA, real(ordA), 1, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr )
         call MPI_FILE_WRITE( fhOrderA, aimag(ordA), 1, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr )
      endif

      ! ==================
      ! Order parameter B
      ! ==================
      ordB  = OrderB(f, 1, 1, 1) 
      ! change the file view
      if( mean_int == 0 ) then
         disp = k * (2*Num + 1) * SIZEOF(TotEn)
         call MPI_FILE_SET_VIEW( fhOrderB, disp, MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION, &
                                 "native", MPI_INFO_NULL, ierr )
         call MPI_FILE_WRITE( fhOrderB, ja, 1, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr )
         call MPI_FILE_WRITE( fhOrderB, real(ordA), 1, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr )
         call MPI_FILE_WRITE( fhOrderB, aimag(ordA), 1, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr )
      else
         disp = (k * (2*Num + 1) + 2*mean_int + 1) * SIZEOF(TotEn)
         call MPI_FILE_SET_VIEW( fhOrderB, disp, MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION, &
                                 "native", MPI_INFO_NULL, ierr )
         call MPI_FILE_WRITE( fhOrderB, real(ordA), 1, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr )
         call MPI_FILE_WRITE( fhOrderB, aimag(ordA), 1, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr )
      endif
      
      ! ============================
      ! atom number fluctuation in A
      ! ============================
      flucA = VarA(f, 1, 1, 1)
      ! change the file view
      if( mean_int == 0 ) then
         disp = k * (Num + 1) * SIZEOF(TotEn)
         call MPI_FILE_SET_VIEW( fhFlucA, disp, MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION, &
                                 "native", MPI_INFO_NULL, ierr )
         call MPI_FILE_WRITE( fhFlucA, ja, 1, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr )
         call MPI_FILE_WRITE( fhFlucA, flucA, 1, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr )
      else
         disp = (k * (Num + 1) + mean_int + 1) * SIZEOF(TotEn)
         call MPI_FILE_SET_VIEW( fhFlucA, disp, MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION, &
                                 "native", MPI_INFO_NULL, ierr )
         call MPI_FILE_WRITE( fhFlucA, flucA, 1, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr )
      endif

      ! ============================
      ! atom number fluctuation in B
      ! ============================
      flucB = VarB(f, 1, 1, 1)
      ! change the file view
      if( mean_int == 0 ) then
         disp = k * (Num + 1) * SIZEOF(TotEn)
         call MPI_FILE_SET_VIEW( fhFlucB, disp, MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION, &
                                 "native", MPI_INFO_NULL, ierr )
         call MPI_FILE_WRITE( fhFlucB, ja, 1, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr )
         call MPI_FILE_WRITE( fhFlucB, flucB, 1, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr )
      else
         disp = (k * (Num + 1) + mean_int + 1) * SIZEOF(TotEn)
         call MPI_FILE_SET_VIEW( fhFlucB, disp, MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION, &
                                 "native", MPI_INFO_NULL, ierr )
         call MPI_FILE_WRITE( fhFlucB, flucB, 1, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr )
      endif
      
      ! ===============
      ! Correlator g^2
      ! ===============
      corr2 = G2correlator(f, 1, 1, 1)
      ! change the file view
      if( mean_int == 0 ) then
         disp = k * (3*Num + 1) * SIZEOF(TotEn)
         call MPI_FILE_SET_VIEW( fhCorr2, disp, MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION, &
                                 "native", MPI_INFO_NULL, ierr )
         call MPI_FILE_WRITE( fhCorr2, ja, 1, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr )
         call MPI_FILE_WRITE( fhCorr2, corr2, 3, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr )
      else
         disp = (k * (3*Num + 1) + 3*mean_int + 1) * SIZEOF(TotEn)
         call MPI_FILE_SET_VIEW( fhCorr2, disp, MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION, &
                                 "native", MPI_INFO_NULL, ierr )
         call MPI_FILE_WRITE( fhCorr2, corr2, 3, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr )
      endif

      ! ===============
      ! Correlator g^3
      ! ===============
      corr3 = G3correlator(f, 1, 1, 1)
      ! change the file view
      if( mean_int == 0 ) then
         disp = k * (4*Num + 1) * SIZEOF(TotEn)
         call MPI_FILE_SET_VIEW( fhCorr3, disp, MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION, &
                                 "native", MPI_INFO_NULL, ierr )
         call MPI_FILE_WRITE( fhCorr3, ja, 1, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr )
         call MPI_FILE_WRITE( fhCorr3, corr3, 4, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr )
      else
         disp = (k * (4*Num + 1) + 4*mean_int + 1) * SIZEOF(TotEn)
         call MPI_FILE_SET_VIEW( fhCorr3, disp, MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION, &
                                 "native", MPI_INFO_NULL, ierr )
         call MPI_FILE_WRITE( fhCorr3, corr3, 4, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr )
      endif

      ja = ja - dj
      jb = jb - dj
      k  = k + 1

   enddo

   mean_int = mean_int + 1
   mea = mea + dmean
   meb = meb - dmean

enddo

call MPI_FILE_CLOSE( fhTotalEnergy, ierr )
call MPI_FILE_CLOSE( fhOrderA, ierr )
call MPI_FILE_CLOSE( fhOrderB, ierr )
call MPI_FILE_CLOSE( fhFlucA, ierr )
call MPI_FILE_CLOSE( fhFlucB, ierr )
call MPI_FILE_CLOSE( fhCorr2, ierr )
call MPI_FILE_CLOSE( fhCorr3, ierr )

call MPI_FINALIZE( ierr )

end program
