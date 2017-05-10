program test01_c3d

   use solvers
   use measures
   use parameters
   use squeezing

   implicit none
   complex( kind = dp ), allocatable :: f(:,:,:,:,:)
   integer M, nmax
   integer i, j, k, na, nb
   real( kind = dp ) ja, jb, ua, ub, uab
   real( kind = dp ) mea, meb, meab, flucta, fluctb, en, sq
   real( kind = dp ) phase_en, phase_sq
   complex( kind = dp ) orda, ordb

   !random numbers
   integer :: time(8)
   integer :: seed(12)
   real( kind = dp ) :: mua, mub, nu

   ! error handling
   integer iostat_int
   character(256) iomsg_char

   ! import from namelist
   namelist /BHparams/ M, ja, jb, ua, ub, uab
   namelist /MEAN/ mea, meb

   open(138, file='BH_input.nml', iostat=iostat_int, iomsg=iomsg_char)
   if(iostat_int /= 0) then
      print *, "Open BH_input.dat failed with iostat = ", iostat_int, " iomsg = "//trim(iomsg_char)
      close(138)
   else
      read(138, BHparams)
      read(138, MEAN)
      close(138)
   endif

   ! redefine energies
   ub  = ub/ua
   uab = uab/ua
   ua  = 1.0_dp

   ! ---------------------------------------------------------------
   ! Compare Random and Uniform initial state for different nmax
   ! ---------------------------------------------------------------
   print *, '---------------------------------------------------------------'
   print *, '                        InitUniformNC                           '
   print *, '---------------------------------------------------------------'
   print *, 'M     = ', M
   print *, 'MeanA = ', mea 
   print *, 'MeanB = ', meb
   print *, 'ja    = ', ja
   print *, 'jb    = ', jb
   print *, 'ua    = ', ua
   print *, 'ub    = ', ub
   print *, 'uab   = ', uab

   do nmax = 2, 12

      allocate( f(0:nmax, 0:nmax, 1:M, 1:M, 1:M) )

      ! check random inital state preparation
      call InitUniformNC( f, mea, meb )

      call GroundStateNC( f, ja, jb, ua, ub, uab, mea, meb )

      ! on-site fluctuation
      flucta = VarA( f, 1, 1, 1 )
      fluctb = VarB( f, 1, 1, 1 )
      ! order parameter
      orda = OrderA( f, 1, 1, 1 )
      ordb = OrderB( f, 1, 1, 1 )
      ! total energy
      en = TotEnergy( f, ja, jb, ua, ub, uab )
      ! Squeezing parameter
      sq = SpinSqueezing( f )

      ! gauge transformation 
      call DATE_AND_TIME(values=time)     ! Get the current time 
      seed(1) = time(4) * (360000*time(5) + 6000*time(6) + 100*time(7) + time(8)) 
      CALL RANDOM_SEED(PUT=seed) 
      CALL RANDOM_NUMBER(HARVEST = mua)
      CALL RANDOM_NUMBER(HARVEST = mub)

      do k=1, ubound(f,5)
      do j=1, ubound(f,4)
      do i=1, ubound(f,3)
         CALL RANDOM_NUMBER(HARVEST = nu)
         do nb=0, ubound(f,2)
            do na=0, ubound(f,1)
               f(na,nb,i,j,k) = f(na,nb,i,j,k) * ZEXP( (0.0_dp, 1.0_dp)* (mua * real(na,dp) + mub * real(nb,dp) + nu ) )
            enddo
         enddo
      enddo
      enddo
      enddo

      ! total energy
      phase_en = TotEnergy( f, ja, jb, ua, ub, uab )
      ! Squeezing parameter
      phase_sq = SpinSqueezing( f )

      print *, '---------------------------------------------------------------'
      print *, 'nmax           = ', nmax
      print *, '---------------------------------------------------------------'
      print *, 'VarA           = ', flucta
      print *, 'VarB           = ', fluctb
      print *, 'OrdA           = ', orda
      print *, 'OrdB           = ', ordb
      print *, 'Energy         = ', en
      print *, 'PhaseEnergy    = ', phase_en
      print *, 'Squeezing      = ', sq
      print *, 'PhaseSqueezing = ', phase_sq

      deallocate(f)

   enddo

   print *, '---------------------------------------------------------------'
   print *, '                        InitRandomNC                          '

   do nmax = 2,12

      allocate( f(0:nmax, 0:nmax, 1:M, 1:M, 1:M) )

      ! check random inital state preparation
      call InitRandomNC( f, mea, meb )

      call GroundStateNC( f, ja, jb, ua, ub, uab, mea, meb )

      ! on-site fluctuation
      flucta = VarA( f, 1, 1, 1 )
      fluctb = VarB( f, 1, 1 ,1 )
      ! order parameter
      orda = OrderA( f, 1, 1 ,1 )
      ordb = OrderB( f, 1, 1 ,1 )
      ! total energy
      en = TotEnergy( f, ja, jb, ua, ub, uab )
      ! Squeezing parameter
      sq = SpinSqueezing( f )

      ! gauge transformation 
      call DATE_AND_TIME(values=time)     ! Get the current time 
      seed(1) = time(4) * (360000*time(5) + 6000*time(6) + 100*time(7) + time(8)) 
      CALL RANDOM_SEED(PUT=seed) 
      CALL RANDOM_NUMBER(HARVEST = mua)
      CALL RANDOM_NUMBER(HARVEST = mub)

      do k=1, ubound(f,5)
      do j=1, ubound(f,4)
      do i=1, ubound(f,3)
         CALL RANDOM_NUMBER(HARVEST = nu)
         do nb=0, ubound(f,2)
            do na=0, ubound(f,1)
               f(na,nb,i,j,k) = f(na,nb,i,j,k) * ZEXP( (0.0_dp, 1.0_dp)* (mua * real(na,dp) + mub * real(nb,dp) + nu ) )
            enddo
         enddo
      enddo
      enddo
      enddo

      ! total energy
      phase_en = TotEnergy( f, ja, jb, ua, ub, uab )
      ! Squeezing parameter
      phase_sq = SpinSqueezing( f )

      print *, '---------------------------------------------------------------'
      print *, 'nmax           = ', nmax
      print *, '---------------------------------------------------------------'
      print *, 'VarA           = ', flucta
      print *, 'VarB           = ', fluctb
      print *, 'OrdA           = ', orda
      print *, 'OrdB           = ', ordb
      print *, 'Energy         = ', en
      print *, 'PhaseEnergy    = ', phase_en
      print *, 'Squeezing      = ', sq
      print *, 'PhaseSqueezing = ', phase_sq

      deallocate(f)

   enddo

end program
