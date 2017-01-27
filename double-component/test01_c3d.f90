program test01_c3d

   use solvers
   use measures
   use parameters

   implicit none
   complex( kind = dp ), allocatable :: f(:,:,:,:,:)
   integer M, nmax
   integer i, na, nb
   real( kind = dp ) ja, jb, ua, ub, uab
   real( kind = dp ) mea, meb, meab, flucta, fluctb, en
   real( kind = dp ) normf, meanaf, meanbf, meanabf, tmp
   complex( kind = dp ) orda, ordb

   ! error handling
   integer iostat_int
   character(256) iomsg_char

   ! import from namelist
   namelist /BHparams/ M, ja, jb, ua, ub, uab, mea, meb

   open(138, file='BH_input.nml', iostat=iostat_int, iomsg=iomsg_char)
   if(iostat_int /= 0) then
      print *, "Open BH_input.dat failed with iostat = ", iostat_int, " iomsg = "//trim(iomsg_char)
      close(138)
   else
      read(138, BHparams)
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


      print *, '---------------------------------------------------------------'
      print *, 'nmax   = ', nmax
      print *, '---------------------------------------------------------------'
      print *, 'VarA   = ', flucta
      print *, 'VarB   = ', fluctb
      print *, 'OrdA   = ', orda
      print *, 'OrdB   = ', ordb
      print *, 'Energy = ', en

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


      print *, '---------------------------------------------------------------'
      print *, 'nmax   = ', nmax
      print *, '---------------------------------------------------------------'
      print *, 'VarA   = ', flucta
      print *, 'VarB   = ', fluctb
      print *, 'OrdA   = ', orda
      print *, 'OrdB   = ', ordb
      print *, 'Energy = ', en

      deallocate(f)

   enddo

end program
