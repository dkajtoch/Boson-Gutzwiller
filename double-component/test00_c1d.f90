program test_c1d

   use solvers
   use measures
   use parameters

   implicit none
   complex( kind = dp ), allocatable :: f(:,:,:)
   integer M, nmax
   integer i, na, nb
   real( kind = dp ) ja, jb, ua, ub, uab
   real( kind = dp ) mea, meb, meab, flucta, fluctb, en
   real( kind = dp ) normf, meanaf, meanbf, meanabf, tmp
   complex( kind = dp ) orda, ordb

   M    = 10
   nmax = 8
   ja   = 0.3_dp
   jb   = 0.3_dp
   ua   = 1.0_dp
   ub   = 95.47_dp/100.44_dp
   uab  = 80.0_dp/100.44_dp

   allocate( f(0:nmax, 0:nmax, 1:M) )

   ! check random inital state preparation
   mea = 0.0_dp * M
   meb = 1.0_dp * M
   call InitUniformNC( f, mea, meb )

   normf = norm( f, 1 )
   meanaf = 0.0_dp
   meanbf = 0.0_dp
   do i = 1, ubound(f,3)
      tmp = MeanA( f, i )
      meanaf = meanaf + tmp
      tmp = MeanB( f, i )
      meanbf = meanbf + tmp
   enddo

   print *, '---------------------------------------------------------------'
   print *, '                        InitRandomNC                           '
   print *, '---------------------------------------------------------------'
   print *, 'MeanA = ', mea, 'MeanB = ', meb
   print *, 'After randomization:'
   print *, 'MeanA = ', meanaf, 'MeanB = ', meanbf
   print *, 'Norm  = ', normf
   print *, 'f(2,3,4) = ', f(2,3,4)
   print *, 'f(2,3,5) = ', f(2,3,5)

   ! after random vector preparation
   call GroundStateNC( f, ja, jb, ua, ub, uab, mea, meb )
   ! mean value
   meanaf  = 0.0_dp
   meanbf  = 0.0_dp
   meanabf = 0.0_dp
   do i = 1, ubound(f,3)
      tmp = MeanA( f, i )
      meanaf = meanaf + tmp
      tmp = MeanB( f, i )
      meanbf = meanbf + tmp
      tmp = MeanAB( f, i )
      meanabf = meanabf + tmp
   enddo
   ! on-site fluctuation
   flucta = VarA( f, 1 )
   fluctb = VarB( f, 1 )
   ! order parameter
   orda = OrderA( f, 1 )
   ordb = OrderB( f, 1 )
   ! total energy
   en = TotEnergy( f, ja, jb, ua, ub, uab )

   print *, '---------------------------------------------------------------'
   print *, '                        GroundStateNC                          '
   print *, '---------------------------------------------------------------'
   print *, 'ja  = ', ja
   print *, 'jb  = ', jb
   print *, 'ua  = ', ua
   print *, 'ub  = ', ub
   print *, 'uab = ', uab
   print *, '---------------------------------------------------------------'
   print *, 'MeanA  = ', meanaf
   print *, 'MeanB  = ', meanbf 
   print *, 'MeanAB = ', meanabf
   print *, 'VarA   = ', flucta
   print *, 'VarB   = ', fluctb
   print *, 'OrdA   = ', orda
   print *, 'OrdB   = ', ordb
   print *, 'Energy = ', en

end program
