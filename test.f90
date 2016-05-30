program test

   use solvers
   use measures
   implicit none
   complex( kind = dp ), allocatable :: f(:,:,:)
   integer M, nmax
   real( kind = dp ) ja, jb, ua, ub, uab
   real( kind = dp ) mea, meb, meab, flua, flub, en1, mna, mnb
   real( kind = dp ) dj
   character(len=1024) filename
   complex( kind = dp ) orda, ordb
   integer i ,na, nb

   M = 10
   nmax = 8
   ja = 0.0_dp
   jb = 0.0_dp
   ua = 1.0_dp
   ub = 0.95_dp
   uab = 2.0_dp

   allocate( f(1:M,0:nmax,0:nmax) )

   dj = 0.01_dp
   filename='constmean_M=10_Na=8_Nb=8_Uab=2p0_Ua=1p0_Ub=0p95.dat' 

   open(21, file=filename, action='write', status='replace')
   close(21)

   mna = 0.2_dp*M
   mnb = 0.8_dp*M

   do while( mna <= 1.0_dp*M )
      ja = 0.0_dp
      jb = 0.0_dp

      do while( ja <= 0.0_dp )

         call InitUniformNC( f, mna,  mnb )
         call GroundStateNC( f, ja, jb, ua, ub, uab, mna, mnb )

         mea = MeanA( f, 1 )
         meb = MeanB( f, 1 )
         meab = MeanAB( f, 1 )
         flua = VarA( f, 1 )
         flub = VarB( f, 1 )
         orda = OrderA( f, 1 )
         ordb = OrderB( f, 1 )
         en1 = TotEnergy( f, ja, jb, ua, ub, uab )


         open(21, file=filename, action='write', position='append')
            write(21,*)  ja, mna/M, mea, meb, meab, flua, flub, real(orda), &
                   real(ordb), en1
         close(21)


         ja = ja + dj
         jb = jb + dj

      enddo

      mna = mna + 0.05_dp*M
      mnb = mnb - 0.05_dp*M

   enddo

end program
