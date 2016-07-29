program test

   use solvers
   use measures
   use parameters
   implicit none
   complex( kind = dp ), allocatable :: f(:,:)
   integer M, nmax
   real( kind = dp ) j, u
   real( kind = dp ) me, flu, en, mn
   real( kind = dp ) dj
   character(len=1024) filename
   complex( kind = dp ) ord
   integer i, na

   M = 10
   nmax = 8
   j = 0.0_dp
   u = 1.0_dp

   allocate( f(0:nmax,1:M) )

   !dj = 0.01_dp
   !filename='constmean_M=10_Na=8_Nb=8_Uab=2p0_Ua=1p0_Ub=0p95.dat' 

   !open(21, file=filename, action='write', status='replace')
   !close(21)

   mn = 0.2_dp*M

   do while( mn <= 1.0_dp*M )
      j = 0.0_dp

      do while( j <= 0.0_dp )

         call InitUniformNC( f, mn )
         call GroundStateNC( f, j, u, mn )

         me  = Mean( f, 1 )
         flu = Var( f, 1 )
         ord = Order( f, 1 )
         en  = TotEnergy( f, j, u )


         !open(21, file=filename, action='write', position='append')
         !   write(21,*)  ja, mna/M, mea, meb, meab, flua, flub, real(orda), &
         !          real(ordb), en1
         !close(21)


         j = j + dj

      enddo

      mn = mn + 0.05_dp*M

   enddo

end program
