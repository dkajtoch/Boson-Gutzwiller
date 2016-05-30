module solvers

   use parameters
   use measures, only: OrderA, OrderB, normalize
   use fixed_mean

   private::solve
   interface solve
      module procedure c1d_solve
   end interface

   interface GroundState
      module procedure c1d_GroundState
   end interface

   interface GroundStateNC
      module procedure c1d_GroundStateNC
   end interface

   interface InitUniform
      module procedure c1d_InitUniform
   end interface

   interface InitUniformNC
      module procedure c1d_InitUniformNC
   end interface

   private::sigma_norm
   interface sigma_norm
      module procedure c1d_sigma_norm
   end interface

   interface const_mean
      module procedure c1d_const_mean
   end interface

   contains

!  | ---------------------------------------------- |
!  | Update Gutzwiller coeffcients using right-hand | 
!  | side of Gutzwiller equation                    |
!  | ---------------------------------------------- |
   subroutine c1d_solve( f, k, ja, jb, ua, ub, uab, mua, mub, mode )

      implicit none
      complex( kind = dp), intent(inout), allocatable :: f(:,:,:)
      complex( kind = dp ), intent(inout), allocatable :: k(:,:,:)
      real( kind = dp ), intent(in) :: ja, jb, ua, ub, uab, mua, mub
      character(len=2), intent(in) :: mode

      ! local variables
      integer i, na, nb
      complex( kind = dp ), allocatable :: ArrayOrderA(:), ArrayOrderB(:)
      complex( kind = dp ) :: ordA, ordB, mlt

      if (mode .eq. 'EV') then
         mlt = im
      else if (mode .eq. 'GS') then
         mlt = re
      else
         mlt = (0.0_dp, 0.0_dp) 
         write(*,*)  "Wrong mode in SOLVE subroutine!"
      endif

      allocate( ArrayOrderA( size(f,1) ), ArrayOrderB( size(f,1) ) )

      ! collect order parameters
      do i = 1, ubound(f,1)
         
         ArrayOrderA(i) = 0.0_dp*re
         ArrayOrderB(i) = 0.0_dp*re

         do na = 0, ubound(f,2)
            do nb = 0, ubound(f,3)
               
               if( na > 0 ) then
                  ArrayOrderA(i) = ArrayOrderA(i) + SQRT( real(na,dp) ) * conjg( f(i,na-1,nb) ) * f(i,na,nb)
               endif
               if( nb > 0 ) then
                  ArrayOrderB(i) = ArrayOrderB(i) + SQRT( real(nb,dp) ) * conjg( f(i,na,nb-1) ) * f(i,na,nb)
               endif

            enddo
         enddo

      enddo

      ! solve coupled system of equations
      do i = 1, ubound(f,1)

         ! calculate nearest-neighbour order parameter (periodic boundary conditions)
         if( i .eq. 1 ) then
            ordA = ArrayOrderA( ubound(f,1) ) + ArrayOrderA(i+1)
            ordB = ArrayOrderB( ubound(f,1) ) + ArrayOrderB(i+1)
         elseif( i .eq. ubound(f,1) ) then
            ordA = ArrayOrderA(1) + ArrayOrderA(i-1)
            ordB = ArrayOrderB(1) + ArrayOrderB(i-1)
         else
            ordA = ArrayOrderA(i+1) + ArrayOrderA(i-1)
            ordB = ArrayOrderB(i+1) + ArrayOrderB(i-1)
         endif
         
         do na = 0, ubound(f,2)
            do nb = 0, ubound(f,3)

               ! chop
               if( abs( real( f(i,na,nb) ) ) < chopCutoff ) then
                  f(i,na,nb) = ( f(i,na,nb) - conjg(f(i,na,nb)) )/2.0_dp
               endif
               if( abs( aimag( f(i,na,nb) ) ) < chopCutoff ) then
                  f(i,na,nb) = ( f(i,na,nb) + conjg(f(i,na,nb)) )/2.0_dp
               endif


               k(i,na,nb) = mlt * ( real(na,dp) * ( mua - ua/2.0_dp * (real(na,dp) - 1.0_dp) ) * f(i,na,nb) + &
                                    real(nB,dp) * ( mub - ub/2.0_dp * (real(nb,dp) - 1.0_dp) ) * f(i,nA,nB) - &
                                    uab * real(na,dp) * real(nb,dp) * f(i,nA,nB) )
               if( na .lt. ubound(f,2) ) then
                  k(i,na,nb) = k(i,na,nb) + mlt * ja * ordA * SQRT( real(na+1,dp) ) * f(i,na+1,nb)
               endif
               if( na .gt. 0 ) then
                  k(i,na,nb) = k(i,na,nb) + mlt * ja * ordA * SQRT( real(na,dp) ) * f(i,na-1,nb)
               endif

               if( nb .lt. ubound(f,3) ) then
                  k(i,na,nb) = k(i,na,nb) + mlt * jb * ordB * SQRT( real(nb+1,dp) ) * f(i,na,nb+1)
               endif
               if( nb .gt. 0 ) then
                  k(i,na,nb) = k(i,na,nb) + mlt * jb * ordB * SQRT( real(nb,dp) ) * f(i,na,nb-1)
               endif

            enddo
         enddo

      enddo

   end subroutine

!  | ---------------------------------------------- |
!  | Propagate in imaginary time to find the ground | 
!  | state configuration                            |
!  | ---------------------------------------------- |
   subroutine c1d_GroundStateNC( f, ja, jb, ua, ub, uab, mna, mnb )

      implicit none
      complex( kind = dp ), intent(inout), allocatable :: f(:,:,:)
      real( kind = dp ), intent(in) :: ja, jb, ua, ub, uab, mna, mnb

      ! | -------------------- |
      ! | For Runge-Kutta step |
      ! | -------------------- |
      complex( kind = dp ), allocatable :: k1(:,:,:)
      complex( kind = dp ), allocatable :: k2(:,:,:)
      complex( kind = dp ), allocatable :: k3(:,:,:)
      complex( kind = dp ), allocatable :: k4(:,:,:)
      complex( kind = dp ), allocatable :: tmp(:,:,:)
      complex( kind = dp ), allocatable :: fpom(:,:,:)
      
      real( kind = dp ) :: error, ierr, mua, mub
      real( kind = dp ) :: tstart, tstop
      integer i, na, nb, cnt, info

      allocate( k1(lbound(f,1):ubound(f,1), lbound(f,2):ubound(f,2), lbound(f,3):ubound(f,3)) )
      allocate( k2(lbound(f,1):ubound(f,1), lbound(f,2):ubound(f,2), lbound(f,3):ubound(f,3)) )
      allocate( k3(lbound(f,1):ubound(f,1), lbound(f,2):ubound(f,2), lbound(f,3):ubound(f,3)) )
      allocate( k4(lbound(f,1):ubound(f,1), lbound(f,2):ubound(f,2), lbound(f,3):ubound(f,3)) )
      allocate( tmp(lbound(f,1):ubound(f,1), lbound(f,2):ubound(f,2), lbound(f,3):ubound(f,3)) )
      allocate( fpom(lbound(f,1):ubound(f,1), lbound(f,2):ubound(f,2), lbound(f,3):ubound(f,3)) )

      fpom(:,:,:) = f(:,:,:)
      cnt = 0
      mua = 0.0_dp; mub = 0.0_dp

      do while( .true. )

         cnt = cnt + 1 
 
         ! 1st Runge-Kutta step
         call solve( f, k1, ja, jb, ua, ub, uab, mua, mub, 'GS' )
         ! update vector to the 2nd step
         tmp(:,:,:) = f(:,:,:) + dt/2.0_dp * k1(:,:,:)
         ! 2nd Runge-Kutta step
         call solve( tmp, k2, ja, jb, ua, ub, uab, mua, mub, 'GS' )
         ! update vector for the 3rd step
         tmp(:,:,:) = f(:,:,:) + dt/2.0_dp * k2(:,:,:)
         ! 3rd Runge-Kutta step
         call solve( tmp, k3, ja, jb, ua, ub, uab, mua, mub, 'GS' )
         ! update vector for the 4th step
         tmp(:,:,:) = f(:,:,:) + dt * k3(:,:,:)
         ! 4th Runge-Kutta step
         call solve( tmp, k4, ja, jb, ua, ua, uab, mua, mub, 'GS' )
         ! update vector
         f(:,:,:) = f(:,:,:) + dt/6.0_dp * ( k1(:,:,:) + 2.0_dp *&
                    k2(:,:,:) + 2.0_dp * k3(:,:,:) + k4(:,:,:) )

         ! normalize vector
         call const_mean( f, mna, mnb, info )
         
         if( info /= 0 ) then
            print *, 'info = ', info
            exit
         endif

         ! calculate error
         if( modulo(cnt, stepsForJudge) .eq. 0 ) then

            print *, cnt
            error = 0.0_dp
            do i = 1, ubound(f,1)
               do na = 0, ubound(f,2)
                  do nb = 0, ubound(f,3)
                     ierr = abs( ( f(i,nA,nB) - fpom(i,nA,nB) ) )
                     if (ierr .gt. error) then
                        error = ierr
                     endif
                  enddo
               enddo
            enddo

            if( error < convCriterion ) then
               exit
            endif

            print *, 'error = ', error
            fpom(:,:,:) = f(:,:,:)

         endif

      end do

   end subroutine

!  | ---------------------------------------------- |
!  | Propagate in imaginary time to find the ground | 
!  | state configuration                            |
!  | ---------------------------------------------- |
   subroutine c1d_GroundState( f, ja, jb, ua, ub, uab, mua, mub )

      implicit none
      complex( kind = dp ), intent(inout), allocatable :: f(:,:,:)
      real( kind = dp ), intent(in) :: ja, jb, ua, ub, uab, mua, mub

      ! | -------------------- |
      ! | For Runge-Kutta step |
      ! | -------------------- |
      complex( kind = dp ), allocatable :: k1(:,:,:)
      complex( kind = dp ), allocatable :: k2(:,:,:)
      complex( kind = dp ), allocatable :: k3(:,:,:)
      complex( kind = dp ), allocatable :: k4(:,:,:)
      complex( kind = dp ), allocatable :: tmp(:,:,:)
      complex( kind = dp ), allocatable :: fpom(:,:,:)
      
      real( kind = dp ) :: error, ierr
      integer i, na, nb, cnt

      allocate( k1(lbound(f,1):ubound(f,1), lbound(f,2):ubound(f,2), lbound(f,3):ubound(f,3)) )
      allocate( k2(lbound(f,1):ubound(f,1), lbound(f,2):ubound(f,2), lbound(f,3):ubound(f,3)) )
      allocate( k3(lbound(f,1):ubound(f,1), lbound(f,2):ubound(f,2), lbound(f,3):ubound(f,3)) )
      allocate( k4(lbound(f,1):ubound(f,1), lbound(f,2):ubound(f,2), lbound(f,3):ubound(f,3)) )
      allocate( tmp(lbound(f,1):ubound(f,1), lbound(f,2):ubound(f,2), lbound(f,3):ubound(f,3)) )
      allocate( fpom(lbound(f,1):ubound(f,1), lbound(f,2):ubound(f,2), lbound(f,3):ubound(f,3)) )

      fpom(:,:,:) = f(:,:,:)
      cnt = 0

      do while( .true. )

         cnt = cnt + 1 

         ! 1st Runge-Kutta step
         call solve( f, k1, ja, jb, ua, ub, uab, mua, mub, 'GS' )
         ! update vector to the 2nd step
         tmp(:,:,:) = f(:,:,:) + dt/2.0_dp * k1(:,:,:)
         ! 2nd Runge-Kutta step
         call solve( tmp, k2, ja, jb, ua, ub, uab, mua, mub, 'GS' )
         ! update vector for the 3rd step
         tmp(:,:,:) = f(:,:,:) + dt/2.0_dp * k2(:,:,:)
         ! 3rd Runge-Kutta step
         call solve( tmp, k3, ja, jb, ua, ub, uab, mua, mub, 'GS' )
         ! update vector for the 4th step
         tmp(:,:,:) = f(:,:,:) + dt * k3(:,:,:)
         ! 4th Runge-Kutta step
         call solve( tmp, k4, ja, jb, ua, ua, uab, mua, mub, 'GS' )
         ! update vector
         f(:,:,:) = f(:,:,:) + dt/6.0_dp * ( k1(:,:,:) + 2.0_dp *&
                    k2(:,:,:) + 2.0_dp * k3(:,:,:) + k4(:,:,:) )

         ! normalize vector
         call normalize( f )

         ! calculate error
         if( modulo(cnt, stepsForJudge) .eq. 0 ) then

            print *, cnt
            error = 0.0_dp
            do i = 1, ubound(f,1)
               do na = 0, ubound(f,2)
                  do nb = 0, ubound(f,3)
                     ierr = abs( ( f(i,nA,nB) - fpom(i,nA,nB) ) )
                     if (ierr .gt. error) then
                        error = ierr
                     endif
                  enddo
               enddo
            enddo

            if( error .lt. convCriterion ) then
               exit
            endif

            print *, 'error = ', error
            fpom(:,:,:) = f(:,:,:)

         endif

      end do

   end subroutine

!  | ---------------------------------------------- |
!  | Initialize uniformly with unit normalization   | 
!  | ---------------------------------------------- |
   subroutine c1d_InitUniform( f )
 
      implicit none
      complex( kind = dp ), intent(inout) :: f(:,:,:)

      real( kind = dp ) nrm
      integer i, na, nb

      nrm = 1.0_dp/SQRT( 2.0_dp*real( size(f,1),dp ) )

      do i = 1, ubound(f,1)
         do na = 0, ubound(f,2)
            do nb = 0, ubound(f,3)
               f(i,na,nb) = (1.0_dp, 1.0_dp)*nrm
            enddo
         enddo
      enddo

   end subroutine

!  | ---------------------------------------------- |
!  | Initialize uniformly with unit normalization   |
!  | and fixed mean number of particles in each     |
!  | ---------------------------------------------- |
   subroutine c1d_InitUniformNC( f, mna, mnb )

      implicit none
      complex( kind = dp ), intent(inout), allocatable :: f(:,:,:)
      real( kind = dp ), intent(in) :: mna, mnb
      integer info

      f(:,:,:) = (1.0_dp, 1.0_dp)
      call const_mean( f, mna, mnb, info )

   end subroutine

!  | ---------------------------------------------- |
!  | Some private methods needed to maintian        |
!  | constant mean number of particles in each comp |
!  | ---------------------------------------------- |
   subroutine c1d_sigma_norm( f, sig, vout)

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:,:)
      real( kind = dp ), intent(in), dimension(2) :: sig
      real( kind = dp ), intent(inout), allocatable :: vout(:)

      integer i, na, nb
      real( kind = dp ) :: powa, powb, tmp

      do i = 1, ubound(f,1)
         powa = 1.0_dp
         tmp = 0.0_dp
         do na = 0, ubound(f,2)
            powb = 1.0_dp
            do nb = 0, ubound(f,3)
               tmp = tmp + powa * powb * abs(f(i,na,nb))**2
               powb = powb * sig(2)
            enddo
            powa = powa * sig(1)
         enddo
         vout(i) = 1.0_dp/tmp
      enddo

   end subroutine

!  | ---------------------------------------------- |
!  | Project Gutzwiller coefficients onto the       |
!  | submanifold of unit normalization and fixed    |
!  | mean number of particles in each component     |
!  | ---------------------------------------------- |
   subroutine c1d_const_mean( f, mna, mnb, info )

      implicit none
      complex( kind = dp ), intent(inout), allocatable :: f(:,:,:)
      real( kind = dp ) :: mna, mnb
      integer :: info

      real( kind = dp ) :: powa, powb
      real( kind = dp ), dimension(2) :: mu
      integer i, na, nb
      real( kind = dp ), allocatable :: norm(:)

      allocate( norm(lbound(f,1): ubound(f,1)) )
      
      info = 0

      if( mna == 0.0_dp ) then

          mu(2) = 1.0_dp
          call roots( mu(2), f, mnb, c1d_GammaB )
          mu(1) = 0.0_dp
          call sigma_norm( f, mu, norm )

      elseif( mnb == 0.0_dp ) then
          
          mu(1) = 1.0_dp
          call roots( mu(1), f, mna, c1d_GammaA )
          mu(2) = 0.0_dp
          call sigma_norm( f,  mu, norm )

      else

         call roots( mu, f, mna, mnb )
         call sigma_norm( f, mu, norm )

      endif

      if( mu(1) < 0.0_dp  .or. mu(2) < 0.0_dp ) info = 1
      if( mu(1) == 1.0_dp .and. mu(2) == 1.0_dp ) info = 2
      if( ISNAN(mu(1)) .or. ISNAN(mu(2)) ) info = 3
 

      if( info  == 0 ) then
         do i = 1, ubound(f,1)
            powa = 1.0_dp
            do na = 0, ubound(f,2)
               powb = 1.0_dp
               do nb = 0, ubound(f,3)
                  f(i,na,nb) = f(i,na,nb) * sqrt(norm(i)) * powa * powb
                  powb = powb * sqrt(mu(2))
               enddo
               powa = powa * sqrt(mu(1))
            enddo
         enddo
      endif

   end subroutine

end module
