module fixed_mean

   use parameters, only: dp, re, im

   private
   public const_mean

   interface const_mean
       module procedure c1d_const_mean
       module procedure c2d_const_mean
       module procedure c3d_const_mean
   end interface
   
contains

! | ---------------------------------------- |
! | Definition of the Gamma function         |
! | needed for finding proper projectors     |
! | ---------------------------------------- | 

   function c1d_Gamma( x, f, mean )

      implicit none
      real( kind = dp ), intent(in) :: x, mean
      complex( kind = dp ), allocatable, intent(in) :: f(:,:)
      real( kind = dp ), dimension(2) :: c1d_Gamma

      ! local variables
      integer i, na
      real( kind = dp) sum1, sum2, sum3, powa

      c1d_Gamma = 0.0_dp

      if( abs(x) <= 1.0_dp ) then

         do  i = 1, ubound(f,2)
            powa = 1.0_dp
            sum1 = 0.0_dp
            sum2 = 0.0_dp
            sum3 = 0.0_dp
            do na = 0, ubound(f,1)
               sum1 = sum1 + powa * abs(f(na,i))**2
               sum2 = sum2 + powa * real(na,dp) * abs(f(na,i))**2
               sum3 = sum3 + powa * real(na*na,dp) * abs(f(na,i))**2
               powa = powa*x
            enddo

            if( x /= 0.0_dp ) then
               c1d_Gamma(1) = c1d_Gamma(1) + sum2/sum1
               c1d_Gamma(2) = c1d_Gamma(2) + (sum3 - sum2**2/sum1)/(x*sum1)
            else
               c1d_Gamma(2) = c1d_Gamma(2) + abs(f(1,i))**2/abs(f(0,i))**2
            endif

         enddo

      else

         do i = 1, ubound(f,2)
            powa = 1.0_dp
            sum1 = 0.0_dp
            sum2 = 0.0_dp
            sum3 = 0.0_dp
            do na = ubound(f,1), 0, -1
               sum1 = sum1 + powa * abs(f(na,i))**2
               sum2 = sum2 + powa * real(na,dp) * abs(f(na,i))**2
               sum3 = sum3 + powa * real(na*na,dp) * abs(f(na,i))**2
               powa = powa/x
            enddo

            c1d_Gamma(1) = c1d_Gamma(1) + sum2/sum1
            c1d_Gamma(2) = c1d_Gamma(2) + ( sum3 - sum2**2/sum1 )/(x*sum1)

         enddo
 
      endif

      c1d_Gamma(1) = c1d_Gamma(1) - mean

   end function

   function c2d_Gamma( x, f, mean )

      implicit none
      real( kind = dp ), intent(in) :: x, mean
      complex( kind = dp ), allocatable, intent(in) :: f(:,:,:)
      real( kind = dp ), dimension(2) :: c2d_Gamma

      ! local variables
      integer i, j, na
      real( kind = dp) sum1, sum2, sum3, powa

      c2d_Gamma = 0.0_dp

      if( abs(x) <= 1.0_dp ) then

         do j = 1, ubound(f,3)
         do i = 1, ubound(f,2)
            powa = 1.0_dp
            sum1 = 0.0_dp
            sum2 = 0.0_dp
            sum3 = 0.0_dp
            do na = 0, ubound(f,1)
               sum1 = sum1 + powa * abs(f(na,i,j))**2
               sum2 = sum2 + powa * real(na,dp) * abs(f(na,i,j))**2
               sum3 = sum3 + powa * real(na*na,dp) * abs(f(na,i,j))**2
               powa = powa*x
            enddo

            if( x /= 0.0_dp ) then
               c2d_Gamma(1) = c2d_Gamma(1) + sum2/sum1
               c2d_Gamma(2) = c2d_Gamma(2) + (sum3 - sum2**2/sum1)/(x*sum1)
            else
               c2d_Gamma(2) = c2d_Gamma(2) + abs(f(1,i,j))**2/abs(f(0,i,j))**2
            endif

         enddo ! j
         enddo ! i

      else

         do j = 1, ubound(f,3)
         do i = 1, ubound(f,2)
            powa = 1.0_dp
            sum1 = 0.0_dp
            sum2 = 0.0_dp
            sum3 = 0.0_dp
            do na = ubound(f,1), 0, -1
               sum1 = sum1 + powa * abs(f(na,i,j))**2
               sum2 = sum2 + powa * real(na,dp) * abs(f(na,i,j))**2
               sum3 = sum3 + powa * real(na*na,dp) * abs(f(na,i,j))**2
               powa = powa/x
            enddo

            c2d_Gamma(1) = c2d_Gamma(1) + sum2/sum1
            c2d_Gamma(2) = c2d_Gamma(2) + ( sum3 - sum2**2/sum1 )/(x*sum1)

         enddo ! j
         enddo ! i
 
      endif

      c2d_Gamma(1) = c2d_Gamma(1) - mean

   end function

   function c3d_Gamma( x, f, mean )

      implicit none
      real( kind = dp ), intent(in) :: x, mean
      complex( kind = dp ), allocatable, intent(in) :: f(:,:,:,:)
      real( kind = dp ), dimension(2) :: c3d_Gamma

      ! local variables
      integer i, j, k, na
      real( kind = dp) sum1, sum2, sum3, powa

      c3d_Gamma = 0.0_dp

      if( abs(x) <= 1.0_dp ) then

         do k = 1, ubound(f,4)
         do j = 1, ubound(f,3)
         do i = 1, ubound(f,2)
            powa = 1.0_dp
            sum1 = 0.0_dp
            sum2 = 0.0_dp
            sum3 = 0.0_dp
            do na = 0, ubound(f,1)
               sum1 = sum1 + powa * abs(f(na,i,j,k))**2
               sum2 = sum2 + powa * real(na,dp) * abs(f(na,i,j,k))**2
               sum3 = sum3 + powa * real(na*na,dp) * abs(f(na,i,j,k))**2
               powa = powa*x
            enddo

            if( x /= 0.0_dp ) then
               c3d_Gamma(1) = c3d_Gamma(1) + sum2/sum1
               c3d_Gamma(2) = c3d_Gamma(2) + (sum3 - sum2**2/sum1)/(x*sum1)
            else
               c3d_Gamma(2) = c3d_Gamma(2) + abs(f(1,i,j,k))**2/abs(f(0,i,j,k))**2
            endif

         enddo ! i
         enddo ! j
         enddo ! k

      else

         do k = 1, ubound(f,4)
         do j = 1, ubound(f,3)
         do i = 1, ubound(f,2)
            powa = 1.0_dp
            sum1 = 0.0_dp
            sum2 = 0.0_dp
            sum3 = 0.0_dp
            do na = ubound(f,1), 0, -1
               sum1 = sum1 + powa * abs(f(na,i,j,k))**2
               sum2 = sum2 + powa * real(na,dp) * abs(f(na,i,j,k))**2
               sum3 = sum3 + powa * real(na*na,dp) * abs(f(na,i,j,k))**2
               powa = powa/x
            enddo

            c3d_Gamma(1) = c3d_Gamma(1) + sum2/sum1
            c3d_Gamma(2) = c3d_Gamma(2) + ( sum3 - sum2**2/sum1 )/(x*sum1)

         enddo ! i
         enddo ! j
         enddo ! k
 
      endif

      c3d_Gamma(1) = c3d_Gamma(1) - mean

   end function

! | ---------------------------------------- |
! | 1d root finding algortihm                |
! | It combines pre-bisection method and     |
! | then refines with newton algorithm       |
! | ---------------------------------------- | 
   subroutine c1d_root1d( x, f, mn, fun )

      implicit none
      real( kind = dp ), intent(inout) :: x
      complex( kind = dp ), allocatable, intent(in) :: f(:,:)
      real( kind = dp ), intent(in) :: mn

      interface
         function fun( x, f, mn )
            import
            real( kind = dp ), intent(in) :: x
            complex( kind = dp ), allocatable, intent(in) :: f(:,:)
            real( kind = dp ), intent(in) :: mn
            real( kind = dp ), dimension(2) :: fun
         end function
      end interface

      ! local variable
      integer iter
      integer, parameter :: maxiter = 100
      real( kind = dp ), parameter :: conv = 1.0E-15
      real( kind = dp ) a, b, c, dx
      real( kind = dp ), dimension(2) :: fa, fb, fc

      ! Use hybrid method for root finding
      ! Start with bisection and refine with secant

      ! bisection algorithm
      ! we assume that fa * fb <= 0 from the beginning
      a = 0.0_dp
      b = 2.0_dp

      ! check boundaries
      fa = fun( a, f, mn )
      fb = fun( a, f, mn )

      if( fa(1) == 0.0_dp ) then
         x = a
      elseif( fb(1) == 0.0_dp ) then
         x = b
      else

         iter = 0
         do while( iter <= 10 )

            iter = iter + 1
            c  = (a+b)/2.0_dp
            fc = fun( c, f, mn )

            if( fc(1) > 0.0_dp ) b = c
            if( fc(1) < 0.0_dp ) a = c
            if( fc(1) == 0.0_dp ) exit

         enddo

         fc = fun( c, f, mn )
         if( abs(fc(1)) < conv ) then
            x = c
         else
         
            ! refine with Newton algorithm
            do while( iter < maxiter )

               fc = fun( c, f, mn )

               if( fc(2) == 0.0_dp ) then
                  print *, 'c1d_root1d: devided by zero!'
                  exit
               else

                  dx = fc(1)/fc(2)
                  c  = c - dx

               endif

                  iter = iter + 1

               if( abs(dx) < conv ) exit

            enddo

            x = c

          endif

       endif

   end subroutine

   subroutine c2d_root1d( x, f, mn, fun )

      implicit none
      real( kind = dp ), intent(inout) :: x
      complex( kind = dp ), allocatable, intent(in) :: f(:,:,:)
      real( kind = dp ), intent(in) :: mn

      interface
         function fun( x, f, mn )
            import
            real( kind = dp ), intent(in) :: x
            complex( kind = dp ), allocatable, intent(in) :: f(:,:,:)
            real( kind = dp ), intent(in) :: mn
            real( kind = dp ), dimension(2) :: fun
         end function
      end interface

      ! local variable
      integer iter
      integer, parameter :: maxiter = 100
      real( kind = dp ), parameter :: conv = 1.0E-15
      real( kind = dp ) a, b, c, dx
      real( kind = dp ), dimension(2) :: fa, fb, fc

      ! Use hybrid method for root finding
      ! Start with bisection and refine with secant

      ! bisection algorithm
      ! we assume that fa * fb <= 0 from the beginning
      a = 0.0_dp
      b = 2.0_dp

      ! check boundaries
      fa = fun( a, f, mn )
      fb = fun( a, f, mn )

      if( fa(1) == 0.0_dp ) then
         x = a
      elseif( fb(1) == 0.0_dp ) then
         x = b
      else

         iter = 0
         do while( iter <= 10 )

            iter = iter + 1
            c  = (a+b)/2.0_dp
            fc = fun( c, f, mn )

            if( fc(1) > 0.0_dp ) b = c
            if( fc(1) < 0.0_dp ) a = c
            if( fc(1) == 0.0_dp ) exit

         enddo

         fc = fun( c, f, mn )
         if( abs(fc(1)) < conv ) then
            x = c
         else
         
            ! refine with Newton algorithm
            do while( iter < maxiter )

               fc = fun( c, f, mn )

               if( fc(2) == 0.0_dp ) then
                  print *, 'c1d_root1d: devided by zero!'
                  exit
               else

                  dx = fc(1)/fc(2)
                  c  = c - dx

               endif

                  iter = iter + 1

               if( abs(dx) < conv ) exit

            enddo

            x = c

          endif

       endif

   end subroutine

   subroutine c3d_root1d( x, f, mn, fun )

      implicit none
      real( kind = dp ), intent(inout) :: x
      complex( kind = dp ), allocatable, intent(in) :: f(:,:,:,:)
      real( kind = dp ), intent(in) :: mn

      interface
         function fun( x, f, mn )
            import
            real( kind = dp ), intent(in) :: x
            complex( kind = dp ), allocatable, intent(in) :: f(:,:,:,:)
            real( kind = dp ), intent(in) :: mn
            real( kind = dp ), dimension(2) :: fun
         end function
      end interface

      ! local variable
      integer iter
      integer, parameter :: maxiter = 100
      real( kind = dp ), parameter :: conv = 1.0E-15
      real( kind = dp ) a, b, c, dx
      real( kind = dp ), dimension(2) :: fa, fb, fc

      ! Use hybrid method for root finding
      ! Start with bisection and refine with secant

      ! bisection algorithm
      ! we assume that fa * fb <= 0 from the beginning
      a = 0.0_dp
      b = 2.0_dp

      ! check boundaries
      fa = fun( a, f, mn )
      fb = fun( a, f, mn )

      if( fa(1) == 0.0_dp ) then
         x = a
      elseif( fb(1) == 0.0_dp ) then
         x = b
      else

         iter = 0
         do while( iter <= 10 )

            iter = iter + 1
            c  = (a+b)/2.0_dp
            fc = fun( c, f, mn )

            if( fc(1) > 0.0_dp ) b = c
            if( fc(1) < 0.0_dp ) a = c
            if( fc(1) == 0.0_dp ) exit

         enddo

         fc = fun( c, f, mn )
         if( abs(fc(1)) < conv ) then
            x = c
         else
         
            ! refine with Newton algorithm
            do while( iter < maxiter )

               fc = fun( c, f, mn )

               if( fc(2) == 0.0_dp ) then
                  print *, 'c1d_root1d: devided by zero!'
                  exit
               else

                  dx = fc(1)/fc(2)
                  c  = c - dx

               endif

                  iter = iter + 1

               if( abs(dx) < conv ) exit

            enddo

            x = c

          endif

       endif

   end subroutine

!  | ---------------------------------------------- |
!  | Some private methods needed to maintian        |
!  | constant mean number of particles in each comp |
!  | ---------------------------------------------- |
   subroutine c1d_sigma_norm( f, sig, vout)

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:)
      real( kind = dp ), intent(in) :: sig
      real( kind = dp ), intent(inout), allocatable :: vout(:)

      integer i, na
      real( kind = dp ) :: powa, tmp

      do i = 1, ubound(f,2)
         powa = 1.0_dp
         tmp = 0.0_dp
         do na = 0, ubound(f,1)
            tmp = tmp + powa *  abs(f(na,i))**2
            powa = powa * sig
         enddo
         vout(i) = 1.0_dp/tmp
      enddo

   end subroutine

   subroutine c2d_sigma_norm( f, sig, vout)

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:,:)
      real( kind = dp ), intent(in) :: sig
      real( kind = dp ), intent(inout), allocatable :: vout(:,:)

      integer i, j, na
      real( kind = dp ) :: powa, tmp

      do j = 1, ubound(f,3)
         do i = 1, ubound(f,2)
            powa = 1.0_dp
            tmp = 0.0_dp
            do na = 0, ubound(f,1)
               tmp = tmp + powa * abs(f(na,i,j))**2
               powa = powa * sig
            enddo
            vout(i,j) = 1.0_dp/tmp
         enddo
      enddo

   end subroutine

   subroutine c3d_sigma_norm( f, sig, vout)

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:,:,:)
      real( kind = dp ), intent(in) :: sig
      real( kind = dp ), intent(inout), allocatable :: vout(:,:,:)

      integer i, j, k, na
      real( kind = dp ) :: powa, tmp

      do k = 1, ubound(f,4)
         do j = 1, ubound(f,3)
            do i = 1, ubound(f,2)
               powa = 1.0_dp
               tmp = 0.0_dp
               do na = 0, ubound(f,1)
                  tmp = tmp + powa * abs(f(na,i,j,k))**2
                  powa = powa * sig
               enddo
               vout(i,j,k) = 1.0_dp/tmp
            enddo
         enddo
      enddo

   end subroutine

!  | ---------------------------------------------- |
!  | Project Gutzwiller coefficients onto the       |
!  | submanifold of unit normalization and fixed    |
!  | mean number of particles in each component     |
!  | ---------------------------------------------- |
   subroutine c1d_const_mean( f, mn, info )

      implicit none
      complex( kind = dp ), intent(inout), allocatable :: f(:,:)
      real( kind = dp ) :: mn
      integer :: info

      real( kind = dp ) :: powa
      real( kind = dp ) :: mu
      integer i, na
      real( kind = dp ), allocatable :: norm(:)

      allocate( norm(lbound(f,2): ubound(f,2)) )
      
      info = 0


       mu = 1.0_dp
       call c1d_root1d( mu, f, mn, c1d_Gamma )
       call c1d_sigma_norm( f, mu, norm )

      if( mu < 0.0_dp  ) info = 1
      if( mu == 1.0_dp ) info = 2
      if( ISNAN(mu) ) info = 3
 

      if( info  == 0 ) then
         do i = 1, ubound(f,2)
            powa = 1.0_dp
            do na = 0, ubound(f,1)
               f(na,i) = f(na,i) * sqrt(norm(i)) * powa
               powa = powa * sqrt(mu)
            enddo
         enddo
      endif

   end subroutine

   subroutine c2d_const_mean( f, mn, info )

      implicit none
      complex( kind = dp ), intent(inout), allocatable :: f(:,:,:)
      real( kind = dp ) :: mn
      integer :: info

      real( kind = dp ) :: powa
      real( kind = dp ) :: mu
      integer i, j, na
      real( kind = dp ), allocatable :: norm(:,:)

      allocate( norm(lbound(f,2): ubound(f,2), lbound(f,3):ubound(f,3)) )
      
      info = 0

      mu = 1.0_dp
      call c2d_root1d( mu, f, mn, c2d_Gamma )
      call c2d_sigma_norm( f, mu, norm )


      if( mu < 0.0_dp ) info = 1
      if( mu == 1.0_dp ) info = 2
      if( ISNAN(mu) ) info = 3
 

      if( info  == 0 ) then
         do j = 1, ubound(f,3)
         do i = 1, ubound(f,2)
            powa = 1.0_dp
            do na = 0, ubound(f,1)
               f(na,i,j) = f(na,i,j) * sqrt(norm(i,j)) * powa
               powa = powa * sqrt(mu)
            enddo
         enddo ! i
         enddo ! j
      endif

   end subroutine

   subroutine c3d_const_mean( f, mn, info )

      implicit none
      complex( kind = dp ), intent(inout), allocatable :: f(:,:,:,:)
      real( kind = dp ) :: mn
      integer :: info

      real( kind = dp ) :: powa
      real( kind = dp ) :: mu
      integer i, j, k, na
      real( kind = dp ), allocatable :: norm(:,:,:)

      allocate( norm(lbound(f,2): ubound(f,2), lbound(f,3):ubound(f,3), lbound(f,4):ubound(f,4)) )
      
      info = 0

      mu = 1.0_dp
      call c3d_root1d( mu, f, mn, c3d_Gamma )
      call c3d_sigma_norm( f, mu, norm )


      if( mu < 0.0_dp ) info = 1
      if( mu == 1.0_dp ) info = 2
      if( ISNAN(mu) ) info = 3
 

      if( info  == 0 ) then
         do k = 1, ubound(f,4)
         do j = 1, ubound(f,3)
         do i = 1, ubound(f,2)
            powa = 1.0_dp
            do na = 0, ubound(f,1)
               f(na,i,j,k) = f(na,i,j,k) * sqrt(norm(i,j,k)) * powa
               powa = powa * sqrt(mu)
            enddo
         enddo ! i
         enddo ! j
         enddo ! k
      endif

   end subroutine

end module
