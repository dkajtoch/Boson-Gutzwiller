module fixed_mean

   use parameters, only: dp

   interface roots
      module procedure c1d_root1d
      module procedure c1d_root2d
   end interface

   interface minimize
      module procedure c1d_nelmin
   end interface


contains

   function c1d_gammaA( x, f, meana )

      implicit none
      real( kind = dp ), intent(in) :: x, meana
      complex( kind = dp ), allocatable, intent(in) :: f(:,:,:)
      real( kind = dp ), dimension(2) :: c1d_gammaA

      ! local variables
      integer i, na
      real( kind = dp) sum1, sum2, sum3, powa

      c1d_gammaA = 0.0_dp

      if( abs(x) <= 1.0_dp ) then

         do  i = 1, ubound(f,1)
            powa = 1.0_dp
            sum1 = 0.0_dp
            sum2 = 0.0_dp
            sum3 = 0.0_dp
            do na = 0, ubound(f,2)
               sum1 = sum1 + powa * abs(f(i,na,0))**2
               sum2 = sum2 + powa * real(na,dp) * abs(f(i,na,0))**2
               sum3 = sum3 + powa * real(na*na,dp) * abs(f(i,na,0))**2
               powa = powa*x
            enddo

            if( x /= 0.0_dp ) then
               c1d_GammaA(1) = c1d_GammaA(1) + sum2/sum1
               c1d_GammaA(2) = c1d_GammaA(2) + (sum3 - sum2**2/sum1)/(x*sum1)
            else
               c1d_GammaA(2) = c1d_GammaA(2) + abs(f(i,1,0))**2/abs(f(i,0,0))**2
            endif

         enddo

      else

         do i = 1, ubound(f,1)
            powa = 1.0_dp
            sum1 = 0.0_dp
            sum2 = 0.0_dp
            sum3 = 0.0_dp
            do na = ubound(f,2), 0, -1
               sum1 = sum1 + powa * abs(f(i,na,0))**2
               sum2 = sum2 + powa * real(na,dp) * abs(f(i,na,0))**2
               sum3 = sum3 + powa * real(na*na,dp) * abs(f(i,na,0))**2
               powa = powa/x
            enddo

            c1d_GammaA(1) = c1d_GammaA(1) + sum2/sum1
            c1d_GammaA(2) = c1d_GammaA(2) + ( sum3 - sum2**2/sum1 )/(x*sum1)

         enddo
 
      endif

      c1d_GammaA(1) = c1d_GammaA(1) - meana

   end function

   function c1d_gammaB( x, f, meanb )

      implicit none
      real( kind = dp ), intent(in) :: x, meanb
      complex( kind = dp ), allocatable, intent(in) :: f(:,:,:)
      real( kind = dp ), dimension(2) :: c1d_gammaB

      ! local variables
      integer i, nb
      real( kind = dp) sum1, sum2, sum3, powa

      c1d_gammaB = 0.0_dp

      if( abs(x) <= 1.0_dp ) then

         do  i = 1, ubound(f,1)
            powa = 1.0_dp
            sum1 = 0.0_dp
            sum2 = 0.0_dp
            sum3 = 0.0_dp
            do nb = 0, ubound(f,3)
               sum1 = sum1 + powa * abs(f(i,0,nb))**2
               sum2 = sum2 + powa * real(nb,dp) * abs(f(i,0,nb))**2
               sum3 = sum3 + powa * real(nb*nb,dp) * abs(f(i,0,nb))**2
               powa = powa*x
            enddo
 
            if( x == 0 ) then
               c1d_GammaB(2) = c1d_GammaB(2) + abs(f(i,0,1))**2/abs(f(i,0,0))**2
            else
               c1d_GammaB(1) = c1d_GammaB(1) + sum2/sum1
               c1d_GammaB(2) = c1d_GammaB(2) + ( sum3 - sum2**2/sum1 )/(x*sum1)
            endif

         enddo

      else

         do i = 1, ubound(f,1)
            powa = 1.0_dp
            sum1 = 0.0_dp
            sum2 = 0.0_dp
            sum3 = 0.0_dp
            do nb = ubound(f,3), 0, -1
               sum1 = sum1 + powa * abs(f(i,0,nb))**2
               sum2 = sum2 + powa * real(nb,dp) * abs(f(i,0,nb))**2
               sum3 = sum3 + powa * real(nb*nb,dp) * abs(f(i,0,nb))**2
               powa = powa/x
            enddo

            c1d_GammaB(1) = c1d_GammaB(1) + sum2/sum1
            c1d_GammaB(2) = c1d_GammaB(2) + ( sum3 - sum2**2/sum1 )/(x*sum1)

         enddo
 
      endif

      c1d_GammaB(1) = c1d_GammaB(1) - meanb

   end function

   function c1d_GammaAB( x, f, mna, mnb )

      implicit none
      real( kind = dp ), intent(inout) :: x(2)
      complex( kind = dp ), allocatable, intent(in) :: f(:,:,:)
      real( kind = dp ), intent(in) :: mna, mnb
      real( kind = dp), dimension(3,2) :: c1d_GammaAB

      ! local variables
      integer i, na, nb
      real( kind = dp ) sum1, sum2, sum3, sum4, sum5, sum6, powa, powb


      c1d_GammaAB = 0.0_dp

      ! separate definiton of functions into 4 different cases
      ! x(1) == x, x(2) == y
      ! 1) x <= 1.0 and y <= 1.0
      if( abs(x(1)) <= 1.0_dp .and. abs(x(2)) <= 1.0_dp ) then

         ! a) x == 0 and y != 0
         if( x(1) == 0.0_dp .and. x(2) /= 0.0_dp ) then

            do i = 1, ubound(f,1)
               sum1 = 0.0_dp
               sum2 = 0.0_dp
               sum3 = 0.0_dp
               sum4 = 0.0_dp
               sum5 = 0.0_dp
               powb = 1.0_dp
               do nb = 0, ubound(f,3)
                  sum1 = sum1 + powb * abs(f(i,0,nb))**2
                  sum2 = sum2 + real(nb,dp) * powb * abs(f(i,0,nb))**2
                  sum3 = sum3 + powb * abs(f(i,1,nb))**2
                  sum4 = sum4 + real(nb,dp)**2 * powb * abs(f(i,0,nb))**2
                  sum5 = sum5 + real(nb,dp) * powb * abs(f(i,1,nb))**2
                  powb = powb * x(2)
               enddo
               
               c1d_GammaAB(1,2) = c1d_GammaAB(1,2) + sum2/sum1
               c1d_GammaAB(2,1) = c1d_GammaAB(2,1) + sum3/sum1
               c1d_GammaAB(3,1) = c1d_GammaAB(3,1) + sum5/sum1
               c1d_GammaAB(3,2) = c1d_GammaAB(3,2) + (sum4 - sum2**2/sum1)/(sum1 * x(2))

            enddo

         ! b) x != 0 and y == 0        
         elseif( x(1) /= 0.0_dp .and. x(2) == 0.0_dp ) then

            do i = 1, ubound(f,1)
               sum1 = 0.0_dp
               sum2 = 0.0_dp
               sum3 = 0.0_dp
               sum4 = 0.0_dp
               sum5 = 0.0_dp
               powa = 1.0_dp
               do na = 0, ubound(f,2)
                  sum1 = sum1 + powa * abs(f(i,na,0))**2
                  sum2 = sum2 + real(na,dp) * powa * abs(f(i,na,0))**2
                  sum3 = sum3 + powa * abs(f(i,na,1))**2
                  sum4 = sum4 + real(na,dp)**2 * powa * abs(f(i,na,0))**2
                  sum5 = sum5 + real(na,dp) * powa * abs(f(i,na,1))**2
                  powa  = powa * x(1)
               enddo

               c1d_GammaAB(1,1) = c1d_GammaAB(1,1) + sum2/sum1
               c1d_GammaAB(2,1) = c1d_GammaAB(2,1) + (sum4 - sum2**2/sum1)/(sum1 * abs(x(1)))
               c1d_GammaAB(2,2) = c1d_GammaAB(2,2) + sum5/sum1
               c1d_GammaAB(3,2) = c1d_GammaAB(3,2) + sum3/sum1

            enddo

         ! c) x == 0 and y == 0
         elseif( x(1) == 0.0_dp .and. x(2) == 0.0_dp ) then

            do i = 0, ubound(f,1)
               c1d_GammaAB(2,1) = c1d_GammaAB(2,1) + abs(f(i,1,0))**2/abs(f(i,0,0))**2
               c1d_GammaAB(3,2) = c1d_GammaAB(3,2) + abs(f(i,0,1))**2/abs(f(i,0,0))**2
            enddo

         ! d) x != 0 and y != 0
         else
           
            do i = 1, ubound(f,1)
               sum1 = 0.0_dp
               sum2 = 0.0_dp
               sum3 = 0.0_dp
               sum4 = 0.0_dp
               sum5 = 0.0_dp
               sum6 = 0.0_dp
               powa = 1.0_dp
               do na = 0, ubound(f,2)
                  powb = 1.0_dp
                  do nb = 0, ubound(f,3)
                     sum1 = sum1 + powa * powb * abs(f(i,na,nb))**2
                     sum2 = sum2 + real(na,dp) * powa * powb * abs(f(i,na,nb))**2
                     sum3 = sum3 + real(nb,dp) * powa * powb * abs(f(i,na,nb))**2
                     sum4 = sum4 + real(na,dp)**2 * powa * powb * abs(f(i,na,nb))**2
                     sum5 = sum5 + real(nb,dp)**2 * powa * powb * abs(f(i,na,nb))**2
                     sum6 = sum6 + real(na,dp) * real(nb,dp) * powa * powb * abs(f(i,na,nb))**2
                     powb  = powb * x(2)
                  enddo
                  powa = powa * x(1)
               enddo

               c1d_GammaAB(1,1) = c1d_GammaAB(1,1) + sum2/sum1
               c1d_GammaAB(1,2) = c1d_GammaAB(1,2) + sum3/sum1
               c1d_GammaAB(2,1) = c1d_GammaAB(2,1) + (sum4 - sum2**2/sum1)/(sum1 * x(1))
               c1d_GammaAB(2,2) = c1d_GammaAB(2,2) + (sum6 - sum2*sum3/sum1)/(sum1 * x(2))
               c1d_GammaAB(3,1) = c1d_GammaAB(3,1) + (sum6 - sum2*sum3/sum1)/(sum1 * x(1))
               c1d_GammaAB(3,2) = c1d_GammaAB(3,2) + (sum5 - sum3**2/sum1)/(sum1 * x(2))

            enddo
            
         endif 

      ! 2) x <= 1.0 and y > 1.0
      elseif( abs(x(1)) <= 1.0_dp .and. abs(x(2)) > 1.0_dp ) then

         ! a) x = 0
         if( x(1) == 0.0_dp ) then

            do i = 1, ubound(f,1)
               sum1 = 0.0_dp
               sum2 = 0.0_dp
               sum3 = 0.0_dp
               sum4 = 0.0_dp
               sum5 = 0.0_dp
               powb = 1.0_dp
               do nb = ubound(f,3), 0, -1
                  sum1 = sum1 + powb * abs(f(i,0,nb))**2
                  sum2 = sum2 + real(nb,dp) * powb * abs(f(i,0,nb))**2
                  sum3 = sum3 + powb * abs(f(i,1,nb))**2
                  sum4 = sum4 + real(nb,dp)**2 * powb * abs(f(i,0,nb))**2
                  sum5 = sum5 + real(nb,dp) * powb * abs(f(i,1,nb))**2
                  powb  = powb/x(2)
               enddo

               c1d_GammaAB(1,2) = c1d_GammaAB(1,2) + sum2/sum1
               c1d_GammaAB(2,1) = c1d_GammaAB(2,1) + sum3/sum1
               c1d_GammaAB(3,1) = c1d_GammaAB(3,1) + sum5/sum1
               c1d_GammaAB(3,2) = c1d_GammaAB(3,2) + (sum4 - sum2**2/sum1)/(sum1 * x(2))

            enddo

         ! b) x != 0
         else

            do i = 1, ubound(f,1)
               sum1 = 0.0_dp
               sum2 = 0.0_dp
               sum3 = 0.0_dp
               sum4 = 0.0_dp
               sum5 = 0.0_dp
               sum6 = 0.0_dp
               powa = 1.0_dp
               do na = 0, ubound(f,2)
                  powb = 1.0_dp
                  do nb = ubound(f,3), 0, -1
                     sum1 = sum1 + powa * powb * abs(f(i,na,nb))**2
                     sum2 = sum2 + real(na,dp) * powa * powb * abs(f(i,na,nb))**2
                     sum3 = sum3 + real(nb,dp) * powa * powb * abs(f(i,na,nb))**2
                     sum4 = sum4 + real(na,dp)**2 * powa * powb * abs(f(i,na,nb))**2
                     sum5 = sum5 + real(nb,dp)**2 * powa * powb * abs(f(i,na,nb))**2
                     sum6 = sum6 + real(na,dp) * real(nb,dp) * powa * powb * abs(f(i,na,nb))**2
                     powb  = powb/x(2)
                  enddo
                  powa = powa * x(1)
               enddo

               c1d_GammaAB(1,1) = c1d_GammaAB(1,1) + sum2/sum1
               c1d_GammaAB(1,2) = c1d_GammaAB(1,2) + sum3/sum1
               c1d_GammaAB(2,1) = c1d_GammaAB(2,1) + (sum4 - sum2**2/sum1)/(sum1 * x(1))
               c1d_GammaAB(2,2) = c1d_GammaAB(2,2) + (sum6 - sum2*sum3/sum1)/(sum1 * x(2))
               c1d_GammaAB(3,1) = c1d_GammaAB(3,1) + (sum6 - sum2*sum3/sum1)/(sum1 * x(1))
               c1d_GammaAB(3,2) = c1d_GammaAB(3,2) + (sum5 - sum3**2/sum1)/(sum1 * x(2))

            enddo

         endif

      ! 3) x > 1.0 and y <= 1.0
      elseif( abs(x(1)) > 1.0_dp .and. abs(x(2)) <= 1.0_dp ) then

         ! a) y = 0
         if( x(2) == 0.0_dp ) then

            do i = 1, ubound(f,1)
               sum1 = 0.0_dp
               sum2 = 0.0_dp
               sum3 = 0.0_dp
               sum4 = 0.0_dp
               sum5 = 0.0_dp
               powa = 1.0_dp
               do na = ubound(f,2), 0, -1
                  sum1 = sum1 + powa * abs(f(i,na,0))**2
                  sum2 = sum2 + real(na,dp) * powa * abs(f(i,na,0))**2
                  sum3 = sum3 + powa * abs(f(i,na,1))**2
                  sum4 = sum4 + real(na,dp)**2 * powa * abs(f(i,na,0))**2
                  sum5 = sum5 + real(na,dp) * powa * abs(f(i,na,1))**2
                  powa  = powa/x(1)
               enddo

               c1d_GammaAB(1,1) = c1d_GammaAB(1,1) + sum2/sum1
               c1d_GammaAB(2,1) = c1d_GammaAB(2,1) + (sum4 - sum2**2/sum1)/(sum1 * x(1))
               c1d_GammaAB(2,2) = c1d_GammaAB(2,2) + sum5/sum1
               c1d_GammaAB(3,2) = c1d_GammaAB(3,2) + sum3/sum1

            enddo

         ! b) y != 0
         else

            do i = 1, ubound(f,1)
               sum1 = 0.0_dp
               sum2 = 0.0_dp
               sum3 = 0.0_dp
               sum4 = 0.0_dp
               sum5 = 0.0_dp
               sum6 = 0.0_dp
               powa = 1.0_dp
               do na = ubound(f,2), 0,  -1
                  powb = 1.0_dp
                  do nb = 0, ubound(f,3)
                     sum1 = sum1 + powa * powb * abs(f(i,na,nb))**2
                     sum2 = sum2 + real(na,dp) * powa * powb * abs(f(i,na,nb))**2
                     sum3 = sum3 + real(nb,dp) * powa * powb * abs(f(i,na,nb))**2
                     sum4 = sum4 + real(na,dp)**2 * powa * powb * abs(f(i,na,nb))**2
                     sum5 = sum5 + real(nb,dp)**2 * powa * powb * abs(f(i,na,nb))**2
                     sum6 = sum6 + real(na,dp) * real(nb,dp) * powa * powb * abs(f(i,na,nb))**2
                     powb  = powb * x(2)
                  enddo
                  powa = powa/x(1)
               enddo

               c1d_GammaAB(1,1) = c1d_GammaAB(1,1) + sum2/sum1
               c1d_GammaAB(1,2) = c1d_GammaAB(1,2) + sum3/sum1
               c1d_GammaAB(2,1) = c1d_GammaAB(2,1) + (sum4 - sum2**2/sum1)/(sum1 * x(1))
               c1d_GammaAB(2,2) = c1d_GammaAB(2,2) + (sum6 - sum2*sum3/sum1)/(sum1 * x(2))
               c1d_GammaAB(3,1) = c1d_GammaAB(3,1) + (sum6 - sum2*sum3/sum1)/(sum1 * x(1))
               c1d_GammaAB(3,2) = c1d_GammaAB(3,2) + (sum5 - sum3**2/sum1)/(sum1 * x(2))

            enddo

         endif

      ! 4) x > 1.0 and y > 1.0
      else 

         do i = 1, ubound(f,1)
            sum1 = 0.0_dp
            sum2 = 0.0_dp
            sum3 = 0.0_dp 
            sum4 = 0.0_dp
            sum5 = 0.0_dp
            sum6 = 0.0_dp
            powa = 1.0_dp
            do na = ubound(f,2), 0, -1
               powb = 1.0_dp
               do nb = ubound(f,3), 0, -1
                  sum1 = sum1 + powa * powb * abs(f(i,na,nb))**2
                  sum2 = sum2 + real(na,dp) * powa * powb * abs(f(i,na,nb))**2
                  sum3 = sum3 + real(nb,dp) * powa * powb * abs(f(i,na,nb))**2
                  sum4 = sum4 + real(na,dp)**2 * powa * powb * abs(f(i,na,nb))**2
                  sum5 = sum5 + real(nb,dp)**2 * powa * powb * abs(f(i,na,nb))**2
                  sum6 = sum6 + real(na,dp) * real(nb,dp) * powa * powb * abs(f(i,na,nb))**2
                  powb = powb/x(2)
               enddo
               powa = powa/x(1)
            enddo

            c1d_GammaAB(1,1) = c1d_GammaAB(1,1) + sum2/sum1
            c1d_GammaAB(1,2) = c1d_GammaAB(1,2) + sum3/sum1
            c1d_GammaAB(2,1) = c1d_GammaAB(2,1) + (sum4 - sum2**2/sum1)/(sum1 * x(1))
            c1d_GammaAB(2,2) = c1d_GammaAB(2,2) + (sum6 - sum2*sum3/sum1)/(sum1 * x(2))
            c1d_GammaAB(3,1) = c1d_GammaAB(3,1) + (sum6 - sum2*sum3/sum1)/(sum1 * x(1))
            c1d_GammaAB(3,2) = c1d_GammaAB(3,2) + (sum5 - sum3**2/sum1)/(sum1 * x(2))

         enddo

      endif

     c1d_GammaAB(1,1) = c1d_GammaAB(1,1) - mna
     c1d_GammaAB(1,2) = c1d_GammaAB(1,2) - mnb

   end function

   function c1d_nelder_mead( x, f, mna, mnb )

      implicit none
      real( kind = dp ), intent(in), dimension(2) :: x
      complex( kind = dp ), intent(in), allocatable :: f(:,:,:)
      real( kind = dp ), intent(in) :: mna, mnb
      real( kind = dp ) :: c1d_nelder_mead

      real( kind = dp ), parameter :: penalty = 1.0E+14
      real( kind = dp ) :: fx(2)
      real( kind = dp ) sum1, sum2, sum3, powa, powb
      integer i, na, nb

      if( x(1) < 0.0_dp .or. x(2) < 0.0_dp .or. x(1) > 2.0_dp .or. x(2) > 2.0_dp ) then
         c1d_nelder_mead = penalty
      else

         fx = 0.0_dp

         if( abs(x(1)) <= 1.0_dp .and. abs(x(2)) <= 1.0_dp ) then
     
            do i = 1, ubound(f,1)
               sum1 = 0.0_dp
               sum2 = 0.0_dp
               sum3 = 0.0_dp
               powa = 1.0_dp
               do na = 0, ubound(f,2)
                  powb = 1.0_dp
                  do nb = 0, ubound(f,3)
                     sum1 = sum1 + powa * powb * abs(f(i,na,nb))**2
                     sum2 = sum2 + powa * powb * real(na,dp)*abs(f(i,na,nb))**2
                     sum3 = sum3 + powa * powb * real(nb,dp)*abs(f(i,na,nb))**2
                     powb = powb*x(2)
                  enddo
                  powa = powa*x(1)
               enddo

               fx(1) = fx(1) + sum2/sum1
               fx(2) = fx(2) + sum3/sum1

            enddo

         elseif( abs(x(1)) <= 1.0_dp .and. abs(x(2)) > 1.0_dp ) then

            do i = 1, ubound(f,1)
               sum1 = 0.0_dp
               sum2 = 0.0_dp
               sum3 = 0.0_dp
               powa = 1.0_dp
               do na = 0, ubound(f,2)
                  powb = 1.0_dp
                  do nb = ubound(f,3), 0, -1
                     sum1 = sum1 + powa * powb * abs(f(i,na,nb))**2
                     sum2 = sum2 + powa * powb * real(na,dp)*abs(f(i,na,nb))**2
                     sum3 = sum3 + powa * powb * real(nb,dp)*abs(f(i,na,nb))**2
                     powb = powb/x(2)
                  enddo
                  powa = powa*x(1)
               enddo

               fx(1) = fx(1) + sum2/sum1
               fx(2) = fx(2) + sum3/sum1

            enddo

         elseif( abs(x(1)) > 1.0_dp .and. abs(x(2)) <= 1.0_dp ) then

            do i = 1, ubound(f,1)
               sum1 = 0.0_dp
               sum2 = 0.0_dp
               sum3 = 0.0_dp
               powa = 1.0_dp
               do na = ubound(f,2), 0, -1
                  powb = 1.0_dp
                  do nb = 0, ubound(f,3)
                     sum1 = sum1 + powa * powb * abs(f(i,na,nb))**2
                     sum2 = sum2 + powa * powb * real(na,dp)*abs(f(i,na,nb))**2
                     sum3 = sum3 + powa * powb * real(nb,dp)*abs(f(i,na,nb))**2
                     powb = powb*x(2)
                  enddo
                  powa = powa/x(1)
               enddo

               fx(1) = fx(1) + sum2/sum1
               fx(2) = fx(2) + sum3/sum1

            enddo

         else

            do i = 1, ubound(f,1)
               sum1 = 0.0_dp
               sum2 = 0.0_dp
               sum3 = 0.0_dp
               powa = 1.0_dp
               do na = ubound(f,2), 0, -1
                  powb = 1.0_dp
                  do nb = ubound(f,3), 0, -1
                     sum1 = sum1 + powa * powb * abs(f(i,na,nb))**2
                     sum2 = sum2 + powa * powb * real(na,dp)*abs(f(i,na,nb))**2
                     sum3 = sum3 + powa * powb * real(nb,dp)*abs(f(i,na,nb))**2
                     powb = powb/x(2)
                  enddo
                  powa = powa/x(1)
               enddo

               fx(1) = fx(1) + sum2/sum1
               fx(2) = fx(2) + sum3/sum1

            enddo

         endif
 
         fx(1) = fx(1) - mna
         fx(2) = fx(2) - mnb

         c1d_nelder_mead = LOG10( fx(1)**2 + fx(2)**2 + 1.0_dp )

      endif
         
   end function

   subroutine c1d_root1d( x, f, mn, fun )

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

   subroutine c1d_root2d( x, f, mna, mnb )

      implicit none
      real( kind = dp ), intent(inout) :: x(2)
      complex( kind = dp ), allocatable, intent(in) :: f(:,:,:)
      real( kind = dp ), intent(in) :: mna, mnb

      ! local parameters
      real( kind = dp ), dimension(3,2) :: fk
      real( kind = dp ), dimension(2,2) :: invj
      real( kind = dp ) det, dx, dy
      real( kind = dp ), parameter :: conv = 1.0E-15
      integer, parameter :: maxiter = 50
      integer iter

      ! start with nelder-mead
      call minimize( c1d_nelder_mead, x, f, mna, mnb )

      ! refine with Newtons's algorithm
      iter = 0

      do while( iter < maxiter )

         iter = iter + 1

         fk = c1d_GammaAB( x, f, mna, mnb )

         ! calculate inverse jacobian matrix
         det = fk(2,1)*fk(3,2) - fk(3,1)*fk(2,2)

         if( det .eq. 0.0_dp ) then
            print *, 'det equal to 0!'
            exit
         endif

         invj(1,1) =  fk(3,2)/det
         invj(1,2) = -fk(2,2)/det
         invj(2,1) = -fk(3,1)/det
         invj(2,2) =  fk(2,1)/det

         dx = invj(1,1)*fk(1,1) + invj(1,2)*fk(1,2)
         dy = invj(2,1)*fk(1,1) + invj(2,2)*fk(1,2)

         x(1) = x(1) - dx
         x(2) = x(2) - dy

         if( abs(dx) < conv .and. abs(dy) < conv ) exit

      enddo

   end subroutine

   subroutine c1d_nelmin ( fn, xmin, f, mna, mnb )

!*****************************************************************************80
!
!! NELMIN minimizes a function using the Nelder-Mead algorithm.
!
!  Discussion:
!
!    This routine seeks the minimum value of a user-specified function.
!
!    Simplex function minimisation procedure due to Nelder and Mead (1965),
!    as implemented by O'Neill(1971, Appl.Statist. 20, 338-45), with
!    subsequent comments by Chambers+Ertel(1974, 23, 250-1), Benyon(1976,
!    25, 97) and Hill(1978, 27, 380-2)
!
!    The function to be minimized must be defined by a function of
!    the form
!
!      function fn ( x, f )
!      real ( kind = 8 ) fn
!      real ( kind = 8 ) x(*)
!
!    and the name of this subroutine must be declared EXTERNAL in the
!    calling routine and passed as the argument FN.
!
!    This routine does not include a termination test using the
!    fitting of a quadratic surface.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 February 2008
!
!  Author:
!
!    Original FORTRAN77 version by R ONeill.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    John Nelder, Roger Mead,
!    A simplex method for function minimization,
!    Computer Journal,
!    Volume 7, 1965, pages 308-313.
!
!    R ONeill,
!    Algorithm AS 47:
!    Function Minimization Using a Simplex Procedure,
!    Applied Statistics,
!    Volume 20, Number 3, 1971, pages 338-345.
!
!  Parameters:
!
!    Input, external FN, the name of the function which evaluates
!    the function to be minimized.
!
!    Input, integer ( kind = 4 ) N, the number of variables.
!    0 < N is required.
!
!    Input/output, real ( kind = 8 ) START(N).  On input, a starting point
!    for the iteration.  On output, this data may have been overwritten.
!
!    Output, real ( kind = 8 ) XMIN(N), the coordinates of the point which
!    is estimated to minimize the function.
!
!    Output, real ( kind = 8 ) YNEWLO, the minimum value of the function.
!
!    Input, real ( kind = 8 ) REQMIN, the terminating limit for the variance
!    of the function values.  0 < REQMIN is required.
!
!    Input, real ( kind = 8 ) STEP(N), determines the size and shape of the
!    initial simplex.  The relative magnitudes of its elements should reflect
!    the units of the variables.
!
!    Input, integer ( kind = 4 ) KONVGE, the convergence check is carried out
!    every KONVGE iterations. 0 < KONVGE is required.
!
!    Input, integer ( kind = 4 ) KCOUNT, the maximum number of function
!    evaluations.
!
!    Output, integer ( kind = 4 ) ICOUNT, the number of function evaluations
!    used.
!
!    Output, integer ( kind = 4 ) NUMRES, the number of restarts.
!
!    Output, integer ( kind = 4 ) IFAULT, error indicator.
!    0, no errors detected.
!    1, REQMIN, N, or KONVGE has an illegal value.
!    2, iteration terminated because KCOUNT was exceeded without convergence.
!
     implicit none

     complex( kind = dp ), intent(in), allocatable :: f(:,:,:)
     real( kind = dp ), intent(in) :: mna, mnb
     interface 
        function fn( x, f, mna, mnb )
           import
           real( kind = dp ), intent(in) :: x(2)
           complex( kind = dp ), allocatable, intent(in) :: f(:,:,:)
           real( kind = dp ), intent(in) :: mna, mnb
           real( kind = dp ) fn
        end function
     end interface

     ! local parameters
     integer( kind = 4 ), parameter :: n=2 ! 2d problem
     real ( kind = dp ), parameter :: ccoeff = 0.5D+00
     real ( kind = dp ) del
     real ( kind = dp ), parameter :: ecoeff = 2.0D+00
     real ( kind = dp ), parameter :: eps = 0.001D+00
     integer ( kind = 4 ) i
     integer ( kind = 4 ) icount
     integer ( kind = 4 ) ifault
     integer ( kind = 4 ) ihi
     integer ( kind = 4 ) ilo
     integer ( kind = 4 ) j
     integer ( kind = 4 ) jcount
     integer ( kind = 4 ) kcount
     integer ( kind = 4 ) konvge
     integer ( kind = 4 ) l
     integer ( kind = 4 ) numres
     real ( kind = dp ) p(n,n+1)
     real ( kind = dp ) p2star(n)
     real ( kind = dp ) pbar(n)
     real ( kind = dp ) pstar(n)
     real ( kind = dp ), parameter :: rcoeff = 1.0D+00
     real ( kind = dp ) reqmin
     real ( kind = dp ) rq
     real ( kind = dp ) start(n)
     real ( kind = dp ) step(n)
     real ( kind = dp ) x
     real ( kind = dp ) xmin(n)
     real ( kind = dp ) y(n+1)
     real ( kind = dp ) y2star
     real ( kind = dp ) ylo
     real ( kind = dp ) ynewlo
     real ( kind = dp ) ystar
     real ( kind = dp ) z

     ! set intial parameters
     konvge = 10
     kcount = 100
     reqmin = 1.0d-05
     start = (/ 1.0_dp, 1.0_dp /)
     step = (/ -0.2_dp, -0.2_dp /)
     !
     !  Initialization.
     !
     icount = 0
     numres = 0
     jcount = konvge
     del = 1.0D+00
     rq = reqmin * real ( n, kind = dp )
     !
     !  Initial or restarted loop.
     !
     do

        p(1:n,n+1) = start(1:n)
        y(n+1) = fn ( start, f, mna, mnb )
        icount = icount + 1
      !
      !  Define the initial simplex.
      !
      do j = 1, n
         x = start(j)
         start(j) = start(j) + step(j) * del
         p(1:n,j) = start(1:n)
         y(j) = fn ( start, f, mna, mnb )
         icount = icount + 1
         start(j) = x
      enddo
      !
      !  Find highest and lowest Y values.  YNEWLO = Y(IHI) indicates
      !  the vertex of the simplex to be replaced.
      !
      ilo = minloc ( y(1:n+1), 1 )
      ylo = y(ilo)
      !
      !  Inner loop.
      !
      do while ( icount < kcount )
      !
      !  YNEWLO is, of course, the HIGHEST value???
      !
         ihi = maxloc ( y(1:n+1), 1 )
         ynewlo = y(ihi)
!
!  Calculate PBAR, the centroid of the simplex vertices
!  excepting the vertex with Y value YNEWLO.
!
      do i = 1, n
        pbar(i) = ( sum ( p(i,1:n+1) ) - p(i,ihi) ) / real ( n, kind = dp )
      end do
!
!  Reflection through the centroid.
!
      pstar(1:n) = pbar(1:n) + rcoeff * ( pbar(1:n) - p(1:n,ihi) )
      ystar = fn ( pstar, f, mna, mnb )
      icount = icount + 1
!
!  Successful reflection, so extension.
!
      if ( ystar < ylo ) then

        p2star(1:n) = pbar(1:n) + ecoeff * ( pstar(1:n) - pbar(1:n) )
        y2star = fn ( p2star, f, mna, mnb )
        icount = icount + 1
!
!  Retain extension or contraction.
!
        if ( ystar < y2star ) then
          p(1:n,ihi) = pstar(1:n)
          y(ihi) = ystar
        else
          p(1:n,ihi) = p2star(1:n)
          y(ihi) = y2star
        end if
!
!  No extension.
!
      else

        l = 0
        do i = 1, n + 1
          if ( ystar < y(i) ) then
            l = l + 1
          end if
        end do

        if ( 1 < l ) then

          p(1:n,ihi) = pstar(1:n)
          y(ihi) = ystar
!
!  Contraction on the Y(IHI) side of the centroid.
!
        else if ( l == 0 ) then

          p2star(1:n) = pbar(1:n) + ccoeff * ( p(1:n,ihi) - pbar(1:n) )
          y2star = fn ( p2star, f, mna, mnb )
          icount = icount + 1
!
!  Contract the whole simplex.
!
          if ( y(ihi) < y2star ) then

            do j = 1, n + 1
              p(1:n,j) = ( p(1:n,j) + p(1:n,ilo) ) * 0.5D+00
              xmin(1:n) = p(1:n,j)
              y(j) = fn ( xmin, f, mna, mnb )
              icount = icount + 1
            end do

            ilo = minloc ( y(1:n+1), 1 )
            ylo = y(ilo)

            cycle
!
!  Retain contraction.
!
          else
            p(1:n,ihi) = p2star(1:n)
            y(ihi) = y2star
          end if
!
!  Contraction on the reflection side of the centroid.
!
        else if ( l == 1 ) then

          p2star(1:n) = pbar(1:n) + ccoeff * ( pstar(1:n) - pbar(1:n) )
          y2star = fn ( p2star, f, mna, mnb )
          icount = icount + 1
!
!  Retain reflection?
!
          if ( y2star <= ystar ) then
            p(1:n,ihi) = p2star(1:n)
            y(ihi) = y2star
          else
            p(1:n,ihi) = pstar(1:n)
            y(ihi) = ystar
          end if

        end if

      end if
!
!  Check if YLO improved.
!
      if ( y(ihi) < ylo ) then
        ylo = y(ihi)
        ilo = ihi
      end if

      jcount = jcount - 1

      if ( 0 < jcount ) then
        cycle
      end if
!
!  Check to see if minimum reached.
!
      if ( icount <= kcount ) then

        jcount = konvge

        x = sum ( y(1:n+1) ) / real ( n + 1, kind = dp )
        z = sum ( ( y(1:n+1) - x )**2 )

        if ( z <= rq ) then
          exit
        end if

      end if

    end do
!
!  Factorial tests to check that YNEWLO is a local minimum.
!
    xmin(1:n) = p(1:n,ilo)
    ynewlo = y(ilo)

    if ( kcount < icount ) then
      ifault = 2
      exit
    end if

    ifault = 0

    do i = 1, n
      del = step(i) * eps
      xmin(i) = xmin(i) + del
      z = fn ( xmin, f, mna, mnb )
      icount = icount + 1
      if ( z < ynewlo ) then
        ifault = 2
        exit
      end if
      xmin(i) = xmin(i) - del - del
      z = fn ( xmin, f, mna, mnb )
      icount = icount + 1
      if ( z < ynewlo ) then
        ifault = 2
        exit
      end if
      xmin(i) = xmin(i) + del
    end do

    if ( ifault == 0 ) then
      exit
    end if
!
!  Restart the procedure.
!
    start(1:n) = xmin(1:n)
    del = eps
    numres = numres + 1

  end do

  return
end subroutine


end module
