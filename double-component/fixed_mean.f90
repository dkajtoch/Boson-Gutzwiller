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
! | Definition of the GammaA function        |
! | needed for finding proper projectors     |
! | when <nb> = 0                            |
! | ---------------------------------------- | 

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

         do  i = 1, ubound(f,3)
            powa = 1.0_dp
            sum1 = 0.0_dp
            sum2 = 0.0_dp
            sum3 = 0.0_dp
            do na = 0, ubound(f,1)
               sum1 = sum1 + powa * abs(f(na,0,i))**2
               sum2 = sum2 + powa * real(na,dp) * abs(f(na,0,i))**2
               sum3 = sum3 + powa * real(na*na,dp) * abs(f(na,0,i))**2
               powa = powa*x
            enddo

            if( x /= 0.0_dp ) then
               c1d_GammaA(1) = c1d_GammaA(1) + sum2/sum1
               c1d_GammaA(2) = c1d_GammaA(2) + (sum3 - sum2**2/sum1)/(x*sum1)
            else
               c1d_GammaA(2) = c1d_GammaA(2) + abs(f(1,0,i))**2/abs(f(0,0,i))**2
            endif

         enddo

      else

         do i = 1, ubound(f,3)
            powa = 1.0_dp
            sum1 = 0.0_dp
            sum2 = 0.0_dp
            sum3 = 0.0_dp
            do na = ubound(f,1), 0, -1
               sum1 = sum1 + powa * abs(f(na,0,i))**2
               sum2 = sum2 + powa * real(na,dp) * abs(f(na,0,i))**2
               sum3 = sum3 + powa * real(na*na,dp) * abs(f(na,0,i))**2
               powa = powa/x
            enddo

            c1d_GammaA(1) = c1d_GammaA(1) + sum2/sum1
            c1d_GammaA(2) = c1d_GammaA(2) + ( sum3 - sum2**2/sum1 )/(x*sum1)

         enddo
 
      endif

      c1d_GammaA(1) = c1d_GammaA(1) - meana

   end function

   function c2d_gammaA( x, f, meana )

      implicit none
      real( kind = dp ), intent(in) :: x, meana
      complex( kind = dp ), allocatable, intent(in) :: f(:,:,:,:)
      real( kind = dp ), dimension(2) :: c2d_gammaA

      ! local variables
      integer i, j, na
      real( kind = dp) sum1, sum2, sum3, powa

      c2d_gammaA = 0.0_dp

      if( abs(x) <= 1.0_dp ) then

         do j = 1, ubound(f,4)
         do i = 1, ubound(f,3)
            powa = 1.0_dp
            sum1 = 0.0_dp
            sum2 = 0.0_dp
            sum3 = 0.0_dp
            do na = 0, ubound(f,1)
               sum1 = sum1 + powa * abs(f(na,0,i,j))**2
               sum2 = sum2 + powa * real(na,dp) * abs(f(na,0,i,j))**2
               sum3 = sum3 + powa * real(na*na,dp) * abs(f(na,0,i,j))**2
               powa = powa*x
            enddo

            if( x /= 0.0_dp ) then
               c2d_GammaA(1) = c2d_GammaA(1) + sum2/sum1
               c2d_GammaA(2) = c2d_GammaA(2) + (sum3 - sum2**2/sum1)/(x*sum1)
            else
               c2d_GammaA(2) = c2d_GammaA(2) + abs(f(1,0,i,j))**2/abs(f(0,0,i,j))**2
            endif

         enddo ! i
         enddo ! j

      else

         do j = 1, ubound(f,4)
         do i = 1, ubound(f,3)
            powa = 1.0_dp
            sum1 = 0.0_dp
            sum2 = 0.0_dp
            sum3 = 0.0_dp
            do na = ubound(f,1), 0, -1
               sum1 = sum1 + powa * abs(f(na,0,i,j))**2
               sum2 = sum2 + powa * real(na,dp) * abs(f(na,0,i,j))**2
               sum3 = sum3 + powa * real(na*na,dp) * abs(f(na,0,i,j))**2
               powa = powa/x
            enddo

            c2d_GammaA(1) = c2d_GammaA(1) + sum2/sum1
            c2d_GammaA(2) = c2d_GammaA(2) + ( sum3 - sum2**2/sum1 )/(x*sum1)

         enddo ! i
         enddo ! j
 
      endif

      c2d_GammaA(1) = c2d_GammaA(1) - meana

   end function

   function c3d_gammaA( x, f, meana )

      implicit none
      real( kind = dp ), intent(in) :: x, meana
      complex( kind = dp ), allocatable, intent(in) :: f(:,:,:,:,:)
      real( kind = dp ), dimension(2) :: c3d_gammaA

      ! local variables
      integer i, j, k, na
      real( kind = dp) sum1, sum2, sum3, powa

      c3d_gammaA = 0.0_dp

      if( abs(x) <= 1.0_dp ) then

         do k = 1, ubound(f,5)
         do j = 1, ubound(f,4)
         do i = 1, ubound(f,3)
            powa = 1.0_dp
            sum1 = 0.0_dp
            sum2 = 0.0_dp
            sum3 = 0.0_dp
            do na = 0, ubound(f,1)
               sum1 = sum1 + powa * abs(f(na,0,i,j,k))**2
               sum2 = sum2 + powa * real(na,dp) * abs(f(na,0,i,j,k))**2
               sum3 = sum3 + powa * real(na*na,dp) * abs(f(na,0,i,j,k))**2
               powa = powa*x
            enddo

            if( x /= 0.0_dp ) then
               c3d_GammaA(1) = c3d_GammaA(1) + sum2/sum1
               c3d_GammaA(2) = c3d_GammaA(2) + (sum3 - sum2**2/sum1)/(x*sum1)
            else
               c3d_GammaA(2) = c3d_GammaA(2) + abs(f(1,0,i,j,k))**2/abs(f(0,0,i,j,k))**2
            endif

         enddo ! i
         enddo ! j
         enddo ! k

      else

         do k = 1, ubound(f,5)
         do j = 1, ubound(f,4)
         do i = 1, ubound(f,3)
            powa = 1.0_dp
            sum1 = 0.0_dp
            sum2 = 0.0_dp
            sum3 = 0.0_dp
            do na = ubound(f,1), 0, -1
               sum1 = sum1 + powa * abs(f(na,0,i,j,k))**2
               sum2 = sum2 + powa * real(na,dp) * abs(f(na,0,i,j,k))**2
               sum3 = sum3 + powa * real(na*na,dp) * abs(f(na,0,i,j,k))**2
               powa = powa/x
            enddo

            c3d_GammaA(1) = c3d_GammaA(1) + sum2/sum1
            c3d_GammaA(2) = c3d_GammaA(2) + ( sum3 - sum2**2/sum1 )/(x*sum1)

         enddo ! i
         enddo ! j
         enddo ! k
 
      endif

      c3d_GammaA(1) = c3d_GammaA(1) - meana

   end function

! | ---------------------------------------- |
! | Definition of the GammaA function        |
! | needed for finding proper projectors     |
! | when <na> = 0                            |
! | ---------------------------------------- | 

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

         do i = 1, ubound(f,3)
            powa = 1.0_dp
            sum1 = 0.0_dp
            sum2 = 0.0_dp
            sum3 = 0.0_dp
            do nb = 0, ubound(f,2)
               sum1 = sum1 + powa * abs(f(0,nb,i))**2
               sum2 = sum2 + powa * real(nb,dp) * abs(f(0,nb,i))**2
               sum3 = sum3 + powa * real(nb*nb,dp) * abs(f(0,nb,i))**2
               powa = powa*x
            enddo
 
            if( x == 0 ) then
               c1d_GammaB(2) = c1d_GammaB(2) + abs(f(0,1,i))**2/abs(f(0,0,i))**2
            else
               c1d_GammaB(1) = c1d_GammaB(1) + sum2/sum1
               c1d_GammaB(2) = c1d_GammaB(2) + ( sum3 - sum2**2/sum1 )/(x*sum1)
            endif

         enddo

      else

         do i = 1, ubound(f,3)
            powa = 1.0_dp
            sum1 = 0.0_dp
            sum2 = 0.0_dp
            sum3 = 0.0_dp
            do nb = ubound(f,2), 0, -1
               sum1 = sum1 + powa * abs(f(0,nb,i))**2
               sum2 = sum2 + powa * real(nb,dp) * abs(f(0,nb,i))**2
               sum3 = sum3 + powa * real(nb*nb,dp) * abs(f(0,nb,i))**2
               powa = powa/x
            enddo

            c1d_GammaB(1) = c1d_GammaB(1) + sum2/sum1
            c1d_GammaB(2) = c1d_GammaB(2) + ( sum3 - sum2**2/sum1 )/(x*sum1)

         enddo
 
      endif

      c1d_GammaB(1) = c1d_GammaB(1) - meanb

   end function

   function c2d_gammaB( x, f, meanb )

      implicit none
      real( kind = dp ), intent(in) :: x, meanb
      complex( kind = dp ), allocatable, intent(in) :: f(:,:,:,:)
      real( kind = dp ), dimension(2) :: c2d_gammaB

      ! local variables
      integer i, j, nb
      real( kind = dp) sum1, sum2, sum3, powa

      c2d_gammaB = 0.0_dp

      if( abs(x) <= 1.0_dp ) then

         do j = 1, ubound(f,4)
         do i = 1, ubound(f,3)
            powa = 1.0_dp
            sum1 = 0.0_dp
            sum2 = 0.0_dp
            sum3 = 0.0_dp
            do nb = 0, ubound(f,2)
               sum1 = sum1 + powa * abs(f(0,nb,i,j))**2
               sum2 = sum2 + powa * real(nb,dp) * abs(f(0,nb,i,j))**2
               sum3 = sum3 + powa * real(nb*nb,dp) * abs(f(0,nb,i,j))**2
               powa = powa*x
            enddo
 
            if( x == 0 ) then
               c2d_GammaB(2) = c2d_GammaB(2) + abs(f(0,1,i,j))**2/abs(f(0,0,i,j))**2
            else
               c2d_GammaB(1) = c2d_GammaB(1) + sum2/sum1
               c2d_GammaB(2) = c2d_GammaB(2) + ( sum3 - sum2**2/sum1 )/(x*sum1)
            endif

         enddo ! i
         enddo ! j

      else

         do j = 1, ubound(f,4)
         do i = 1, ubound(f,3)
            powa = 1.0_dp
            sum1 = 0.0_dp
            sum2 = 0.0_dp
            sum3 = 0.0_dp
            do nb = ubound(f,2), 0, -1
               sum1 = sum1 + powa * abs(f(0,nb,i,j))**2
               sum2 = sum2 + powa * real(nb,dp) * abs(f(0,nb,i,j))**2
               sum3 = sum3 + powa * real(nb*nb,dp) * abs(f(0,nb,i,j))**2
               powa = powa/x
            enddo

            c2d_GammaB(1) = c2d_GammaB(1) + sum2/sum1
            c2d_GammaB(2) = c2d_GammaB(2) + ( sum3 - sum2**2/sum1 )/(x*sum1)

         enddo ! i
         enddo ! j
 
      endif

      c2d_GammaB(1) = c2d_GammaB(1) - meanb

   end function

   function c3d_gammaB( x, f, meanb )

      implicit none
      real( kind = dp ), intent(in) :: x, meanb
      complex( kind = dp ), allocatable, intent(in) :: f(:,:,:,:,:)
      real( kind = dp ), dimension(2) :: c3d_gammaB

      ! local variables
      integer i, j, k, nb
      real( kind = dp) sum1, sum2, sum3, powa

      c3d_gammaB = 0.0_dp

      if( abs(x) <= 1.0_dp ) then

         do k = 1, ubound(f,5)
         do j = 1, ubound(f,4)
         do i = 1, ubound(f,3)
            powa = 1.0_dp
            sum1 = 0.0_dp
            sum2 = 0.0_dp
            sum3 = 0.0_dp
            do nb = 0, ubound(f,2)
               sum1 = sum1 + powa * abs(f(0,nb,i,j,k))**2
               sum2 = sum2 + powa * real(nb,dp) * abs(f(0,nb,i,j,k))**2
               sum3 = sum3 + powa * real(nb*nb,dp) * abs(f(0,nb,i,j,k))**2
               powa = powa*x
            enddo
 
            if( x == 0 ) then
               c3d_GammaB(2) = c3d_GammaB(2) + abs(f(0,1,i,j,k))**2/abs(f(0,0,i,j,k))**2
            else
               c3d_GammaB(1) = c3d_GammaB(1) + sum2/sum1
               c3d_GammaB(2) = c3d_GammaB(2) + ( sum3 - sum2**2/sum1 )/(x*sum1)
            endif

         enddo ! i
         enddo ! j
         enddo ! k

      else

         do k = 1, ubound(f,5)
         do j = 1, ubound(f,4)
         do i = 1, ubound(f,3)
            powa = 1.0_dp
            sum1 = 0.0_dp
            sum2 = 0.0_dp
            sum3 = 0.0_dp
            do nb = ubound(f,2), 0, -1
               sum1 = sum1 + powa * abs(f(0,nb,i,j,k))**2
               sum2 = sum2 + powa * real(nb,dp) * abs(f(0,nb,i,j,k))**2
               sum3 = sum3 + powa * real(nb*nb,dp) * abs(f(0,nb,i,j,k))**2
               powa = powa/x
            enddo

            c3d_GammaB(1) = c3d_GammaB(1) + sum2/sum1
            c3d_GammaB(2) = c3d_GammaB(2) + ( sum3 - sum2**2/sum1 )/(x*sum1)

         enddo ! i
         enddo ! j
         enddo ! k
 
      endif

      c3d_GammaB(1) = c3d_GammaB(1) - meanb

   end function

! | ---------------------------------------- |
! | Definition of the GammaAB function       |
! | needed for finding proper projectors     |
! | when both  <na>, <nb> are non-zero       |
! | ---------------------------------------- | 

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

            do i = 1, ubound(f,3)
               sum1 = 0.0_dp
               sum2 = 0.0_dp
               sum3 = 0.0_dp
               sum4 = 0.0_dp
               sum5 = 0.0_dp
               powb = 1.0_dp
               do nb = 0, ubound(f,2)
                  sum1 = sum1 + powb * abs(f(0,nb,i))**2
                  sum2 = sum2 + real(nb,dp) * powb * abs(f(0,nb,i))**2
                  sum3 = sum3 + powb * abs(f(1,nb,i))**2
                  sum4 = sum4 + real(nb,dp)**2 * powb * abs(f(0,nb,i))**2
                  sum5 = sum5 + real(nb,dp) * powb * abs(f(1,nb,i))**2
                  powb = powb * x(2)
               enddo
               
               c1d_GammaAB(1,2) = c1d_GammaAB(1,2) + sum2/sum1
               c1d_GammaAB(2,1) = c1d_GammaAB(2,1) + sum3/sum1
               c1d_GammaAB(3,1) = c1d_GammaAB(3,1) + sum5/sum1
               c1d_GammaAB(3,2) = c1d_GammaAB(3,2) + (sum4 - sum2**2/sum1)/(sum1 * x(2))

            enddo

         ! b) x != 0 and y == 0        
         elseif( x(1) /= 0.0_dp .and. x(2) == 0.0_dp ) then

            do i = 1, ubound(f,3)
               sum1 = 0.0_dp
               sum2 = 0.0_dp
               sum3 = 0.0_dp
               sum4 = 0.0_dp
               sum5 = 0.0_dp
               powa = 1.0_dp
               do na = 0, ubound(f,1)
                  sum1 = sum1 + powa * abs(f(na,0,i))**2
                  sum2 = sum2 + real(na,dp) * powa * abs(f(na,0,i))**2
                  sum3 = sum3 + powa * abs(f(na,1,i))**2
                  sum4 = sum4 + real(na,dp)**2 * powa * abs(f(na,0,i))**2
                  sum5 = sum5 + real(na,dp) * powa * abs(f(na,1,i))**2
                  powa  = powa * x(1)
               enddo

               c1d_GammaAB(1,1) = c1d_GammaAB(1,1) + sum2/sum1
               c1d_GammaAB(2,1) = c1d_GammaAB(2,1) + (sum4 - sum2**2/sum1)/(sum1 * abs(x(1)))
               c1d_GammaAB(2,2) = c1d_GammaAB(2,2) + sum5/sum1
               c1d_GammaAB(3,2) = c1d_GammaAB(3,2) + sum3/sum1

            enddo

         ! c) x == 0 and y == 0
         elseif( x(1) == 0.0_dp .and. x(2) == 0.0_dp ) then

            do i = 1, ubound(f,3)
               c1d_GammaAB(2,1) = c1d_GammaAB(2,1) + abs(f(1,0,i))**2/abs(f(0,0,i))**2
               c1d_GammaAB(3,2) = c1d_GammaAB(3,2) + abs(f(0,1,i))**2/abs(f(0,0,i))**2
            enddo

         ! d) x != 0 and y != 0
         else
           
            do i = 1, ubound(f,3)
               sum1 = 0.0_dp
               sum2 = 0.0_dp
               sum3 = 0.0_dp
               sum4 = 0.0_dp
               sum5 = 0.0_dp
               sum6 = 0.0_dp
               powb = 1.0_dp
               do nb = 0, ubound(f,2)
                  powa = 1.0_dp
                  do na = 0, ubound(f,1)
                     sum1 = sum1 + powa * powb * abs(f(na,nb,i))**2
                     sum2 = sum2 + real(na,dp) * powa * powb * abs(f(na,nb,i))**2
                     sum3 = sum3 + real(nb,dp) * powa * powb * abs(f(na,nb,i))**2
                     sum4 = sum4 + real(na,dp)**2 * powa * powb * abs(f(na,nb,i))**2
                     sum5 = sum5 + real(nb,dp)**2 * powa * powb * abs(f(na,nb,i))**2
                     sum6 = sum6 + real(na,dp) * real(nb,dp) * powa * powb * abs(f(na,nb,i))**2
                     powa  = powa * x(1)
                  enddo
                  powb = powb * x(2)
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

            do i = 1, ubound(f,3)
               sum1 = 0.0_dp
               sum2 = 0.0_dp
               sum3 = 0.0_dp
               sum4 = 0.0_dp
               sum5 = 0.0_dp
               powb = 1.0_dp
               do nb = ubound(f,2), 0, -1
                  sum1 = sum1 + powb * abs(f(0,nb,i))**2
                  sum2 = sum2 + real(nb,dp) * powb * abs(f(0,nb,i))**2
                  sum3 = sum3 + powb * abs(f(1,nb,i))**2
                  sum4 = sum4 + real(nb,dp)**2 * powb * abs(f(0,nb,i))**2
                  sum5 = sum5 + real(nb,dp) * powb * abs(f(1,nb,i))**2
                  powb  = powb/x(2)
               enddo

               c1d_GammaAB(1,2) = c1d_GammaAB(1,2) + sum2/sum1
               c1d_GammaAB(2,1) = c1d_GammaAB(2,1) + sum3/sum1
               c1d_GammaAB(3,1) = c1d_GammaAB(3,1) + sum5/sum1
               c1d_GammaAB(3,2) = c1d_GammaAB(3,2) + (sum4 - sum2**2/sum1)/(sum1 * x(2))

            enddo

         ! b) x != 0
         else

            do i = 1, ubound(f,3)
               sum1 = 0.0_dp
               sum2 = 0.0_dp
               sum3 = 0.0_dp
               sum4 = 0.0_dp
               sum5 = 0.0_dp
               sum6 = 0.0_dp
               powb = 1.0_dp
               do nb = ubound(f,2), 0, -1
                  powa = 1.0_dp
                  do na = 0, ubound(f,1)
                     sum1 = sum1 + powa * powb * abs(f(na,nb,i))**2
                     sum2 = sum2 + real(na,dp) * powa * powb * abs(f(na,nb,i))**2
                     sum3 = sum3 + real(nb,dp) * powa * powb * abs(f(na,nb,i))**2
                     sum4 = sum4 + real(na,dp)**2 * powa * powb * abs(f(na,nb,i))**2
                     sum5 = sum5 + real(nb,dp)**2 * powa * powb * abs(f(na,nb,i))**2
                     sum6 = sum6 + real(na,dp) * real(nb,dp) * powa * powb * abs(f(na,nb,i))**2
                     powa  = powa * x(1)
                  enddo
                  powb = powb / x(2)
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

            do i = 1, ubound(f,3)
               sum1 = 0.0_dp
               sum2 = 0.0_dp
               sum3 = 0.0_dp
               sum4 = 0.0_dp
               sum5 = 0.0_dp
               powa = 1.0_dp
               do na = ubound(f,1), 0, -1
                  sum1 = sum1 + powa * abs(f(na,0,i))**2
                  sum2 = sum2 + real(na,dp) * powa * abs(f(na,0,i))**2
                  sum3 = sum3 + powa * abs(f(na,1,i))**2
                  sum4 = sum4 + real(na,dp)**2 * powa * abs(f(na,0,i))**2
                  sum5 = sum5 + real(na,dp) * powa * abs(f(na,1,i))**2
                  powa  = powa/x(1)
               enddo

               c1d_GammaAB(1,1) = c1d_GammaAB(1,1) + sum2/sum1
               c1d_GammaAB(2,1) = c1d_GammaAB(2,1) + (sum4 - sum2**2/sum1)/(sum1 * x(1))
               c1d_GammaAB(2,2) = c1d_GammaAB(2,2) + sum5/sum1
               c1d_GammaAB(3,2) = c1d_GammaAB(3,2) + sum3/sum1

            enddo

         ! b) y != 0
         else

            do i = 1, ubound(f,3)
               sum1 = 0.0_dp
               sum2 = 0.0_dp
               sum3 = 0.0_dp
               sum4 = 0.0_dp
               sum5 = 0.0_dp
               sum6 = 0.0_dp
               powb = 1.0_dp
               do nb = 0, ubound(f,2)
                  powa = 1.0_dp
                  do na = ubound(f,1), 0, -1
                     sum1 = sum1 + powa * powb * abs(f(na,nb,i))**2
                     sum2 = sum2 + real(na,dp) * powa * powb * abs(f(na,nb,i))**2
                     sum3 = sum3 + real(nb,dp) * powa * powb * abs(f(na,nb,i))**2
                     sum4 = sum4 + real(na,dp)**2 * powa * powb * abs(f(na,nb,i))**2
                     sum5 = sum5 + real(nb,dp)**2 * powa * powb * abs(f(na,nb,i))**2
                     sum6 = sum6 + real(na,dp) * real(nb,dp) * powa * powb * abs(f(na,nb,i))**2
                     powa  = powa / x(1)
                  enddo
                  powb = powb * x(2)
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

         do i = 1, ubound(f,3)
            sum1 = 0.0_dp
            sum2 = 0.0_dp
            sum3 = 0.0_dp 
            sum4 = 0.0_dp
            sum5 = 0.0_dp
            sum6 = 0.0_dp
            powb = 1.0_dp
            do nb = ubound(f,2), 0, -1
               powa = 1.0_dp
               do na = ubound(f,1), 0, -1
                  sum1 = sum1 + powa * powb * abs(f(na,nb,i))**2
                  sum2 = sum2 + real(na,dp) * powa * powb * abs(f(na,nb,i))**2
                  sum3 = sum3 + real(nb,dp) * powa * powb * abs(f(na,nb,i))**2
                  sum4 = sum4 + real(na,dp)**2 * powa * powb * abs(f(na,nb,i))**2
                  sum5 = sum5 + real(nb,dp)**2 * powa * powb * abs(f(na,nb,i))**2
                  sum6 = sum6 + real(na,dp) * real(nb,dp) * powa * powb * abs(f(na,nb,i))**2
                  powa = powa / x(2)
               enddo
               powb = powb / x(1)
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

   function c2d_GammaAB( x, f, mna, mnb )

      implicit none
      real( kind = dp ), intent(inout) :: x(2)
      complex( kind = dp ), allocatable, intent(in) :: f(:,:,:,:)
      real( kind = dp ), intent(in) :: mna, mnb
      real( kind = dp), dimension(3,2) :: c2d_GammaAB

      ! local variables
      integer i, j, na, nb
      real( kind = dp ) sum1, sum2, sum3, sum4, sum5, sum6, powa, powb


      c2d_GammaAB = 0.0_dp

      ! separate definiton of functions into 4 different cases
      ! x(1) == x, x(2) == y
      ! 1) x <= 1.0 and y <= 1.0
      if( abs(x(1)) <= 1.0_dp .and. abs(x(2)) <= 1.0_dp ) then

         ! a) x == 0 and y != 0
         if( x(1) == 0.0_dp .and. x(2) /= 0.0_dp ) then

            do j = 1, ubound(f,4)
            do i = 1, ubound(f,3)
               sum1 = 0.0_dp
               sum2 = 0.0_dp
               sum3 = 0.0_dp
               sum4 = 0.0_dp
               sum5 = 0.0_dp
               powb = 1.0_dp
               do nb = 0, ubound(f,2)
                  sum1 = sum1 + powb * abs(f(0,nb,i,j))**2
                  sum2 = sum2 + real(nb,dp) * powb * abs(f(0,nb,i,j))**2
                  sum3 = sum3 + powb * abs(f(1,nb,i,j))**2
                  sum4 = sum4 + real(nb,dp)**2 * powb * abs(f(0,nb,i,j))**2
                  sum5 = sum5 + real(nb,dp) * powb * abs(f(1,nb,i,j))**2
                  powb = powb * x(2)
               enddo
               
               c2d_GammaAB(1,2) = c2d_GammaAB(1,2) + sum2/sum1
               c2d_GammaAB(2,1) = c2d_GammaAB(2,1) + sum3/sum1
               c2d_GammaAB(3,1) = c2d_GammaAB(3,1) + sum5/sum1
               c2d_GammaAB(3,2) = c2d_GammaAB(3,2) + (sum4 - sum2**2/sum1)/(sum1 * x(2))

            enddo ! i
            enddo ! j

         ! b) x != 0 and y == 0        
         elseif( x(1) /= 0.0_dp .and. x(2) == 0.0_dp ) then

            do j = 1, ubound(f,4)
            do i = 1, ubound(f,3)
               sum1 = 0.0_dp
               sum2 = 0.0_dp
               sum3 = 0.0_dp
               sum4 = 0.0_dp
               sum5 = 0.0_dp
               powa = 1.0_dp
               do na = 0, ubound(f,1)
                  sum1 = sum1 + powa * abs(f(na,0,i,j))**2
                  sum2 = sum2 + real(na,dp) * powa * abs(f(na,0,i,j))**2
                  sum3 = sum3 + powa * abs(f(na,1,i,j))**2
                  sum4 = sum4 + real(na,dp)**2 * powa * abs(f(na,0,i,j))**2
                  sum5 = sum5 + real(na,dp) * powa * abs(f(na,1,i,j))**2
                  powa  = powa * x(1)
               enddo

               c2d_GammaAB(1,1) = c2d_GammaAB(1,1) + sum2/sum1
               c2d_GammaAB(2,1) = c2d_GammaAB(2,1) + (sum4 - sum2**2/sum1)/(sum1 * abs(x(1)))
               c2d_GammaAB(2,2) = c2d_GammaAB(2,2) + sum5/sum1
               c2d_GammaAB(3,2) = c2d_GammaAB(3,2) + sum3/sum1

            enddo ! i
            enddo ! j

         ! c) x == 0 and y == 0
         elseif( x(1) == 0.0_dp .and. x(2) == 0.0_dp ) then

            do j = 1, ubound(f,4)
            do i = 1, ubound(f,3)
               c2d_GammaAB(2,1) = c2d_GammaAB(2,1) + abs(f(1,0,i,j))**2/abs(f(0,0,i,j))**2
               c2d_GammaAB(3,2) = c2d_GammaAB(3,2) + abs(f(0,1,i,j))**2/abs(f(0,0,i,j))**2
            enddo ! i
            enddo ! j

         ! d) x != 0 and y != 0
         else
           
            do j = 1, ubound(f,4)
            do i = 1, ubound(f,3)
               sum1 = 0.0_dp
               sum2 = 0.0_dp
               sum3 = 0.0_dp
               sum4 = 0.0_dp
               sum5 = 0.0_dp
               sum6 = 0.0_dp
               powb = 1.0_dp
               do nb = 0, ubound(f,2)
                  powa = 1.0_dp
                  do na = 0, ubound(f,1)
                     sum1 = sum1 + powa * powb * abs(f(na,nb,i,j))**2
                     sum2 = sum2 + real(na,dp) * powa * powb * abs(f(na,nb,i,j))**2
                     sum3 = sum3 + real(nb,dp) * powa * powb * abs(f(na,nb,i,j))**2
                     sum4 = sum4 + real(na,dp)**2 * powa * powb * abs(f(na,nb,i,j))**2
                     sum5 = sum5 + real(nb,dp)**2 * powa * powb * abs(f(na,nb,i,j))**2
                     sum6 = sum6 + real(na,dp) * real(nb,dp) * powa * powb * abs(f(na,nb,i,j))**2
                     powa  = powa * x(1)
                  enddo
                  powb = powb * x(2)
               enddo

               c2d_GammaAB(1,1) = c2d_GammaAB(1,1) + sum2/sum1
               c2d_GammaAB(1,2) = c2d_GammaAB(1,2) + sum3/sum1
               c2d_GammaAB(2,1) = c2d_GammaAB(2,1) + (sum4 - sum2**2/sum1)/(sum1 * x(1))
               c2d_GammaAB(2,2) = c2d_GammaAB(2,2) + (sum6 - sum2*sum3/sum1)/(sum1 * x(2))
               c2d_GammaAB(3,1) = c2d_GammaAB(3,1) + (sum6 - sum2*sum3/sum1)/(sum1 * x(1))
               c2d_GammaAB(3,2) = c2d_GammaAB(3,2) + (sum5 - sum3**2/sum1)/(sum1 * x(2))

            enddo ! i
            enddo ! j
            
         endif 

      ! 2) x <= 1.0 and y > 1.0
      elseif( abs(x(1)) <= 1.0_dp .and. abs(x(2)) > 1.0_dp ) then

         ! a) x = 0
         if( x(1) == 0.0_dp ) then

            do j = 1, ubound(f,4)
            do i = 1, ubound(f,3)
               sum1 = 0.0_dp
               sum2 = 0.0_dp
               sum3 = 0.0_dp
               sum4 = 0.0_dp
               sum5 = 0.0_dp
               powb = 1.0_dp
               do nb = ubound(f,2), 0, -1
                  sum1 = sum1 + powb * abs(f(0,nb,i,j))**2
                  sum2 = sum2 + real(nb,dp) * powb * abs(f(0,nb,i,j))**2
                  sum3 = sum3 + powb * abs(f(1,nb,i,j))**2
                  sum4 = sum4 + real(nb,dp)**2 * powb * abs(f(0,nb,i,j))**2
                  sum5 = sum5 + real(nb,dp) * powb * abs(f(1,nb,i,j))**2
                  powb  = powb/x(2)
               enddo

               c2d_GammaAB(1,2) = c2d_GammaAB(1,2) + sum2/sum1
               c2d_GammaAB(2,1) = c2d_GammaAB(2,1) + sum3/sum1
               c2d_GammaAB(3,1) = c2d_GammaAB(3,1) + sum5/sum1
               c2d_GammaAB(3,2) = c2d_GammaAB(3,2) + (sum4 - sum2**2/sum1)/(sum1 * x(2))

            enddo ! i
            enddo ! j

         ! b) x != 0
         else

            do j = 1, ubound(f,4)
            do i = 1, ubound(f,3)
               sum1 = 0.0_dp
               sum2 = 0.0_dp
               sum3 = 0.0_dp
               sum4 = 0.0_dp
               sum5 = 0.0_dp
               sum6 = 0.0_dp
               powb = 1.0_dp
               do nb = ubound(f,2), 0, -1
                  powa = 1.0_dp
                  do na = 0, ubound(f,1)
                     sum1 = sum1 + powa * powb * abs(f(na,nb,i,j))**2
                     sum2 = sum2 + real(na,dp) * powa * powb * abs(f(na,nb,i,j))**2
                     sum3 = sum3 + real(nb,dp) * powa * powb * abs(f(na,nb,i,j))**2
                     sum4 = sum4 + real(na,dp)**2 * powa * powb * abs(f(na,nb,i,j))**2
                     sum5 = sum5 + real(nb,dp)**2 * powa * powb * abs(f(na,nb,i,j))**2
                     sum6 = sum6 + real(na,dp) * real(nb,dp) * powa * powb * abs(f(na,nb,i,j))**2
                     powa  = powa * x(1)
                  enddo
                  powb = powb / x(2)
               enddo

               c2d_GammaAB(1,1) = c2d_GammaAB(1,1) + sum2/sum1
               c2d_GammaAB(1,2) = c2d_GammaAB(1,2) + sum3/sum1
               c2d_GammaAB(2,1) = c2d_GammaAB(2,1) + (sum4 - sum2**2/sum1)/(sum1 * x(1))
               c2d_GammaAB(2,2) = c2d_GammaAB(2,2) + (sum6 - sum2*sum3/sum1)/(sum1 * x(2))
               c2d_GammaAB(3,1) = c2d_GammaAB(3,1) + (sum6 - sum2*sum3/sum1)/(sum1 * x(1))
               c2d_GammaAB(3,2) = c2d_GammaAB(3,2) + (sum5 - sum3**2/sum1)/(sum1 * x(2))

            enddo ! i
            enddo ! j

         endif

      ! 3) x > 1.0 and y <= 1.0
      elseif( abs(x(1)) > 1.0_dp .and. abs(x(2)) <= 1.0_dp ) then

         ! a) y = 0
         if( x(2) == 0.0_dp ) then

            do j = 1, ubound(f,4)
            do i = 1, ubound(f,3)
               sum1 = 0.0_dp
               sum2 = 0.0_dp
               sum3 = 0.0_dp
               sum4 = 0.0_dp
               sum5 = 0.0_dp
               powa = 1.0_dp
               do na = ubound(f,1), 0, -1
                  sum1 = sum1 + powa * abs(f(na,0,i,j))**2
                  sum2 = sum2 + real(na,dp) * powa * abs(f(na,0,i,j))**2
                  sum3 = sum3 + powa * abs(f(na,1,i,j))**2
                  sum4 = sum4 + real(na,dp)**2 * powa * abs(f(na,0,i,j))**2
                  sum5 = sum5 + real(na,dp) * powa * abs(f(na,1,i,j))**2
                  powa  = powa/x(1)
               enddo

               c2d_GammaAB(1,1) = c2d_GammaAB(1,1) + sum2/sum1
               c2d_GammaAB(2,1) = c2d_GammaAB(2,1) + (sum4 - sum2**2/sum1)/(sum1 * x(1))
               c2d_GammaAB(2,2) = c2d_GammaAB(2,2) + sum5/sum1
               c2d_GammaAB(3,2) = c2d_GammaAB(3,2) + sum3/sum1

            enddo ! i
            enddo ! j

         ! b) y != 0
         else

            do j = 1, ubound(f,4)
            do i = 1, ubound(f,3)
               sum1 = 0.0_dp
               sum2 = 0.0_dp
               sum3 = 0.0_dp
               sum4 = 0.0_dp
               sum5 = 0.0_dp
               sum6 = 0.0_dp
               powb = 1.0_dp
               do nb = 0, ubound(f,2)
                  powa = 1.0_dp
                  do na = ubound(f,1), 0, -1
                     sum1 = sum1 + powa * powb * abs(f(na,nb,i,j))**2
                     sum2 = sum2 + real(na,dp) * powa * powb * abs(f(na,nb,i,j))**2
                     sum3 = sum3 + real(nb,dp) * powa * powb * abs(f(na,nb,i,j))**2
                     sum4 = sum4 + real(na,dp)**2 * powa * powb * abs(f(na,nb,i,j))**2
                     sum5 = sum5 + real(nb,dp)**2 * powa * powb * abs(f(na,nb,i,j))**2
                     sum6 = sum6 + real(na,dp) * real(nb,dp) * powa * powb * abs(f(na,nb,i,j))**2
                     powa  = powa / x(1)
                  enddo
                  powb = powb * x(2)
               enddo

               c2d_GammaAB(1,1) = c2d_GammaAB(1,1) + sum2/sum1
               c2d_GammaAB(1,2) = c2d_GammaAB(1,2) + sum3/sum1
               c2d_GammaAB(2,1) = c2d_GammaAB(2,1) + (sum4 - sum2**2/sum1)/(sum1 * x(1))
               c2d_GammaAB(2,2) = c2d_GammaAB(2,2) + (sum6 - sum2*sum3/sum1)/(sum1 * x(2))
               c2d_GammaAB(3,1) = c2d_GammaAB(3,1) + (sum6 - sum2*sum3/sum1)/(sum1 * x(1))
               c2d_GammaAB(3,2) = c2d_GammaAB(3,2) + (sum5 - sum3**2/sum1)/(sum1 * x(2))

            enddo ! i
            enddo ! j

         endif

      ! 4) x > 1.0 and y > 1.0
      else 

         do j = 1, ubound(f,4)
         do i = 1, ubound(f,3)
            sum1 = 0.0_dp
            sum2 = 0.0_dp
            sum3 = 0.0_dp 
            sum4 = 0.0_dp
            sum5 = 0.0_dp
            sum6 = 0.0_dp
            powb = 1.0_dp
            do nb = ubound(f,2), 0, -1
               powa = 1.0_dp
               do na = ubound(f,1), 0, -1
                  sum1 = sum1 + powa * powb * abs(f(na,nb,i,j))**2
                  sum2 = sum2 + real(na,dp) * powa * powb * abs(f(na,nb,i,j))**2
                  sum3 = sum3 + real(nb,dp) * powa * powb * abs(f(na,nb,i,j))**2
                  sum4 = sum4 + real(na,dp)**2 * powa * powb * abs(f(na,nb,i,j))**2
                  sum5 = sum5 + real(nb,dp)**2 * powa * powb * abs(f(na,nb,i,j))**2
                  sum6 = sum6 + real(na,dp) * real(nb,dp) * powa * powb * abs(f(na,nb,i,j))**2
                  powa = powa / x(1)
               enddo
               powb = powb / x(2)
            enddo

            c2d_GammaAB(1,1) = c2d_GammaAB(1,1) + sum2/sum1
            c2d_GammaAB(1,2) = c2d_GammaAB(1,2) + sum3/sum1
            c2d_GammaAB(2,1) = c2d_GammaAB(2,1) + (sum4 - sum2**2/sum1)/(sum1 * x(1))
            c2d_GammaAB(2,2) = c2d_GammaAB(2,2) + (sum6 - sum2*sum3/sum1)/(sum1 * x(2))
            c2d_GammaAB(3,1) = c2d_GammaAB(3,1) + (sum6 - sum2*sum3/sum1)/(sum1 * x(1))
            c2d_GammaAB(3,2) = c2d_GammaAB(3,2) + (sum5 - sum3**2/sum1)/(sum1 * x(2))

         enddo ! i
         enddo ! j

      endif

     c2d_GammaAB(1,1) = c2d_GammaAB(1,1) - mna
     c2d_GammaAB(1,2) = c2d_GammaAB(1,2) - mnb

   end function

   function c3d_GammaAB( x, f, mna, mnb )

      implicit none
      real( kind = dp ), intent(inout) :: x(2)
      complex( kind = dp ), allocatable, intent(in) :: f(:,:,:,:,:)
      real( kind = dp ), intent(in) :: mna, mnb
      real( kind = dp), dimension(3,2) :: c3d_GammaAB

      ! local variables
      integer i, j, k, na, nb
      real( kind = dp ) sum1, sum2, sum3, sum4, sum5, sum6, powa, powb


      c3d_GammaAB = 0.0_dp

      ! separate definiton of functions into 4 different cases
      ! x(1) == x, x(2) == y
      ! 1) x <= 1.0 and y <= 1.0
      if( abs(x(1)) <= 1.0_dp .and. abs(x(2)) <= 1.0_dp ) then

         ! a) x == 0 and y != 0
         if( x(1) == 0.0_dp .and. x(2) /= 0.0_dp ) then

            do k = 1, ubound(f,5)
            do j = 1, ubound(f,4)
            do i = 1, ubound(f,3)
               sum1 = 0.0_dp
               sum2 = 0.0_dp
               sum3 = 0.0_dp
               sum4 = 0.0_dp
               sum5 = 0.0_dp
               powb = 1.0_dp
               do nb = 0, ubound(f,2)
                  sum1 = sum1 + powb * abs(f(0,nb,i,j,k))**2
                  sum2 = sum2 + real(nb,dp) * powb * abs(f(0,nb,i,j,k))**2
                  sum3 = sum3 + powb * abs(f(1,nb,i,j,k))**2
                  sum4 = sum4 + real(nb,dp)**2 * powb * abs(f(0,nb,i,j,k))**2
                  sum5 = sum5 + real(nb,dp) * powb * abs(f(1,nb,i,j,k))**2
                  powb = powb * x(2)
               enddo
               
               c3d_GammaAB(1,2) = c3d_GammaAB(1,2) + sum2/sum1
               c3d_GammaAB(2,1) = c3d_GammaAB(2,1) + sum3/sum1
               c3d_GammaAB(3,1) = c3d_GammaAB(3,1) + sum5/sum1
               c3d_GammaAB(3,2) = c3d_GammaAB(3,2) + (sum4 - sum2**2/sum1)/(sum1 * x(2))

            enddo ! i
            enddo ! j
            enddo ! k

         ! b) x != 0 and y == 0        
         elseif( x(1) /= 0.0_dp .and. x(2) == 0.0_dp ) then

            do k = 1, ubound(f,5)
            do j = 1, ubound(f,4)
            do i = 1, ubound(f,3)
               sum1 = 0.0_dp
               sum2 = 0.0_dp
               sum3 = 0.0_dp
               sum4 = 0.0_dp
               sum5 = 0.0_dp
               powa = 1.0_dp
               do na = 0, ubound(f,1)
                  sum1 = sum1 + powa * abs(f(na,0,i,j,k))**2
                  sum2 = sum2 + real(na,dp) * powa * abs(f(na,0,i,j,k))**2
                  sum3 = sum3 + powa * abs(f(na,1,i,j,k))**2
                  sum4 = sum4 + real(na,dp)**2 * powa * abs(f(na,0,i,j,k))**2
                  sum5 = sum5 + real(na,dp) * powa * abs(f(na,1,i,j,k))**2
                  powa  = powa * x(1)
               enddo

               c3d_GammaAB(1,1) = c3d_GammaAB(1,1) + sum2/sum1
               c3d_GammaAB(2,1) = c3d_GammaAB(2,1) + (sum4 - sum2**2/sum1)/(sum1 * abs(x(1)))
               c3d_GammaAB(2,2) = c3d_GammaAB(2,2) + sum5/sum1
               c3d_GammaAB(3,2) = c3d_GammaAB(3,2) + sum3/sum1

            enddo ! i
            enddo ! j
            enddo ! k

         ! c) x == 0 and y == 0
         elseif( x(1) == 0.0_dp .and. x(2) == 0.0_dp ) then

            do k = 1, ubound(f,5)
            do j = 1, ubound(f,4)
            do i = 1, ubound(f,3)
               c3d_GammaAB(2,1) = c3d_GammaAB(2,1) + abs(f(1,0,i,j,k))**2/abs(f(0,0,i,j,k))**2
               c3d_GammaAB(3,2) = c3d_GammaAB(3,2) + abs(f(0,1,i,j,k))**2/abs(f(0,0,i,j,k))**2
            enddo ! i
            enddo ! j
            enddo ! k

         ! d) x != 0 and y != 0
         else
           
            do k = 1, ubound(f,5)
            do j = 1, ubound(f,4)
            do i = 1, ubound(f,3)
               sum1 = 0.0_dp
               sum2 = 0.0_dp
               sum3 = 0.0_dp
               sum4 = 0.0_dp
               sum5 = 0.0_dp
               sum6 = 0.0_dp
               powb = 1.0_dp
               do nb = 0, ubound(f,2)
                  powa = 1.0_dp
                  do na = 0, ubound(f,1)
                     sum1 = sum1 + powa * powb * abs(f(na,nb,i,j,k))**2
                     sum2 = sum2 + real(na,dp) * powa * powb * abs(f(na,nb,i,j,k))**2
                     sum3 = sum3 + real(nb,dp) * powa * powb * abs(f(na,nb,i,j,k))**2
                     sum4 = sum4 + real(na,dp)**2 * powa * powb * abs(f(na,nb,i,j,k))**2
                     sum5 = sum5 + real(nb,dp)**2 * powa * powb * abs(f(na,nb,i,j,k))**2
                     sum6 = sum6 + real(na,dp) * real(nb,dp) * powa * powb * abs(f(na,nb,i,j,k))**2
                     powa = powa * x(1)
                  enddo
                  powb = powb * x(2)
               enddo

               c3d_GammaAB(1,1) = c3d_GammaAB(1,1) + sum2/sum1
               c3d_GammaAB(1,2) = c3d_GammaAB(1,2) + sum3/sum1
               c3d_GammaAB(2,1) = c3d_GammaAB(2,1) + (sum4 - sum2**2/sum1)/(sum1 * x(1))
               c3d_GammaAB(2,2) = c3d_GammaAB(2,2) + (sum6 - sum2*sum3/sum1)/(sum1 * x(2))
               c3d_GammaAB(3,1) = c3d_GammaAB(3,1) + (sum6 - sum2*sum3/sum1)/(sum1 * x(1))
               c3d_GammaAB(3,2) = c3d_GammaAB(3,2) + (sum5 - sum3**2/sum1)/(sum1 * x(2))

            enddo ! i
            enddo ! j
            enddo ! k
            
         endif 

      ! 2) x <= 1.0 and y > 1.0
      elseif( abs(x(1)) <= 1.0_dp .and. abs(x(2)) > 1.0_dp ) then

         ! a) x = 0
         if( x(1) == 0.0_dp ) then

            do k = 1, ubound(f,5)
            do j = 1, ubound(f,4)
            do i = 1, ubound(f,3)
               sum1 = 0.0_dp
               sum2 = 0.0_dp
               sum3 = 0.0_dp
               sum4 = 0.0_dp
               sum5 = 0.0_dp
               powb = 1.0_dp
               do nb = ubound(f,2), 0, -1
                  sum1 = sum1 + powb * abs(f(0,nb,i,j,k))**2
                  sum2 = sum2 + real(nb,dp) * powb * abs(f(0,nb,i,j,k))**2
                  sum3 = sum3 + powb * abs(f(1,nb,i,j,k))**2
                  sum4 = sum4 + real(nb,dp)**2 * powb * abs(f(0,nb,i,j,k))**2
                  sum5 = sum5 + real(nb,dp) * powb * abs(f(1,nb,i,j,k))**2
                  powb  = powb/x(2)
               enddo

               c3d_GammaAB(1,2) = c3d_GammaAB(1,2) + sum2/sum1
               c3d_GammaAB(2,1) = c3d_GammaAB(2,1) + sum3/sum1
               c3d_GammaAB(3,1) = c3d_GammaAB(3,1) + sum5/sum1
               c3d_GammaAB(3,2) = c3d_GammaAB(3,2) + (sum4 - sum2**2/sum1)/(sum1 * x(2))

            enddo ! i
            enddo ! j
            enddo ! k

         ! b) x != 0
         else

            do k = 1, ubound(f,5)
            do j = 1, ubound(f,4)
            do i = 1, ubound(f,3)
               sum1 = 0.0_dp
               sum2 = 0.0_dp
               sum3 = 0.0_dp
               sum4 = 0.0_dp
               sum5 = 0.0_dp
               sum6 = 0.0_dp
               powb = 1.0_dp
               do nb = ubound(f,2), 0, -1
                  powa = 1.0_dp
                  do na = 0, ubound(f,1)
                     sum1 = sum1 + powa * powb * abs(f(na,nb,i,j,k))**2
                     sum2 = sum2 + real(na,dp) * powa * powb * abs(f(na,nb,i,j,k))**2
                     sum3 = sum3 + real(nb,dp) * powa * powb * abs(f(na,nb,i,j,k))**2
                     sum4 = sum4 + real(na,dp)**2 * powa * powb * abs(f(na,nb,i,j,k))**2
                     sum5 = sum5 + real(nb,dp)**2 * powa * powb * abs(f(na,nb,i,j,k))**2
                     sum6 = sum6 + real(na,dp) * real(nb,dp) * powa * powb * abs(f(na,nb,i,j,k))**2
                     powa = powa * x(1)
                  enddo
                  powb = powb / x(2)
               enddo

               c3d_GammaAB(1,1) = c3d_GammaAB(1,1) + sum2/sum1
               c3d_GammaAB(1,2) = c3d_GammaAB(1,2) + sum3/sum1
               c3d_GammaAB(2,1) = c3d_GammaAB(2,1) + (sum4 - sum2**2/sum1)/(sum1 * x(1))
               c3d_GammaAB(2,2) = c3d_GammaAB(2,2) + (sum6 - sum2*sum3/sum1)/(sum1 * x(2))
               c3d_GammaAB(3,1) = c3d_GammaAB(3,1) + (sum6 - sum2*sum3/sum1)/(sum1 * x(1))
               c3d_GammaAB(3,2) = c3d_GammaAB(3,2) + (sum5 - sum3**2/sum1)/(sum1 * x(2))

            enddo ! i
            enddo ! j
            enddo ! k

         endif

      ! 3) x > 1.0 and y <= 1.0
      elseif( abs(x(1)) > 1.0_dp .and. abs(x(2)) <= 1.0_dp ) then

         ! a) y = 0
         if( x(2) == 0.0_dp ) then

            do k = 1, ubound(f,5)
            do j = 1, ubound(f,4)
            do i = 1, ubound(f,3)
               sum1 = 0.0_dp
               sum2 = 0.0_dp
               sum3 = 0.0_dp
               sum4 = 0.0_dp
               sum5 = 0.0_dp
               powa = 1.0_dp
               do na = ubound(f,1), 0, -1
                  sum1 = sum1 + powa * abs(f(na,0,i,j,k))**2
                  sum2 = sum2 + real(na,dp) * powa * abs(f(na,0,i,j,k))**2
                  sum3 = sum3 + powa * abs(f(na,1,i,j,k))**2
                  sum4 = sum4 + real(na,dp)**2 * powa * abs(f(na,0,i,j,k))**2
                  sum5 = sum5 + real(na,dp) * powa * abs(f(na,1,i,j,k))**2
                  powa  = powa/x(1)
               enddo

               c3d_GammaAB(1,1) = c3d_GammaAB(1,1) + sum2/sum1
               c3d_GammaAB(2,1) = c3d_GammaAB(2,1) + (sum4 - sum2**2/sum1)/(sum1 * x(1))
               c3d_GammaAB(2,2) = c3d_GammaAB(2,2) + sum5/sum1
               c3d_GammaAB(3,2) = c3d_GammaAB(3,2) + sum3/sum1

            enddo ! i
            enddo ! j
            enddo ! k

         ! b) y != 0
         else

            do k = 1, ubound(f,5)
            do j = 1, ubound(f,4)
            do i = 1, ubound(f,3)
               sum1 = 0.0_dp
               sum2 = 0.0_dp
               sum3 = 0.0_dp
               sum4 = 0.0_dp
               sum5 = 0.0_dp
               sum6 = 0.0_dp
               powb = 1.0_dp
               do nb = 0, ubound(f,2)
                  powa = 1.0_dp
                  do na = ubound(f,1), 0, -1
                     sum1 = sum1 + powa * powb * abs(f(na,nb,i,j,k))**2
                     sum2 = sum2 + real(na,dp) * powa * powb * abs(f(na,nb,i,j,k))**2
                     sum3 = sum3 + real(nb,dp) * powa * powb * abs(f(na,nb,i,j,k))**2
                     sum4 = sum4 + real(na,dp)**2 * powa * powb * abs(f(na,nb,i,j,k))**2
                     sum5 = sum5 + real(nb,dp)**2 * powa * powb * abs(f(na,nb,i,j,k))**2
                     sum6 = sum6 + real(na,dp) * real(nb,dp) * powa * powb * abs(f(na,nb,i,j,k))**2
                     powa = powa / x(1)
                  enddo
                  powb = powb * x(2)
               enddo

               c3d_GammaAB(1,1) = c3d_GammaAB(1,1) + sum2/sum1
               c3d_GammaAB(1,2) = c3d_GammaAB(1,2) + sum3/sum1
               c3d_GammaAB(2,1) = c3d_GammaAB(2,1) + (sum4 - sum2**2/sum1)/(sum1 * x(1))
               c3d_GammaAB(2,2) = c3d_GammaAB(2,2) + (sum6 - sum2*sum3/sum1)/(sum1 * x(2))
               c3d_GammaAB(3,1) = c3d_GammaAB(3,1) + (sum6 - sum2*sum3/sum1)/(sum1 * x(1))
               c3d_GammaAB(3,2) = c3d_GammaAB(3,2) + (sum5 - sum3**2/sum1)/(sum1 * x(2))

            enddo ! i
            enddo ! j
            enddo ! k

         endif

      ! 4) x > 1.0 and y > 1.0
      else 

         do k = 1, ubound(f,5)
         do j = 1, ubound(f,4)
         do i = 1, ubound(f,3)
            sum1 = 0.0_dp
            sum2 = 0.0_dp
            sum3 = 0.0_dp 
            sum4 = 0.0_dp
            sum5 = 0.0_dp
            sum6 = 0.0_dp
            powb = 1.0_dp
            do nb = ubound(f,2), 0, -1
               powa = 1.0_dp
               do na = ubound(f,1), 0, -1
                  sum1 = sum1 + powa * powb * abs(f(na,nb,i,j,k))**2
                  sum2 = sum2 + real(na,dp) * powa * powb * abs(f(na,nb,i,j,k))**2
                  sum3 = sum3 + real(nb,dp) * powa * powb * abs(f(na,nb,i,j,k))**2
                  sum4 = sum4 + real(na,dp)**2 * powa * powb * abs(f(na,nb,i,j,k))**2
                  sum5 = sum5 + real(nb,dp)**2 * powa * powb * abs(f(na,nb,i,j,k))**2
                  sum6 = sum6 + real(na,dp) * real(nb,dp) * powa * powb * abs(f(na,nb,i,j,k))**2
                  powa = powa / x(1)
               enddo
               powb = powb / x(2)
            enddo

            c3d_GammaAB(1,1) = c3d_GammaAB(1,1) + sum2/sum1
            c3d_GammaAB(1,2) = c3d_GammaAB(1,2) + sum3/sum1
            c3d_GammaAB(2,1) = c3d_GammaAB(2,1) + (sum4 - sum2**2/sum1)/(sum1 * x(1))
            c3d_GammaAB(2,2) = c3d_GammaAB(2,2) + (sum6 - sum2*sum3/sum1)/(sum1 * x(2))
            c3d_GammaAB(3,1) = c3d_GammaAB(3,1) + (sum6 - sum2*sum3/sum1)/(sum1 * x(1))
            c3d_GammaAB(3,2) = c3d_GammaAB(3,2) + (sum5 - sum3**2/sum1)/(sum1 * x(2))

         enddo ! i
         enddo ! j
         enddo ! k

      endif

     c3d_GammaAB(1,1) = c3d_GammaAB(1,1) - mna
     c3d_GammaAB(1,2) = c3d_GammaAB(1,2) - mnb

   end function

! | ---------------------------------------- |
! | scalar function nedded for nelder-mead   |
! | minimization with penalty outside domain |
! | of interest                              |
! | ---------------------------------------- |
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
     
            do i = 1, ubound(f,3)
               sum1 = 0.0_dp
               sum2 = 0.0_dp
               sum3 = 0.0_dp
               powb = 1.0_dp
               do nb = 0, ubound(f,2)
                  powa = 1.0_dp
                  do na = 0, ubound(f,1)
                     sum1 = sum1 + powa * powb * abs(f(na,nb,i))**2
                     sum2 = sum2 + powa * powb * real(na,dp)*abs(f(na,nb,i))**2
                     sum3 = sum3 + powa * powb * real(nb,dp)*abs(f(na,nb,i))**2
                     powa = powa * x(1)
                  enddo
                  powb = powb * x(2)
               enddo

               fx(1) = fx(1) + sum2/sum1
               fx(2) = fx(2) + sum3/sum1

            enddo

         elseif( abs(x(1)) <= 1.0_dp .and. abs(x(2)) > 1.0_dp ) then

            do i = 1, ubound(f,3)
               sum1 = 0.0_dp
               sum2 = 0.0_dp
               sum3 = 0.0_dp
               powb = 1.0_dp
               do nb = ubound(f,2), 0, -1
                  powa = 1.0_dp
                  do na = 0, ubound(f,1)
                     sum1 = sum1 + powa * powb * abs(f(na,nb,i))**2
                     sum2 = sum2 + powa * powb * real(na,dp)*abs(f(na,nb,i))**2
                     sum3 = sum3 + powa * powb * real(nb,dp)*abs(f(na,nb,i))**2
                     powa = powa * x(1)
                  enddo
                  powb = powb / x(2)
               enddo

               fx(1) = fx(1) + sum2/sum1
               fx(2) = fx(2) + sum3/sum1

            enddo

         elseif( abs(x(1)) > 1.0_dp .and. abs(x(2)) <= 1.0_dp ) then

            do i = 1, ubound(f,3)
               sum1 = 0.0_dp
               sum2 = 0.0_dp
               sum3 = 0.0_dp
               powb = 1.0_dp
               do nb = 0, ubound(f,2)
                  powa = 1.0_dp
                  do na = ubound(f,1), 0, -1
                     sum1 = sum1 + powa * powb * abs(f(na,nb,i))**2
                     sum2 = sum2 + powa * powb * real(na,dp)*abs(f(na,nb,i))**2
                     sum3 = sum3 + powa * powb * real(nb,dp)*abs(f(na,nb,i))**2
                     powa = powa / x(1)
                  enddo
                  powb = powb * x(2)
               enddo

               fx(1) = fx(1) + sum2/sum1
               fx(2) = fx(2) + sum3/sum1

            enddo

         else

            do i = 1, ubound(f,3)
               sum1 = 0.0_dp
               sum2 = 0.0_dp
               sum3 = 0.0_dp
               powb = 1.0_dp
               do nb = ubound(f,2), 0, -1
                  powa = 1.0_dp
                  do na = ubound(f,1), 0, -1
                     sum1 = sum1 + powa * powb * abs(f(na,nb,i))**2
                     sum2 = sum2 + powa * powb * real(na,dp)*abs(f(na,nb,i))**2
                     sum3 = sum3 + powa * powb * real(nb,dp)*abs(f(na,nb,i))**2
                     powa = powa / x(1)
                  enddo
                  powb = powb / x(2)
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

   function c2d_nelder_mead( x, f, mna, mnb )

      implicit none
      real( kind = dp ), intent(in), dimension(2) :: x
      complex( kind = dp ), intent(in), allocatable :: f(:,:,:,:)
      real( kind = dp ), intent(in) :: mna, mnb
      real( kind = dp ) :: c2d_nelder_mead

      real( kind = dp ), parameter :: penalty = 1.0E+14
      real( kind = dp ) :: fx(2)
      real( kind = dp ) sum1, sum2, sum3, powa, powb
      integer i, j, na, nb

      if( x(1) < 0.0_dp .or. x(2) < 0.0_dp .or. x(1) > 2.0_dp .or. x(2) > 2.0_dp ) then
         c2d_nelder_mead = penalty
      else

         fx = 0.0_dp

         if( abs(x(1)) <= 1.0_dp .and. abs(x(2)) <= 1.0_dp ) then
     
            do j = 1, ubound(f,4)
            do i = 1, ubound(f,3)
               sum1 = 0.0_dp
               sum2 = 0.0_dp
               sum3 = 0.0_dp
               powb = 1.0_dp
               do nb = 0, ubound(f,2)
                  powa = 1.0_dp
                  do na = 0, ubound(f,1)
                     sum1 = sum1 + powa * powb * abs(f(na,nb,i,j))**2
                     sum2 = sum2 + powa * powb * real(na,dp)*abs(f(na,nb,i,j))**2
                     sum3 = sum3 + powa * powb * real(nb,dp)*abs(f(na,nb,i,j))**2
                     powa = powa * x(1)
                  enddo
                  powb = powb * x(2)
               enddo

               fx(1) = fx(1) + sum2/sum1
               fx(2) = fx(2) + sum3/sum1

            enddo ! i
            enddo ! j

         elseif( abs(x(1)) <= 1.0_dp .and. abs(x(2)) > 1.0_dp ) then

            do j = 1, ubound(f,4)
            do i = 1, ubound(f,3)
               sum1 = 0.0_dp
               sum2 = 0.0_dp
               sum3 = 0.0_dp
               powb = 1.0_dp
               do nb = ubound(f,2), 0, -1
                  powa = 1.0_dp
                  do na = 0, ubound(f,1)
                     sum1 = sum1 + powa * powb * abs(f(na,nb,i,j))**2
                     sum2 = sum2 + powa * powb * real(na,dp)*abs(f(na,nb,i,j))**2
                     sum3 = sum3 + powa * powb * real(nb,dp)*abs(f(na,nb,i,j))**2
                     powa = powa * x(1)
                  enddo
                  powb = powb / x(2)
               enddo

               fx(1) = fx(1) + sum2/sum1
               fx(2) = fx(2) + sum3/sum1

            enddo ! i
            enddo ! j

         elseif( abs(x(1)) > 1.0_dp .and. abs(x(2)) <= 1.0_dp ) then

            do j = 1, ubound(f,4)
            do i = 1, ubound(f,3)
               sum1 = 0.0_dp
               sum2 = 0.0_dp
               sum3 = 0.0_dp
               powb = 1.0_dp
               do nb = 0, ubound(f,2)
                  powa = 1.0_dp
                  do na = ubound(f,1), 0, -1
                     sum1 = sum1 + powa * powb * abs(f(na,nb,i,j))**2
                     sum2 = sum2 + powa * powb * real(na,dp)*abs(f(na,nb,i,j))**2
                     sum3 = sum3 + powa * powb * real(nb,dp)*abs(f(na,nb,i,j))**2
                     powa = powa / x(1)
                  enddo
                  powb = powb * x(2)
               enddo

               fx(1) = fx(1) + sum2/sum1
               fx(2) = fx(2) + sum3/sum1

            enddo ! i
            enddo ! j

         else

            do j = 1, ubound(f,4)
            do i = 1, ubound(f,3)
               sum1 = 0.0_dp
               sum2 = 0.0_dp
               sum3 = 0.0_dp
               powb = 1.0_dp
               do nb = ubound(f,2), 0, -1
                  powa = 1.0_dp
                  do na = ubound(f,1), 0, -1
                     sum1 = sum1 + powa * powb * abs(f(na,nb,i,j))**2
                     sum2 = sum2 + powa * powb * real(na,dp)*abs(f(na,nb,i,j))**2
                     sum3 = sum3 + powa * powb * real(nb,dp)*abs(f(na,nb,i,j))**2
                     powa = powa / x(1)
                  enddo
                  powb = powb / x(2)
               enddo

               fx(1) = fx(1) + sum2/sum1
               fx(2) = fx(2) + sum3/sum1

            enddo ! i
            enddo ! j

         endif
 
         fx(1) = fx(1) - mna
         fx(2) = fx(2) - mnb

         c2d_nelder_mead = LOG10( fx(1)**2 + fx(2)**2 + 1.0_dp )

      endif
         
   end function

   function c3d_nelder_mead( x, f, mna, mnb )

      implicit none
      real( kind = dp ), intent(in), dimension(2) :: x
      complex( kind = dp ), intent(in), allocatable :: f(:,:,:,:,:)
      real( kind = dp ), intent(in) :: mna, mnb
      real( kind = dp ) :: c3d_nelder_mead

      real( kind = dp ), parameter :: penalty = 1.0E+14
      real( kind = dp ) :: fx(2)
      real( kind = dp ) sum1, sum2, sum3, powa, powb
      integer i, j, k,  na, nb

      if( x(1) < 0.0_dp .or. x(2) < 0.0_dp .or. x(1) > 2.0_dp .or. x(2) > 2.0_dp ) then
         c3d_nelder_mead = penalty
      else

         fx = 0.0_dp

         if( abs(x(1)) <= 1.0_dp .and. abs(x(2)) <= 1.0_dp ) then
     
            do k = 1, ubound(f,5)
            do j = 1, ubound(f,4)
            do i = 1, ubound(f,3)
               sum1 = 0.0_dp
               sum2 = 0.0_dp
               sum3 = 0.0_dp
               powb = 1.0_dp
               do nb = 0, ubound(f,2)
                  powa = 1.0_dp
                  do na = 0, ubound(f,1)
                     sum1 = sum1 + powa * powb * abs(f(na,nb,i,j,k))**2
                     sum2 = sum2 + powa * powb * real(na,dp)*abs(f(na,nb,i,j,k))**2
                     sum3 = sum3 + powa * powb * real(nb,dp)*abs(f(na,nb,i,j,k))**2
                     powa = powa * x(1)
                  enddo
                  powb = powb * x(2)
               enddo

               fx(1) = fx(1) + sum2/sum1
               fx(2) = fx(2) + sum3/sum1

            enddo ! i
            enddo ! j
            enddo ! k

         elseif( abs(x(1)) <= 1.0_dp .and. abs(x(2)) > 1.0_dp ) then

            do k = 1, ubound(f,5)
            do j = 1, ubound(f,4)
            do i = 1, ubound(f,3)
               sum1 = 0.0_dp
               sum2 = 0.0_dp
               sum3 = 0.0_dp
               powb = 1.0_dp
               do nb = ubound(f,2), 0, -1
                  powa = 1.0_dp
                  do na = 0, ubound(f,1)
                     sum1 = sum1 + powa * powb * abs(f(na,nb,i,j,k))**2
                     sum2 = sum2 + powa * powb * real(na,dp)*abs(f(na,nb,i,j,k))**2
                     sum3 = sum3 + powa * powb * real(nb,dp)*abs(f(na,nb,i,j,k))**2
                     powa = powa * x(1)
                  enddo
                  powb = powb / x(2)
               enddo

               fx(1) = fx(1) + sum2/sum1
               fx(2) = fx(2) + sum3/sum1

            enddo ! i
            enddo ! j
            enddo ! k

         elseif( abs(x(1)) > 1.0_dp .and. abs(x(2)) <= 1.0_dp ) then

            do k = 1, ubound(f,5)
            do j = 1, ubound(f,4)
            do i = 1, ubound(f,3)
               sum1 = 0.0_dp
               sum2 = 0.0_dp
               sum3 = 0.0_dp
               powb = 1.0_dp
               do nb = 0, ubound(f,2)
                  powa = 1.0_dp
                  do na = ubound(f,1), 0, -1
                     sum1 = sum1 + powa * powb * abs(f(na,nb,i,j,k))**2
                     sum2 = sum2 + powa * powb * real(na,dp)*abs(f(na,nb,i,j,k))**2
                     sum3 = sum3 + powa * powb * real(nb,dp)*abs(f(na,nb,i,j,k))**2
                     powa = powa / x(1)
                  enddo
                  powb = powb * x(2)
               enddo

               fx(1) = fx(1) + sum2/sum1
               fx(2) = fx(2) + sum3/sum1

            enddo ! i
            enddo ! j
            enddo ! k

         else

            do k = 1, ubound(f,5)
            do j = 1, ubound(f,4)
            do i = 1, ubound(f,3)
               sum1 = 0.0_dp
               sum2 = 0.0_dp
               sum3 = 0.0_dp
               powb = 1.0_dp
               do nb = ubound(f,2), 0, -1
                  powa = 1.0_dp
                  do na = ubound(f,1), 0, -1
                     sum1 = sum1 + powa * powb * abs(f(na,nb,i,j,k))**2
                     sum2 = sum2 + powa * powb * real(na,dp)*abs(f(na,nb,i,j,k))**2
                     sum3 = sum3 + powa * powb * real(nb,dp)*abs(f(na,nb,i,j,k))**2
                     powa = powa / x(1)
                  enddo
                  powb = powb / x(2)
               enddo

               fx(1) = fx(1) + sum2/sum1
               fx(2) = fx(2) + sum3/sum1

            enddo ! i
            enddo ! j
            enddo ! k

         endif
 
         fx(1) = fx(1) - mna
         fx(2) = fx(2) - mnb

         c3d_nelder_mead = LOG10( fx(1)**2 + fx(2)**2 + 1.0_dp )

      endif
         
   end function

! | ---------------------------------------- |
! | 1d root finding algortihm                |
! | It combines pre-bisection method and     |
! | then refines with newton algorithm       |
! | ---------------------------------------- | 
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

   subroutine c2d_root1d( x, f, mn, fun )

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

   subroutine c3d_root1d( x, f, mn, fun )

      implicit none
      real( kind = dp ), intent(inout) :: x
      complex( kind = dp ), allocatable, intent(in) :: f(:,:,:,:,:)
      real( kind = dp ), intent(in) :: mn

      interface
         function fun( x, f, mn )
            import
            real( kind = dp ), intent(in) :: x
            complex( kind = dp ), allocatable, intent(in) :: f(:,:,:,:,:)
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

! | ---------------------------------------- |
! | 2d root finding algorithm                |
! | combination of initial contrained        |
! | nelder-mead minimization and newton      |
! | ---------------------------------------- | 
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
      call c1d_nelmin( c1d_nelder_mead, x, f, mna, mnb )

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

   subroutine c2d_root2d( x, f, mna, mnb )

      implicit none
      real( kind = dp ), intent(inout) :: x(2)
      complex( kind = dp ), allocatable, intent(in) :: f(:,:,:,:)
      real( kind = dp ), intent(in) :: mna, mnb

      ! local parameters
      real( kind = dp ), dimension(3,2) :: fk
      real( kind = dp ), dimension(2,2) :: invj
      real( kind = dp ) det, dx, dy
      real( kind = dp ), parameter :: conv = 1.0E-15
      integer, parameter :: maxiter = 50
      integer iter

      ! start with nelder-mead
      call c2d_nelmin( c2d_nelder_mead, x, f, mna, mnb )

      ! refine with Newtons's algorithm
      iter = 0

      do while( iter < maxiter )

         iter = iter + 1

         fk = c2d_GammaAB( x, f, mna, mnb )

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

   subroutine c3d_root2d( x, f, mna, mnb )

      implicit none
      real( kind = dp ), intent(inout) :: x(2)
      complex( kind = dp ), allocatable, intent(in) :: f(:,:,:,:,:)
      real( kind = dp ), intent(in) :: mna, mnb

      ! local parameters
      real( kind = dp ), dimension(3,2) :: fk
      real( kind = dp ), dimension(2,2) :: invj
      real( kind = dp ) det, dx, dy
      real( kind = dp ), parameter :: conv = 1.0E-15
      integer, parameter :: maxiter = 50
      integer iter

      ! start with nelder-mead
      call c3d_nelmin( c3d_nelder_mead, x, f, mna, mnb )

      ! refine with Newtons's algorithm
      iter = 0

      do while( iter < maxiter )

         iter = iter + 1

         fk = c3d_GammaAB( x, f, mna, mnb )

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

! | ---------------------------------------- |
! | Nelder-Mead mimimization algorithm       |
! | ---------------------------------------- | 
   subroutine c1d_nelmin ( fn, xmin, f, mna, mnb )

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

   subroutine c2d_nelmin ( fn, xmin, f, mna, mnb )

     implicit none

     complex( kind = dp ), intent(in), allocatable :: f(:,:,:,:)
     real( kind = dp ), intent(in) :: mna, mnb
     interface 
        function fn( x, f, mna, mnb )
           import
           real( kind = dp ), intent(in) :: x(2)
           complex( kind = dp ), allocatable, intent(in) :: f(:,:,:,:)
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

   subroutine c3d_nelmin ( fn, xmin, f, mna, mnb )

     implicit none

     complex( kind = dp ), intent(in), allocatable :: f(:,:,:,:,:)
     real( kind = dp ), intent(in) :: mna, mnb
     interface 
        function fn( x, f, mna, mnb )
           import
           real( kind = dp ), intent(in) :: x(2)
           complex( kind = dp ), allocatable, intent(in) :: f(:,:,:,:,:)
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

      do i = 1, ubound(f,3)
         powb = 1.0_dp
         tmp = 0.0_dp
         do nb = 0, ubound(f,2)
            powa = 1.0_dp
            do na = 0, ubound(f,1)
               tmp = tmp + powa * powb * abs(f(na,nb,i))**2
               powa = powa * sig(1)
            enddo
            powb = powb * sig(2)
         enddo
         vout(i) = 1.0_dp/tmp
      enddo

   end subroutine

   subroutine c2d_sigma_norm( f, sig, vout)

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:,:,:)
      real( kind = dp ), intent(in), dimension(2) :: sig
      real( kind = dp ), intent(inout), allocatable :: vout(:,:)

      integer i, j, na, nb
      real( kind = dp ) :: powa, powb, tmp

      do j = 1, ubound(f,4)
         do i = 1, ubound(f,3)
            powb = 1.0_dp
            tmp = 0.0_dp
            do nb = 0, ubound(f,2)
               powa = 1.0_dp
               do na = 0, ubound(f,1)
                  tmp = tmp + powa * powb * abs(f(na,nb,i,j))**2
                  powa = powa * sig(1)
               enddo
               powb = powb * sig(2)
            enddo
            vout(i,j) = 1.0_dp/tmp
         enddo
      enddo

   end subroutine

   subroutine c3d_sigma_norm( f, sig, vout)

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:,:,:,:)
      real( kind = dp ), intent(in), dimension(2) :: sig
      real( kind = dp ), intent(inout), allocatable :: vout(:,:,:)

      integer i, j, k, na, nb
      real( kind = dp ) :: powa, powb, tmp

      do k = 1, ubound(f,5)
         do  j = 1, ubound(f,4)
            do i = 1, ubound(f,3)
               powb = 1.0_dp
               tmp = 0.0_dp
               do nb = 0, ubound(f,2)
                  powa = 1.0_dp
                  do na = 0, ubound(f,1)
                     tmp = tmp + powa * powb * abs(f(na,nb,i,j,k))**2
                     powa = powa * sig(1)
                  enddo
                  powb = powb * sig(2)
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
   subroutine c1d_const_mean( f, mna, mnb, info )

      implicit none
      complex( kind = dp ), intent(inout), allocatable :: f(:,:,:)
      real( kind = dp ) :: mna, mnb
      integer :: info

      real( kind = dp ) :: powa, powb
      real( kind = dp ), dimension(2) :: mu
      integer i, na, nb
      real( kind = dp ), allocatable :: norm(:)

      allocate( norm(lbound(f,3): ubound(f,3)) )
      
      info = 0

      if( mna == 0.0_dp ) then

          mu(2) = 1.0_dp
          call c1d_root1d( mu(2), f, mnb, c1d_GammaB )
          mu(1) = 0.0_dp
          call c1d_sigma_norm( f, mu, norm )

      elseif( mnb == 0.0_dp ) then
          
          mu(1) = 1.0_dp
          call c1d_root1d( mu(1), f, mna, c1d_GammaA )
          mu(2) = 0.0_dp
          call c1d_sigma_norm( f, mu, norm )

      else

         call c1d_root2d( mu, f, mna, mnb )
         call c1d_sigma_norm( f, mu, norm )

      endif

      if( mu(1) < 0.0_dp  .or. mu(2) < 0.0_dp ) info = 1
      if( mu(1) == 1.0_dp .and. mu(2) == 1.0_dp ) info = 2
      if( ISNAN(mu(1)) .or. ISNAN(mu(2)) ) info = 3
 

      if( info  == 0 ) then
         do i = 1, ubound(f,3)
            powb = 1.0_dp
            do nb = 0, ubound(f,2)
               powa = 1.0_dp
               do na = 0, ubound(f,1)
                  f(na,nb,i) = f(na,nb,i) * sqrt(norm(i)) * powa * powb
                  powa = powa * sqrt(mu(1))
               enddo
               powb = powb * sqrt(mu(2))
            enddo
         enddo
      endif

   end subroutine

   subroutine c2d_const_mean( f, mna, mnb, info )

      implicit none
      complex( kind = dp ), intent(inout), allocatable :: f(:,:,:,:)
      real( kind = dp ) :: mna, mnb
      integer :: info

      real( kind = dp ) :: powa, powb
      real( kind = dp ), dimension(2) :: mu
      integer i, j, na, nb
      real( kind = dp ), allocatable :: norm(:,:)

      allocate( norm(lbound(f,3): ubound(f,3), lbound(f,4):ubound(f,4)) )
      
      info = 0

      if( mna == 0.0_dp ) then

          mu(2) = 1.0_dp
          call c2d_root1d( mu(2), f, mnb, c2d_GammaB )
          mu(1) = 0.0_dp
          call c2d_sigma_norm( f, mu, norm )

      elseif( mnb == 0.0_dp ) then
          
          mu(1) = 1.0_dp
          call c2d_root1d( mu(1), f, mna, c2d_GammaA )
          mu(2) = 0.0_dp
          call c2d_sigma_norm( f, mu, norm )

      else

         call c2d_root2d( mu, f, mna, mnb )
         call c2d_sigma_norm( f, mu, norm )

      endif

      if( mu(1) < 0.0_dp  .or. mu(2) < 0.0_dp ) info = 1
      if( mu(1) == 1.0_dp .and. mu(2) == 1.0_dp ) info = 2
      if( ISNAN(mu(1)) .or. ISNAN(mu(2)) ) info = 3
 

      if( info  == 0 ) then
         do j = 1, ubound(f,4)
         do i = 1, ubound(f,3)
            powb = 1.0_dp
            do nb = 0, ubound(f,2)
               powa = 1.0_dp
               do na = 0, ubound(f,1)
                  f(na,nb,i,j) = f(na,nb,i,j) * sqrt(norm(i,j)) * powa * powb
                  powa = powa * sqrt(mu(1))
               enddo
               powb = powb * sqrt(mu(2))
            enddo
         enddo ! i
         enddo ! j
      endif

   end subroutine

   subroutine c3d_const_mean( f, mna, mnb, info )

      implicit none
      complex( kind = dp ), intent(inout), allocatable :: f(:,:,:,:,:)
      real( kind = dp ) :: mna, mnb
      integer :: info

      real( kind = dp ) :: powa, powb
      real( kind = dp ), dimension(2) :: mu
      integer i, j, k, na, nb
      real( kind = dp ), allocatable :: norm(:,:,:)

      allocate( norm(lbound(f,3): ubound(f,3), lbound(f,4):ubound(f,4), lbound(f,5):ubound(f,5)) )
      
      info = 0

      if( mna == 0.0_dp ) then

          mu(2) = 1.0_dp
          call c3d_root1d( mu(2), f, mnb, c3d_GammaB )
          mu(1) = 0.0_dp
          call c3d_sigma_norm( f, mu, norm )

      elseif( mnb == 0.0_dp ) then
          
          mu(1) = 1.0_dp
          call c3d_root1d( mu(1), f, mna, c3d_GammaA )
          mu(2) = 0.0_dp
          call c3d_sigma_norm( f, mu, norm )

      else

         call c3d_root2d( mu, f, mna, mnb )
         call c3d_sigma_norm( f, mu, norm )

      endif

      if( mu(1) < 0.0_dp  .or. mu(2) < 0.0_dp ) info = 1
      if( mu(1) == 1.0_dp .and. mu(2) == 1.0_dp ) info = 2
      if( ISNAN(mu(1)) .or. ISNAN(mu(2)) ) info = 3
 

      if( info  == 0 ) then
         do k = 1, ubound(f,5)
         do j = 1, ubound(f,4)
         do i = 1, ubound(f,3)
            powb = 1.0_dp
            do nb = 0, ubound(f,2)
               powa = 1.0_dp
               do na = 0, ubound(f,1)
                  f(na,nb,i,j,k) = f(na,nb,i,j,k) * sqrt(norm(i,j,k)) * powa * powb
                  powa = powa * sqrt(mu(1))
               enddo
               powb = powb * sqrt(mu(2))
            enddo
         enddo ! i
         enddo ! j
         enddo ! k
      endif

   end subroutine

end module
