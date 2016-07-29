module measures

   use parameters, only: dp, re, im

   private
   public Mean, Var, norm, normalize, Order, TotEnergy

   ! | ---------------------------------------- |
   ! | For better memeory layout it is better   |
   ! | to put Fock basis index at the beginning |
   ! | ---------------------------------------- |
   interface Mean
      module procedure c1d_Mean
      module procedure c2d_Mean
      module procedure c3d_Mean
   end interface

   interface Var
      module procedure c1d_Var
      module procedure c2d_Var
      module procedure c3d_Var
   end interface 

   interface norm
      module procedure c1d_norm
      module procedure c2d_norm
      module procedure c3d_norm
   end interface

   interface normalize
      module procedure c1d_normalize
      module procedure c2d_normalize
      module procedure c3d_normalize
   end interface

   interface Order
      module procedure c1d_Order
      module procedure c2d_Order
      module procedure c3d_Order
   end interface

   interface TotEnergy
      module procedure c1d_TotEnergy
      module procedure c2d_TotEnergy
      module procedure c3d_TotEnergy
   end interface  

   interface ChemPot
      module procedure c1d_ChemPot
      module procedure c2d_ChemPot
      module procedure c3d_ChemPot
   end interface

contains

!  | ----------------------------------- |
!  | Mean number of particles in comp. A |
!  | ----------------------------------- |
   function c1d_Mean( f, i )

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:)
      integer, intent(in) :: i
      real( kind = dp ) :: c1d_Mean

      ! local variables
      integer  na

      c1d_Mean = 0.0_dp
      do na = 0, ubound(f,1)
         c1d_Mean = c1d_Mean + real( na,dp ) * abs(f(na,i))**2
      enddo

   end function

   function c2d_Mean( f, i, j )

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:,:)
      integer, intent(in) :: i, j
      real( kind = dp ) :: c2d_Mean

      ! local variables
      integer  na

      c2d_Mean = 0.0_dp
      do na = 0, ubound(f,1)
         c2d_Mean = c2d_Mean + real( na,dp ) * abs(f(na,i,j))**2
      enddo

   end function

   function c3d_Mean( f, i, j, k )

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:,:,:)
      integer, intent(in) :: i, j, k
      real( kind = dp ) :: c3d_Mean

      ! local variables
      integer  na

      c3d_Mean = 0.0_dp
      do na = 0, ubound(f,1)
         c3d_Mean = c3d_Mean + real( na,dp ) * abs(f(na,i,j,k))**2
      enddo

   end function

!  | ----------------------------------- |
!  | Variance of particle number in A    |
!  | ----------------------------------- |
   function c1d_Var( f, i )

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:)
      integer, intent(in) :: i
      real( kind = dp ) :: c1d_Var

      ! local variables
      integer  na
      real( kind = dp ) suma

      c1d_Var = 0.0_dp
      suma   = 0.0_dp

      do na = 0, ubound(f,1)
         suma   = suma + real( na,dp ) * abs(f(na,i))**2
         c1d_Var = c1d_Var + real( na*na,dp ) * abs(f(na,i))**2
      enddo

      c1d_Var = c1d_Var - suma**2

   end function

   function c2d_Var( f, i, j )

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:,:)
      integer, intent(in) :: i, j
      real( kind = dp ) :: c2d_Var

      ! local variables
      integer  na
      real( kind = dp ) suma

      c2d_Var = 0.0_dp
      suma   = 0.0_dp

      do na = 0, ubound(f,1)
         suma   = suma + real( na,dp ) * abs(f(na,i,j))**2
         c2d_Var = c2d_Var + real( na*na,dp ) * abs(f(na,i,j))**2
      enddo

      c2d_Var = c2d_Var - suma**2

   end function

   function c3d_Var( f, i, j, k )

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:,:,:)
      integer, intent(in) :: i, j, k
      real( kind = dp ) :: c3d_Var

      ! local variables
      integer  na
      real( kind = dp ) suma

      c3d_Var = 0.0_dp
      suma   = 0.0_dp

      do na = 0, ubound(f,1)
         suma   = suma + real( na,dp ) * abs(f(na,i,j,k))**2
         c3d_Var = c3d_Var + real( na*na,dp ) * abs(f(na,i,j,k))**2
      enddo

      c3d_Var = c3d_Var - suma**2

   end function

!  | ----------------------------------- |
!  | Local norm of Gutzwiller wavefun    |
!  | ----------------------------------- |
   function c1d_norm( f, i )

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:)
      integer, intent(in) :: i
      real( kind = dp ) :: c1d_norm

      ! local variables
      integer  na

      c1d_norm = 0.0_dp

      do na = 0, ubound(f,1)
         c1d_norm = c1d_norm + abs(f(na,i))**2
      enddo

   end function

   function c2d_norm( f, i, j )

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:,:)
      integer, intent(in) :: i, j
      real( kind = dp ) :: c2d_norm

      ! local variables
      integer  na

      c2d_norm = 0.0_dp

      do na = 0, ubound(f,1)
         c2d_norm = c2d_norm + abs(f(na,i,j))**2
      enddo

   end function

   function c3d_norm( f, i, j, k )

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:,:,:)
      integer, intent(in) :: i, j, k
      real( kind = dp ) :: c3d_norm

      ! local variables
      integer  na

      c3d_norm = 0.0_dp

      do na = 0, ubound(f,1)
         c3d_norm = c3d_norm + abs(f(na,i,j,k))**2
      enddo

   end function

!  | ----------------------------------- |
!  | normalize Gutzwiller wavefun        |
!  | ----------------------------------- |
   subroutine c1d_normalize( f )

      implicit none
      complex( kind = dp ), intent(inout), allocatable :: f(:,:)

      ! local variables
      integer i
      real( kind = dp) nout

      do i = 1, ubound(f,2)
         nout = c1d_norm( f, i )
         f(:,i) = f(:,i)*1.0_dp/SQRT( nout )
      enddo

   end subroutine

   subroutine c2d_normalize( f )

      implicit none
      complex( kind = dp ), intent(inout), allocatable :: f(:,:,:)

      ! local variables
      integer i, j
      real( kind = dp) nout

      do i = 1, ubound(f,2)
         do j = 1, ubound(f,3)
            nout = c2d_norm( f, i, j )
            f(:,i,j) = f(:,i,j)*1.0_dp/SQRT( nout )
         enddo
      enddo

   end subroutine

   subroutine c3d_normalize( f )

      implicit none
      complex( kind = dp ), intent(inout), allocatable :: f(:,:,:,:)

      ! local variables
      integer i, j, k
      real( kind = dp) nout

      do i = 1, ubound(f,2)
         do j = 1, ubound(f,3)
            do k = 1, ubound(f,4)
               nout = c3d_norm( f, i, j, k )
               f(:,i,j,k) = f(:,i,j,k)*1.0_dp/SQRT( nout )
            enddo
         enddo
      enddo

   end subroutine

!  | ----------------------------------- |
!  | Local order parameter of comp. A    |
!  | ----------------------------------- |
   function c1d_Order( f, i )

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:)
      integer, intent(in) :: i
      complex( kind = dp ) :: c1d_Order

      ! local variables
      integer  na

      c1d_Order = 0.0_dp*re

      do na = 1, ubound(f,1)
         c1d_Order = c1d_Order + SQRT( real( na,dp ) ) * conjg(f(na-1,i)) * f(na,i)
      enddo

   end function

   function c2d_Order( f, i, j )

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:,:)
      integer, intent(in) :: i, j
      complex( kind = dp ) :: c2d_Order

      ! local variables
      integer  na

      c2d_Order = 0.0_dp*re

      do na = 1, ubound(f,1)
         c2d_Order = c2d_Order + SQRT( real( na,dp ) ) * conjg(f(na-1,i,j)) * f(na,i,j)
      enddo

   end function

   function c3d_Order( f, i, j, k )

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:,:,:)
      integer, intent(in) :: i, j, k
      complex( kind = dp ) :: c3d_Order

      ! local variables
      integer  na

      c3d_Order = 0.0_dp*re

      do na = 1, ubound(f,1)
         c3d_Order = c3d_Order + SQRT( real( na,dp ) ) * conjg(f(na-1,i,j,k)) * f(na,i,j,k)
      enddo

   end function

!  | ----------------------------------- |
!  | Total Energy                        |
!  | ----------------------------------- |
   function c1d_TotEnergy( f, ja, u )

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:)
      real( kind = dp ), intent(in) :: ja, u
      real( kind = dp ) :: c1d_TotEnergy

      ! local variables
      integer i, na
      complex( kind = dp ) ord
      complex( kind = dp ), allocatable :: ArrayOrder(:)

      allocate( ArrayOrder( size(f,2) ) )

      c1d_TotEnergy = 0.0_dp

      do i = 1, ubound(f,2)

         ! calculate order parameters
         ord = Order( f, i )
         ArrayOrder(i) = ord

         do na = 0, ubound(f,1)
            c1d_TotEnergy = c1d_TotEnergy + u/2.0_dp * real( na*(na-1),dp ) * abs(f(na,i))**2
         enddo

      enddo

      do i = 1, ubound(f,2)

         if( i < ubound(f,2) ) then

            c1d_TotEnergy = c1d_TotEnergy - 2.0_dp * ja * real( conjg(ArrayOrder(i)) * ArrayOrder(i+1) )

         else

            c1d_TotEnergy = c1d_TotEnergy - 2.0_dp * ja * real( conjg(ArrayOrder( ubound(f,2) )) * ArrayOrder(1) )

         endif

      enddo

   end function

   function c2d_TotEnergy( f, ja, u)

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:,:)
      real( kind = dp ), intent(in) :: ja, u
      real( kind = dp ) :: c2d_TotEnergy

      ! local variables
      integer i, j, na
      complex( kind = dp ) ord
      complex( kind = dp ), allocatable :: ArrayOrder(:,:)

      allocate( ArrayOrder( size(f,2), size(f,3) ) )

      c2d_TotEnergy = 0.0_dp

      do i = 1, ubound(f,2)
         do j = 1, ubound(f,3)

            ! calculate order parameters
            ord = Order( f, i, j )
            ArrayOrder(i,j) = ord

            do na = 0, ubound(f,1)
               c2d_TotEnergy = c2d_TotEnergy + u/2.0_dp * real( na*(na-1),dp ) * abs(f(na,i,j))**2
            enddo

         enddo
      enddo

      do i = 1, ubound(f,2)
         do j = 1, ubound(f,3)

            if( i < ubound(f,2) ) then

               c2d_TotEnergy = c2d_TotEnergy - 2.0_dp * ja * real( conjg(ArrayOrder(i,j)) * ArrayOrder(i+1,j) )

            else

               c2d_TotEnergy = c2d_TotEnergy - 2.0_dp * ja * real( conjg(ArrayOrder( ubound(f,2), j )) * ArrayOrder(1,j) )

            endif

            if( j < ubound(f,3) ) then

               c2d_TotEnergy = c2d_TotEnergy - 2.0_dp * ja * real( conjg(ArrayOrder(i,j)) * ArrayOrder(i,j+1) )

           else

               c2d_TotEnergy = c2d_TotEnergy - 2.0_dp * ja * real( conjg(ArrayOrder( i, ubound(f,3) )) * ArrayOrder(i,1) )

           endif

         enddo
      enddo

   end function

   function c3d_TotEnergy( f, ja, u )

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:,:,:)
      real( kind = dp ), intent(in) :: ja, u
      real( kind = dp ) :: c3d_TotEnergy

      ! local variables
      integer i, j, k, na
      complex( kind = dp ) ord
      complex( kind = dp ), allocatable :: ArrayOrder(:,:,:)

      allocate( ArrayOrder( size(f,2), size(f,3), size(f,4) ) )

      c3d_TotEnergy = 0.0_dp

      do i = 1, ubound(f,2)
         do j = 1, ubound(f,3)
            do k = 1, ubound(f,4)

               ! calculate order parameters
               ord = Order( f, i, j, k )
               ArrayOrder(i,j,k) = ord

               do na = 0, ubound(f,1)
                     c3d_TotEnergy = c3d_TotEnergy + u/2.0_dp * real( na*(na-1),dp ) * abs(f(na,i,j,k))**2
               enddo

            enddo
         enddo
      enddo

      do i = 1, ubound(f,2)
         do j = 1, ubound(f,3)
            do k = 1, ubound(f,4)

               if( i < ubound(f,2) ) then

                  c3d_TotEnergy = c3d_TotEnergy - 2.0_dp * ja * real( conjg(ArrayOrder(i,j,k)) * ArrayOrder(i+1,j,k) )

               else

                  c3d_TotEnergy = c3d_TotEnergy - 2.0_dp * ja * real( conjg(ArrayOrder( ubound(f,2),j,k) ) * ArrayOrder(1,j,k) )

               endif

               if( j < ubound(f,3) ) then

                  c3d_TotEnergy = c3d_TotEnergy - 2.0_dp * ja * real( conjg(ArrayOrder(i,j,k)) * ArrayOrder(i,j+1,k) )

               else

                  c3d_TotEnergy = c3d_TotEnergy - 2.0_dp * ja * real( conjg(ArrayOrder( i,ubound(f,3),k) ) * ArrayOrder(i,1,k) )

               endif

               if( k < ubound(f,4) ) then

                  c3d_TotEnergy = c3d_TotEnergy - 2.0_dp * ja * real( conjg(ArrayOrder(i,j,k)) * ArrayOrder(i,j,k+1) )

               else

                  c3d_TotEnergy = c3d_TotEnergy - 2.0_dp * ja * real( conjg(ArrayOrder( i,j,ubound(f,4)) ) * ArrayOrder(i,j,1) )

               endif

            enddo
         enddo
      enddo

   end function

! | -------------------------------------------------------- |
! | Definition of the chemical potential based on Gutzwiller |
! | -------------------------------------------------------- |
   function c1d_ChemPot( f, ja, u )

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:)
      real( kind = dp ), intent(in) :: ja, u
      real( kind = dp ) :: c1d_ChemPot

      integer :: i, n
      complex( kind = dp ) :: ord, phi, xi
      complex( kind = dp ), allocatable :: ArrayOrder(:)
      real( kind = dp ) :: sum1, sum2, sum3, tot1, tot2, tot3

      allocate( ArrayOrder( size(f,2) ) )

      ! collect order parameters
      do i = 1, ubound(f,2)
         ArrayOrder(i) = c1d_Order( f, i )
      enddo

      tot1 = 0.0_dp
      tot2 = 0.0_dp
      tot3 = 0.0_dp

      do i = lbound(f,2), ubound(f,2)
         ! nearest neighbour order parameter
         if( i == 1 ) then
            ord = ArrayOrder( i+1 ) + ArrayOrder( ubound(f,2) )
         elseif( i == ubound(f,2) ) then
            ord = ArrayOrder( 1 ) + ArrayOrder( i-1 )
         else 
            ord = ArrayOrder( i-1 ) + ArrayOrder( i+1 )
         endif

         sum1 = 0.0_dp; sum2 = 0.0_dp; sum3 = 0.0_dp; phi = 0.0_dp; xi = 0.0_dp
         do n = 0, ubound(f,1)
            sum1 = sum1 + real(n,dp) * abs(f(n,i))**2
            sum2 = sum2 + real(n,dp)**2 * abs(f(n,i))**2
            sum3 = sum3 + real(n,dp)**3 * abs(f(n,i))**2
            if (n > 0) then
               phi = phi + sqrt(real(n,dp)) * f(n,i) * conjg(f(n-1,i))
               xi  = xi + real(n,dp)*sqrt(real(n,dp)) * f(n,i) * conjg(f(n-1,i))
            endif
         enddo
         tot1 = tot1 + (sum2 - sum1**2)
         tot2 = tot2 + (sum3 + sum1**2 - sum2 - sum2*sum1)
         tot3 = tot3 + real(conjg(ord)*(phi + 2.0_dp*phi*sum1 - 2.0_dp*xi))
      enddo

      c1d_ChemPot = 1.0_dp/tot1 * (u/2.0_dp*tot2 + ja*tot3)

   end function

   function c2d_ChemPot( f, ja, u )

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:,:)
      real( kind = dp ), intent(in) :: ja, u
      real( kind = dp ) :: c2d_ChemPot

      integer :: i, j, n
      complex( kind = dp ) :: ord, phi, xi
      complex( kind = dp ), allocatable :: ArrayOrder(:,:)
      real( kind = dp ) :: sum1, sum2, sum3, tot1, tot2, tot3

      allocate( ArrayOrder( size(f,2), size(f,3) ) )

      ! collect order parameters
      do i = 1, ubound(f,2)
         do j = 1, ubound(f,3)
            ArrayOrder(i,j) = c2d_Order( f, i, j )
         enddo
      enddo

      tot1 = 0.0_dp
      tot2 = 0.0_dp
      tot3 = 0.0_dp

      do i = 1, ubound(f,2)
      do j = 1, ubound(f,3)

         ! nearest neighbour order parameter
         if( i == 1 ) then
            ord = ArrayOrder( i+1,j ) + ArrayOrder( ubound(f,2),j )
         elseif( i == ubound(f,2) ) then
            ord = ArrayOrder( 1,j ) + ArrayOrder( i-1,j )
         else 
            ord = ArrayOrder( i-1,j ) + ArrayOrder( i+1,j )
         endif

         if( j == 1 ) then
            ord = ord + ArrayOrder( i, j+1 ) + ArrayOrder( i, ubound(f,3) )
         elseif( j == ubound(f,3) ) then
            ord = ord + ArrayOrder( i, 1 ) + ArrayOrder( i, j-1 )
         else 
            ord = ord + ArrayOrder( i, j-1 ) + ArrayOrder( i, j+1 )
         endif

         sum1 = 0.0_dp; sum2 = 0.0_dp; sum3 = 0.0_dp; phi = 0.0_dp; xi = 0.0_dp
         do n = 0, ubound(f,1)
            sum1 = sum1 + real(n,dp) * abs(f(n,i,j))**2
            sum2 = sum2 + real(n,dp)**2 * abs(f(n,i,j))**2
            sum3 = sum3 + real(n,dp)**3 * abs(f(n,i,j))**2
            if (n > 0) then
               phi = phi + sqrt(real(n,dp)) * f(n,i,j) * conjg(f(n-1,i,j))
               xi  = xi + real(n,dp)*sqrt(real(n,dp)) * f(n,i,j) * conjg(f(n-1,i,j))
            endif
         enddo
         tot1 = tot1 + (sum2 - sum1**2)
         tot2 = tot2 + (sum3 + sum1**2 - sum2 - sum2*sum1)
         tot3 = tot3 + real(conjg(ord)*(phi + 2.0_dp*phi*sum1 - 2.0_dp*xi))
      enddo ! i
      enddo ! j

      c2d_ChemPot = 1.0_dp/tot1 * (u/2.0_dp*tot2 + ja*tot3)

   end function

   function c3d_ChemPot( f, ja, u )

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:,:,:)
      real( kind = dp ), intent(in) :: ja, u
      real( kind = dp ) :: c3d_ChemPot

      integer :: i, j, k, n
      complex( kind = dp ) :: ord, phi, xi
      complex( kind = dp ), allocatable :: ArrayOrder(:,:,:)
      real( kind = dp ) :: sum1, sum2, sum3, tot1, tot2, tot3

      allocate( ArrayOrder( size(f,2), size(f,3), size(f,4) ) )

      ! collect order parameters
      do i = 1, ubound(f,2)
         do j = 1, ubound(f,3)
            do k = 1, ubound(f,4)
               ArrayOrder(i,j,k) = c3d_Order( f, i, j, k )
            enddo
         enddo
      enddo

      tot1 = 0.0_dp
      tot2 = 0.0_dp
      tot3 = 0.0_dp

      do i = 1, ubound(f,2)
      do j = 1, ubound(f,3)
      do k = 1, ubound(f,4)

         ! nearest neighbour order parameter
         if( i == 1 ) then
            ord = ArrayOrder( i+1,j,k ) + ArrayOrder( ubound(f,2),j,k )
         elseif( i == ubound(f,2) ) then
            ord = ArrayOrder( 1,j,k ) + ArrayOrder( i-1,j,k )
         else 
            ord = ArrayOrder( i-1,j,k ) + ArrayOrder( i+1,j,k )
         endif

         if( j == 1 ) then
            ord = ord + ArrayOrder( i, j+1,k ) + ArrayOrder( i, ubound(f,3),k )
         elseif( j == ubound(f,3) ) then
            ord = ord + ArrayOrder( i, 1,k ) + ArrayOrder( i, j-1,k )
         else 
            ord = ord + ArrayOrder( i, j-1,k ) + ArrayOrder( i, j+1,k )
         endif

         if( k == 1 ) then
            ord = ord + ArrayOrder( i,j,k+1 ) + ArrayOrder( i,j,ubound(f,4) )
         elseif( k == ubound(f,4) ) then
            ord = ord + ArrayOrder( i,j,1 ) + ArrayOrder( i,j,k-1 )
         else 
            ord = ord + ArrayOrder( i,j,k-1 ) + ArrayOrder( i,j,k+1 )
         endif

         sum1 = 0.0_dp; sum2 = 0.0_dp; sum3 = 0.0_dp; phi = 0.0_dp; xi = 0.0_dp
         do n = 0, ubound(f,1)
            sum1 = sum1 + real(n,dp) * abs(f(n,i,j,k))**2
            sum2 = sum2 + real(n,dp)**2 * abs(f(n,i,j,k))**2
            sum3 = sum3 + real(n,dp)**3 * abs(f(n,i,j,k))**2
            if (n > 0) then
               phi = phi + sqrt(real(n,dp)) * f(n,i,j,k) * conjg(f(n-1,i,j,k))
               xi  = xi + real(n,dp)*sqrt(real(n,dp)) * f(n,i,j,k) * conjg(f(n-1,i,j,k))
            endif
         enddo
         tot1 = tot1 + (sum2 - sum1**2)
         tot2 = tot2 + (sum3 + sum1**2 - sum2 - sum2*sum1)
         tot3 = tot3 + real(conjg(ord)*(phi + 2.0_dp*phi*sum1 - 2.0_dp*xi))
      enddo ! i
      enddo ! j
      enddo ! k

      c3d_ChemPot = 1.0_dp/tot1 * (u/2.0_dp*tot2 + ja*tot3)

   end function

end module
