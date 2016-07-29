module measures

   use parameters, only: dp, re, im

   private
   public MeanA, MeanB, MeanAB, VarA, VarB, norm, normalize, OrderA, OrderB, TotEnergy

   ! | ---------------------------------------- |
   ! | For better memeory layout it is better   |
   ! | to put Fock basis index at the beginning |
   ! | ---------------------------------------- |
   interface MeanA
      module procedure c1d_MeanA
      module procedure c2d_MeanA
      module procedure c3d_MeanA
   end interface

   interface MeanB
      module procedure c1d_MeanB
      module procedure c2d_MeanB
      module procedure c3d_MeanB
   end interface

   interface MeanAB
      module procedure c1d_MeanAB
      module procedure c2d_MeanAB
      module procedure c3d_MeanAB
   end interface

   interface VarA
      module procedure c1d_VarA
      module procedure c2d_VarA
      module procedure c3d_VarA
   end interface 

   interface VarB
      module procedure c1d_VarB
      module procedure c2d_VarB
      module procedure c3d_VarB
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

   interface OrderA
      module procedure c1d_OrderA
      module procedure c2d_OrderA
      module procedure c3d_OrderA
   end interface

   interface OrderB
      module procedure c1d_OrderB
      module procedure c2d_OrderB
      module procedure c3d_OrderB
   end interface

   interface TotEnergy
      module procedure c1d_TotEnergy
      module procedure c2d_TotEnergy
      module procedure c3d_TotEnergy
   end interface  

contains

!  | ----------------------------------- |
!  | Mean number of particles in comp. A |
!  | ----------------------------------- |
   function c1d_MeanA( f, i )

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:,:)
      integer, intent(in) :: i
      real( kind = dp ) :: c1d_MeanA

      ! local variables
      integer  na, nb

      c1d_MeanA = 0.0_dp
      do na = 0, ubound(f,1)
         do nb = 0, ubound(f,2)
            c1d_MeanA = c1d_MeanA + real( na,dp ) * abs(f(na,nb,i))**2
         enddo
      enddo

   end function

   function c2d_MeanA( f, i, j )

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:,:,:)
      integer, intent(in) :: i, j
      real( kind = dp ) :: c2d_MeanA

      ! local variables
      integer  na, nb

      c2d_MeanA = 0.0_dp
      do na = 0, ubound(f,1)
         do nb = 0, ubound(f,2)
            c2d_MeanA = c2d_MeanA + real( na,dp ) * abs(f(na,nb,i,j))**2
         enddo
      enddo

   end function

   function c3d_MeanA( f, i, j, k )

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:,:,:,:)
      integer, intent(in) :: i, j, k
      real( kind = dp ) :: c3d_MeanA

      ! local variables
      integer  na, nb

      c3d_MeanA = 0.0_dp
      do na = 0, ubound(f,1)
         do nb = 0, ubound(f,2)
            c3d_MeanA = c3d_MeanA + real( na,dp ) * abs(f(na,nb,i,j,k))**2
         enddo
      enddo

   end function

!  | ----------------------------------- |
!  | Mean number of particles in comp. B |
!  | ----------------------------------- |
   function c1d_MeanB( f, i )

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:,:)
      integer, intent(in) :: i
      real( kind = dp ) :: c1d_MeanB

      ! local variables
      integer  na, nb

      c1d_MeanB = 0.0_dp
      do na = 0, ubound(f,1)
         do nb = 0, ubound(f,2)
            c1d_MeanB = c1d_MeanB + real( nb,dp ) * abs(f(na,nb,i))**2
         enddo
      enddo

   end function

   function c2d_MeanB( f, i, j )

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:,:,:)
      integer, intent(in) :: i, j
      real( kind = dp ) :: c2d_MeanB

      ! local variables
      integer  na, nb

      c2d_MeanB = 0.0_dp
      do na = 0, ubound(f,1)
         do nb = 0, ubound(f,2)
            c2d_MeanB = c2d_MeanB + real( nb,dp ) * abs(f(na,nb,i,j))**2
         enddo
      enddo

   end function

   function c3d_MeanB( f, i, j, k )

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:,:,:,:)
      integer, intent(in) :: i, j, k
      real( kind = dp ) :: c3d_MeanB

      ! local variables
      integer  na, nb

      c3d_MeanB = 0.0_dp
      do na = 0, ubound(f,1)
         do nb = 0, ubound(f,2)
            c3d_MeanB = c3d_MeanB + real( nb,dp ) * abs(f(na,nb,i,j,k))**2
         enddo
      enddo

   end function

!  | ----------------------------------- |
!  | Mean number <na*nb>                 |
!  | ----------------------------------- |
   function c1d_MeanAB( f, i )

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:,:)
      integer, intent(in) :: i
      real( kind = dp ) :: c1d_MeanAB

      ! local variables
      integer  na, nb

      c1d_MeanAB = 0.0_dp
      do na = 0, ubound(f,1)
         do nb = 0, ubound(f,2)
            c1d_MeanAB = c1d_MeanAB + real( na*nb,dp ) * abs(f(na,nb,i))**2
         enddo
      enddo

   end function

   function c2d_MeanAB( f, i, j )

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:,:,:)
      integer, intent(in) :: i, j
      real( kind = dp ) :: c2d_MeanAB

      ! local variables
      integer  na, nb

      c2d_MeanAB = 0.0_dp
      do na = 0, ubound(f,1)
         do nb = 0, ubound(f,2)
            c2d_MeanAB = c2d_MeanAB + real( na*nb,dp ) * abs(f(na,nb,i,j))**2
         enddo
      enddo

   end function

   function c3d_MeanAB( f, i, j, k )

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:,:,:,:)
      integer, intent(in) :: i, j, k
      real( kind = dp ) :: c3d_MeanAB

      ! local variables
      integer  na, nb

      c3d_MeanAB = 0.0_dp
      do na = 0, ubound(f,1)
         do nb = 0, ubound(f,2)
            c3d_MeanAB = c3d_MeanAB + real( na*nb,dp ) * abs(f(na,nb,i,j,k))**2
         enddo
      enddo

   end function

!  | ----------------------------------- |
!  | Variance of particle number in A    |
!  | ----------------------------------- |
   function c1d_VarA( f, i )

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:,:)
      integer, intent(in) :: i
      real( kind = dp ) :: c1d_VarA

      ! local variables
      integer  na, nb
      real( kind = dp ) suma

      c1d_VarA = 0.0_dp
      suma   = 0.0_dp

      do na = 0, ubound(f,1)
         do nb = 0, ubound(f,2)
            suma   = suma + real( na,dp ) * abs(f(na,nb,i))**2
            c1d_VarA = c1d_VarA + real( na*na,dp ) * abs(f(na,nb,i))**2
         enddo
      enddo

      c1d_VarA = c1d_VarA - suma**2

   end function

   function c2d_VarA( f, i, j )

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:,:,:)
      integer, intent(in) :: i, j
      real( kind = dp ) :: c2d_VarA

      ! local variables
      integer  na, nb
      real( kind = dp ) suma

      c2d_VarA = 0.0_dp
      suma   = 0.0_dp

      do na = 0, ubound(f,1)
         do nb = 0, ubound(f,2)
            suma   = suma + real( na,dp ) * abs(f(na,nb,i,j))**2
            c2d_VarA = c2d_VarA + real( na*na,dp ) * abs(f(na,nb,i,j))**2
         enddo
      enddo

      c2d_VarA = c2d_VarA - suma**2

   end function

   function c3d_VarA( f, i, j, k )

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:,:,:,:)
      integer, intent(in) :: i, j, k
      real( kind = dp ) :: c3d_VarA

      ! local variables
      integer  na, nb
      real( kind = dp ) suma

      c3d_VarA = 0.0_dp
      suma   = 0.0_dp

      do na = 0, ubound(f,1)
         do nb = 0, ubound(f,2)
            suma   = suma + real( na,dp ) * abs(f(na,nb,i,j,k))**2
            c3d_VarA = c3d_VarA + real( na*na,dp ) * abs(f(na,nb,i,j,k))**2
         enddo
      enddo

      c3d_VarA = c3d_VarA - suma**2

   end function

!  | ----------------------------------- |
!  | Variance of particle number in B    |
!  | ----------------------------------- |
   function c1d_VarB( f, i )

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:,:)
      integer, intent(in) :: i
      real( kind = dp ) :: c1d_VarB

      ! local variables
      integer  na, nb
      real( kind = dp ) suma

      c1d_VarB = 0.0_dp
      suma   = 0.0_dp

      do na = 0, ubound(f,1)
         do nb = 0, ubound(f,2)
            suma   = suma + real( nb,dp ) * abs(f(na,nb,i))**2
            c1d_VarB = c1d_VarB + real( nb*nb,dp ) * abs(f(na,nb,i))**2
         enddo
      enddo

      c1d_VarB = c1d_VarB - suma**2

   end function

   function c2d_VarB( f, i, j )

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:,:,:)
      integer, intent(in) :: i, j
      real( kind = dp ) :: c2d_VarB

      ! local variables
      integer  na, nb
      real( kind = dp ) suma

      c2d_VarB = 0.0_dp
      suma   = 0.0_dp

      do na = 0, ubound(f,1)
         do nb = 0, ubound(f,2)
            suma   = suma + real( nb,dp ) * abs(f(na,nb,i,j))**2
            c2d_VarB = c2d_VarB + real( nb*nb,dp ) * abs(f(na,nb,i,j))**2
         enddo
      enddo

      c2d_VarB = c2d_VarB - suma**2

   end function

   function c3d_VarB( f, i, j, k )

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:,:,:,:)
      integer, intent(in) :: i, j, k
      real( kind = dp ) :: c3d_VarB

      ! local variables
      integer  na, nb
      real( kind = dp ) suma

      c3d_VarB = 0.0_dp
      suma   = 0.0_dp

      do na = 0, ubound(f,1)
         do nb = 0, ubound(f,2)
            suma   = suma + real( nb,dp ) * abs(f(na,nb,i,j,k))**2
            c3d_VarB = c3d_VarB + real( nb*nb,dp ) * abs(f(na,nb,i,j,k))**2
         enddo
      enddo

      c3d_VarB = c3d_VarB - suma**2

   end function

!  | ----------------------------------- |
!  | Local norm of Gutzwiller wavefun    |
!  | ----------------------------------- |
   function c1d_norm( f, i )

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:,:)
      integer, intent(in) :: i
      real( kind = dp ) :: c1d_norm

      ! local variables
      integer  na, nb

      c1d_norm = 0.0_dp

      do na = 0, ubound(f,1)
         do nb = 0, ubound(f,2)
            c1d_norm = c1d_norm + abs(f(na,nb,i))**2
         enddo
      enddo

   end function

   function c2d_norm( f, i, j )

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:,:,:)
      integer, intent(in) :: i, j
      real( kind = dp ) :: c2d_norm

      ! local variables
      integer  na, nb

      c2d_norm = 0.0_dp

      do na = 0, ubound(f,1)
         do nb = 0, ubound(f,2)
            c2d_norm = c2d_norm + abs(f(na,nb,i,j))**2
         enddo
      enddo

   end function

   function c3d_norm( f, i, j, k )

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:,:,:,:)
      integer, intent(in) :: i, j, k
      real( kind = dp ) :: c3d_norm

      ! local variables
      integer  na, nb

      c3d_norm = 0.0_dp

      do na = 0, ubound(f,1)
         do nb = 0, ubound(f,2)
            c3d_norm = c3d_norm + abs(f(na,nb,i,j,k))**2
         enddo
      enddo

   end function

!  | ----------------------------------- |
!  | normalize Gutzwiller wavefun        |
!  | ----------------------------------- |
   subroutine c1d_normalize( f )

      implicit none
      complex( kind = dp ), intent(inout), allocatable :: f(:,:,:)

      ! local variables
      integer i
      real( kind = dp) nout

      do i = 1, ubound(f,3)
         nout = c1d_norm( f, i )
         f(:,:,i) = f(:,:,i)*1.0_dp/SQRT( nout )
      enddo

   end subroutine

   subroutine c2d_normalize( f )

      implicit none
      complex( kind = dp ), intent(inout), allocatable :: f(:,:,:,:)

      ! local variables
      integer i, j
      real( kind = dp) nout

      do i = 1, ubound(f,3)
         do j = 1, ubound(f,4)
            nout = c2d_norm( f, i, j )
            f(:,:,i,j) = f(:,:,i,j)*1.0_dp/SQRT( nout )
         enddo
      enddo

   end subroutine

   subroutine c3d_normalize( f )

      implicit none
      complex( kind = dp ), intent(inout), allocatable :: f(:,:,:,:,:)

      ! local variables
      integer i, j, k
      real( kind = dp) nout

      do i = 1, ubound(f,3)
         do j = 1, ubound(f,4)
            do k = 1, ubound(f,5)
               nout = c3d_norm( f, i, j, k )
               f(:,:,i,j,k) = f(:,:,i,j,k)*1.0_dp/SQRT( nout )
            enddo
         enddo
      enddo

   end subroutine
!  | ----------------------------------- |
!  | Local order parameter of comp. A    |
!  | ----------------------------------- |
   function c1d_orderA( f, i )

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:,:)
      integer, intent(in) :: i
      complex( kind = dp ) :: c1d_orderA

      ! local variables
      integer  na, nb

      c1d_orderA = 0.0_dp*re

      do na = 1, ubound(f,1)
         do nb = 0, ubound(f,2)
            c1d_orderA = c1d_orderA + SQRT( real( na,dp ) ) * conjg(f(na-1,nb,i)) * f(na,nb,i)
         enddo
      enddo

   end function

   function c2d_orderA( f, i, j )

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:,:,:)
      integer, intent(in) :: i, j
      complex( kind = dp ) :: c2d_orderA

      ! local variables
      integer  na, nb

      c2d_orderA = 0.0_dp*re

      do na = 1, ubound(f,1)
         do nb = 0, ubound(f,2)
            c2d_orderA = c2d_orderA + SQRT( real( na,dp ) ) * conjg(f(na-1,nb,i,j)) * f(na,nb,i,j)
         enddo
      enddo

   end function

   function c3d_orderA( f, i, j, k )

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:,:,:,:)
      integer, intent(in) :: i, j, k
      complex( kind = dp ) :: c3d_orderA

      ! local variables
      integer  na, nb

      c3d_orderA = 0.0_dp*re

      do na = 1, ubound(f,1)
         do nb = 0, ubound(f,2)
            c3d_orderA = c3d_orderA + SQRT( real( na,dp ) ) * conjg(f(na-1,nb,i,j,k)) * f(na,nb,i,j,k)
         enddo
      enddo

   end function

!  | ----------------------------------- |
!  | Local order parameter of comp. B    |
!  | ----------------------------------- |
   function c1d_orderB( f, i )

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:,:)
      integer, intent(in) :: i
      complex( kind = dp ) :: c1d_orderB

      ! local variables
      integer  na, nb

      c1d_orderB = 0.0_dp*re

      do na = 0, ubound(f,1)
         do nb = 1, ubound(f,2)
            c1d_orderB = c1d_orderB + SQRT( real( nb,dp ) ) * conjg(f(na,nb-1,i)) * f(na,nb,i)
         enddo
      enddo

   end function

   function c2d_orderB( f, i, j )

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:,:,:)
      integer, intent(in) :: i, j
      complex( kind = dp ) :: c2d_orderB

      ! local variables
      integer  na, nb

      c2d_orderB = 0.0_dp*re

      do na = 0, ubound(f,1)
         do nb = 1, ubound(f,2)
            c2d_orderB = c2d_orderB + SQRT( real( nb,dp ) ) * conjg(f(na,nb-1,i,j)) * f(na,nb,i,j)
         enddo
      enddo

   end function

   function c3d_orderB( f, i, j, k )

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:,:,:,:)
      integer, intent(in) :: i, j, k
      complex( kind = dp ) :: c3d_orderB

      ! local variables
      integer  na, nb

      c3d_orderB = 0.0_dp*re

      do na = 0, ubound(f,1)
         do nb = 1, ubound(f,2)
            c3d_orderB = c3d_orderB + SQRT( real( nb,dp ) ) * conjg(f(na,nb-1,i,j,k)) * f(na,nb,i,j,k)
         enddo
      enddo

   end function

!  | ----------------------------------- |
!  | Total Energy                        |
!  | ----------------------------------- |
   function c1d_TotEnergy( f, ja, jb, ua, ub, uab )

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:,:)
      real( kind = dp ), intent(in) :: ja, jb, ua, ub, uab
      real( kind = dp ) :: c1d_TotEnergy

      ! local variables
      integer i, na, nb
      complex( kind = dp ) ord
      complex( kind = dp ), allocatable :: ArrayOrderA(:), ArrayOrderB(:)

      allocate( ArrayOrderA( size(f,3) ), ArrayOrderB( size(f,3) ) )

      c1d_TotEnergy = 0.0_dp

      do i = 1, ubound(f,3)

         ! calculate order parameters
         ord = OrderA( f, i )
         ArrayOrderA(i) = ord
         ord = OrderB( f, i )
         ArrayOrderB(i) = ord

         do na = 0, ubound(f,1)
            do nb = 0, ubound(f,2)

               c1d_TotEnergy = c1d_TotEnergy + ua/2.0_dp * real( na*(na-1),dp ) * abs(f(na,nb,i))**2
               c1d_TotEnergy = c1d_TotEnergy + ub/2.0_dp * real( nb*(nb-1),dp ) * abs(f(na,nb,i))**2
               c1d_TotEnergy = c1d_TotEnergy + uab * real( na*nb,dp ) * abs(f(na,nb,i))**2

            enddo
         enddo

      enddo

      do i = 1, ubound(f,3)

         if( i < ubound(f,3) ) then

            c1d_TotEnergy = c1d_TotEnergy - 2.0_dp * ja * real( conjg(ArrayOrderA(i)) * ArrayOrderA(i+1) )
            c1d_TotEnergy = c1d_TotEnergy - 2.0_dp * jb * real( conjg(ArrayOrderB(i)) * ArrayOrderB(i+1) )

         else

            c1d_TotEnergy = c1d_TotEnergy - 2.0_dp * ja * real( conjg(ArrayOrderA( ubound(f,3) )) * ArrayOrderA(1) )
            c1d_TotEnergy = c1d_TotEnergy - 2.0_dp * jb * real( conjg(ArrayOrderB( ubound(f,3) )) * ArrayOrderB(1) )

         endif

      enddo

   end function

   function c2d_TotEnergy( f, ja, jb, ua, ub, uab )

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:,:,:)
      real( kind = dp ), intent(in) :: ja, jb, ua, ub, uab
      real( kind = dp ) :: c2d_TotEnergy

      ! local variables
      integer i, j, na, nb
      complex( kind = dp ) ord
      complex( kind = dp ), allocatable :: ArrayOrderA(:,:), ArrayOrderB(:,:)

      allocate( ArrayOrderA( size(f,3), size(f,4) ), ArrayOrderB( size(f,3), size(f,4) ) )

      c2d_TotEnergy = 0.0_dp

      do i = 1, ubound(f,3)
         do j = 1, ubound(f,4)

            ! calculate order parameters
            ord = OrderA( f, i, j )
            ArrayOrderA(i,j) = ord
            ord = OrderB( f, i, j )
            ArrayOrderB(i,j) = ord

            do na = 0, ubound(f,1)
               do nb = 0, ubound(f,2)

                  c2d_TotEnergy = c2d_TotEnergy + ua/2.0_dp * real( na*(na-1),dp ) * abs(f(na,nb,i,j))**2
                  c2d_TotEnergy = c2d_TotEnergy + ub/2.0_dp * real( nb*(nb-1),dp ) * abs(f(na,nb,i,j))**2
                  c2d_TotEnergy = c2d_TotEnergy + uab * real( na*nb,dp ) * abs(f(na,nb,i,j))**2

               enddo
            enddo

         enddo
      enddo

      do i = 1, ubound(f,3)
         do j = 1, ubound(f,4)

            if( i < ubound(f,3 ) ) then

               c2d_TotEnergy = c2d_TotEnergy - 2.0_dp * ja * real( conjg(ArrayOrderA(i,j)) * ArrayOrderA(i+1,j) )
               c2d_TotEnergy = c2d_TotEnergy - 2.0_dp * jb * real( conjg(ArrayOrderB(i,j)) * ArrayOrderB(i+1,j) )

            else

               c2d_TotEnergy = c2d_TotEnergy - 2.0_dp * ja * real( conjg(ArrayOrderA( ubound(f,3), j )) * ArrayOrderA(1,j) )
               c2d_TotEnergy = c2d_TotEnergy - 2.0_dp * jb * real( conjg(ArrayOrderB( ubound(f,3), j )) * ArrayOrderB(1,j) )

            endif

            if( j < ubound(f,4) ) then

               c2d_TotEnergy = c2d_TotEnergy - 2.0_dp * ja * real( conjg(ArrayOrderA(i,j)) * ArrayOrderA(i,j+1) )
               c2d_TotEnergy = c2d_TotEnergy - 2.0_dp * jb * real( conjg(ArrayOrderB(i,j)) * ArrayOrderB(i,j+1) )

           else

               c2d_TotEnergy = c2d_TotEnergy - 2.0_dp * ja * real( conjg(ArrayOrderA( i, ubound(f,4) )) * ArrayOrderA(i,1) )
               c2d_TotEnergy = c2d_TotEnergy - 2.0_dp * jb * real( conjg(ArrayOrderB( i, ubound(f,4) )) * ArrayOrderB(i,1) )

           endif

         enddo
      enddo

   end function

   function c3d_TotEnergy( f, ja, jb, ua, ub, uab )

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:,:,:,:)
      real( kind = dp ), intent(in) :: ja, jb, ua, ub, uab
      real( kind = dp ) :: c3d_TotEnergy

      ! local variables
      integer i, j, k, na, nb
      complex( kind = dp ) ord
      complex( kind = dp ), allocatable :: ArrayOrderA(:,:,:), ArrayOrderB(:,:,:)

      allocate( ArrayOrderA( size(f,3), size(f,4), size(f,5) ), ArrayOrderB( size(f,3), size(f,4), size(f,5) ) )

      c3d_TotEnergy = 0.0_dp

      do i = 1, ubound(f,3)
         do j = 1, ubound(f,4)
            do k = 1, ubound(f,5)

               ! calculate order parameters
               ord = OrderA( f, i, j, k )
               ArrayOrderA(i,j,k) = ord
               ord = OrderB( f, i, j, k )
               ArrayOrderB(i,j,k) = ord

               do na = 0, ubound(f,1)
                  do nb = 0, ubound(f,2)

                     c3d_TotEnergy = c3d_TotEnergy + ua/2.0_dp * real( na*(na-1),dp ) * abs(f(na,nb,i,j,k))**2
                     c3d_TotEnergy = c3d_TotEnergy + ub/2.0_dp * real( nb*(nb-1),dp ) * abs(f(na,nb,i,j,k))**2
                     c3d_TotEnergy = c3d_TotEnergy + uab * real( na*nb,dp ) * abs(f(na,nb,i,j,k))**2

                  enddo
               enddo

            enddo
         enddo
      enddo

      do i = 1, ubound(f,3)
         do j = 1, ubound(f,4)
            do k = 1, ubound(f,5)

               if( i < ubound(f,3) ) then

                  c3d_TotEnergy = c3d_TotEnergy - 2.0_dp * ja * real( conjg(ArrayOrderA(i,j,k)) * ArrayOrderA(i+1,j,k) )
                  c3d_TotEnergy = c3d_TotEnergy - 2.0_dp * jb * real( conjg(ArrayOrderB(i,j,k)) * ArrayOrderB(i+1,j,k) )

               else

                  c3d_TotEnergy = c3d_TotEnergy - 2.0_dp * ja * real( conjg(ArrayOrderA( ubound(f,3),j,k) ) * ArrayOrderA(1,j,k) )
                  c3d_TotEnergy = c3d_TotEnergy - 2.0_dp * jb * real( conjg(ArrayOrderB( ubound(f,3),j,k) ) * ArrayOrderB(1,j,k) )

               endif

               if( j < ubound(f,4) ) then

                  c3d_TotEnergy = c3d_TotEnergy - 2.0_dp * ja * real( conjg(ArrayOrderA(i,j,k)) * ArrayOrderA(i,j+1,k) )
                  c3d_TotEnergy = c3d_TotEnergy - 2.0_dp * jb * real( conjg(ArrayOrderB(i,j,k)) * ArrayOrderB(i,j+1,k) )

               else

                  c3d_TotEnergy = c3d_TotEnergy - 2.0_dp * ja * real( conjg(ArrayOrderA( i,ubound(f,4),k) ) * ArrayOrderA(i,1,k) )
                  c3d_TotEnergy = c3d_TotEnergy - 2.0_dp * jb * real( conjg(ArrayOrderB( i,ubound(f,4),k) ) * ArrayOrderB(i,1,k) )

               endif

               if( k < ubound(f,5) ) then

                  c3d_TotEnergy = c3d_TotEnergy - 2.0_dp * ja * real( conjg(ArrayOrderA(i,j,k)) * ArrayOrderA(i,j,k+1) )
                  c3d_TotEnergy = c3d_TotEnergy - 2.0_dp * jb * real( conjg(ArrayOrderB(i,j,k)) * ArrayOrderB(i,j,k+1) )

               else

                  c3d_TotEnergy = c3d_TotEnergy - 2.0_dp * ja * real( conjg(ArrayOrderA( i,j,ubound(f,5)) ) * ArrayOrderA(i,j,1) )
                  c3d_TotEnergy = c3d_TotEnergy - 2.0_dp * jb * real( conjg(ArrayOrderB( i,j,ubound(f,5)) ) * ArrayOrderB(i,j,1) )

               endif

            enddo
         enddo
      enddo

   end function

end module
