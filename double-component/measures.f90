module measures

   use parameters, only: dp, re, im

   private
   public MeanA, MeanB, MeanAB, VarA, VarB, norm, normalize, OrderA, OrderB, TotEnergy, ChemPot

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

   interface ChemPot
      module procedure c1d_ChemPot
      module procedure c2d_ChemPot
      module procedure c3d_ChemPot
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
      do nb = 0, ubound(f,2)
         do na = 0, ubound(f,1)
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
      do nb = 0, ubound(f,2)
         do na = 0, ubound(f,1)
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
      do nb = 0, ubound(f,2)
         do na = 0, ubound(f,1)
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
      do nb = 0, ubound(f,2)
         do na = 0, ubound(f,1)
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
      do nb = 0, ubound(f,2)
         do na = 0, ubound(f,1)
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
      do nb = 0, ubound(f,2)
         do na = 0, ubound(f,1)
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
      do nb = 0, ubound(f,2)
         do na = 0, ubound(f,1)
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
      do nb = 0, ubound(f,2)
         do na = 0, ubound(f,1)
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
      do nb = 0, ubound(f,2)
         do na = 0, ubound(f,1)
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

      do nb = 0, ubound(f,2)
         do na = 0, ubound(f,1)
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

      do nb = 0, ubound(f,2)
         do na = 0, ubound(f,1)
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

      do nb = 0, ubound(f,2)
         do na = 0, ubound(f,1)
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

      do nb = 0, ubound(f,2)
         do na = 0, ubound(f,1)
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

      do nb = 0, ubound(f,2)
         do na = 0, ubound(f,1)
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

      do nb = 0, ubound(f,2)
         do na = 0, ubound(f,1)
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

      do nb = 0, ubound(f,2)
         do na = 0, ubound(f,1)
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

      do nb = 0, ubound(f,2)
         do na = 0, ubound(f,1)
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

      do nb = 0, ubound(f,2)
         do na = 0, ubound(f,1)
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

      do j = 1, ubound(f,4)
         do i = 1, ubound(f,3)
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

      do k = 1, ubound(f,5)
         do j = 1, ubound(f,4)
            do i = 1, ubound(f,3)
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

      do nb = 0, ubound(f,2)
         do na = 1, ubound(f,1)
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

      do nb = 0, ubound(f,2)
         do na = 1, ubound(f,1)
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

      do nb = 0, ubound(f,2)
         do na = 1, ubound(f,1)
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

      do nb = 1, ubound(f,2)
         do na = 0, ubound(f,1)
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

      do nb = 1, ubound(f,2)
         do na = 0, ubound(f,1)
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

      do nb = 1, ubound(f,2)
         do na = 0, ubound(f,1)
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

         do nb = 0, ubound(f,2)
            do na = 0, ubound(f,1)

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

      do j = 1, ubound(f,4)
         do i = 1, ubound(f,3)

            ! calculate order parameters
            ord = OrderA( f, i, j )
            ArrayOrderA(i,j) = ord
            ord = OrderB( f, i, j )
            ArrayOrderB(i,j) = ord

            do nb = 0, ubound(f,2)
               do na = 0, ubound(f,1)

                  c2d_TotEnergy = c2d_TotEnergy + ua/2.0_dp * real( na*(na-1),dp ) * abs(f(na,nb,i,j))**2
                  c2d_TotEnergy = c2d_TotEnergy + ub/2.0_dp * real( nb*(nb-1),dp ) * abs(f(na,nb,i,j))**2
                  c2d_TotEnergy = c2d_TotEnergy + uab * real( na*nb,dp ) * abs(f(na,nb,i,j))**2

               enddo
            enddo

         enddo
      enddo

      do j = 1, ubound(f,4)
         do i = 1, ubound(f,3)

            if( i < ubound(f,3) ) then

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

      do k = 1, ubound(f,5)
         do j = 1, ubound(f,4)
            do i = 1, ubound(f,3)

               ! calculate order parameters
               ord = OrderA( f, i, j, k )
               ArrayOrderA(i,j,k) = ord
               ord = OrderB( f, i, j, k )
               ArrayOrderB(i,j,k) = ord

               do nb = 0, ubound(f,2)
                  do na = 0, ubound(f,1)

                     c3d_TotEnergy = c3d_TotEnergy + ua/2.0_dp * real( na*(na-1),dp ) * abs(f(na,nb,i,j,k))**2
                     c3d_TotEnergy = c3d_TotEnergy + ub/2.0_dp * real( nb*(nb-1),dp ) * abs(f(na,nb,i,j,k))**2
                     c3d_TotEnergy = c3d_TotEnergy + uab * real( na*nb,dp ) * abs(f(na,nb,i,j,k))**2

                  enddo
               enddo

            enddo
         enddo
      enddo

      do k = 1, ubound(f,5)
         do j = 1, ubound(f,4)
            do i = 1, ubound(f,3)

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

! | ----------------------------------- |
! | Chemical Potential                  |
! | ----------------------------------- | 

! | helper functions
! | --------------------------------------------------
   function c1d_xiA( f, i )

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:,:)
      integer, intent(in) :: i
      complex( kind = dp ) :: c1d_xiA
      integer na, nb

      c1d_xiA = (0.0_dp, 0.0_dp)
      do nb = 0, ubound(f,2)
         do na = 1, ubound(f,1)
            c1d_xiA = c1d_xiA + real(na,dp)*sqrt(real(na,dp)) * conjg(f(na-1,nb,i))*f(na,nb,i)
         enddo
      enddo

   end function

   function c1d_xiB( f, i )

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:,:)
      integer, intent(in) :: i
      complex( kind = dp ) :: c1d_xiB
      integer na, nb

      c1d_xiB = (0.0_dp, 0.0_dp)
      do nb = 1, ubound(f,2)
         do na = 0, ubound(f,1)
            c1d_xiB = c1d_xiB + real(nb,dp)*sqrt(real(nb,dp)) * conjg(f(na,nb-1,i))*f(na,nb,i)
         enddo
      enddo

   end function

   function c1d_etaA( f, i )

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:,:)
      integer, intent(in) :: i
      complex( kind = dp ) :: c1d_etaA
      integer na, nb

      c1d_etaA = (0.0_dp, 0.0_dp)
      do nb = 1, ubound(f,2)
         do na = 0, ubound(f,1)
            c1d_etaA = c1d_etaA + real(na,dp)*sqrt(real(nb,dp)) * conjg(f(na,nb-1,i))*f(na,nb,i)
         enddo
      enddo

   end function

   function c1d_etaB( f, i )

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:,:)
      integer, intent(in) :: i
      complex( kind = dp ) :: c1d_etaB
      integer na, nb

      c1d_etaB = (0.0_dp, 0.0_dp)
      do nb = 0, ubound(f,2)
         do na = 1, ubound(f,1)
            c1d_etaB = c1d_etaB + real(nb,dp)*sqrt(real(na,dp)) * conjg(f(na-1,nb,i))*f(na,nb,i)
         enddo
      enddo

   end function

! | --------------------------------------------------
   function c2d_xiA( f, i, j )

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:,:,:)
      integer, intent(in) :: i, j
      complex( kind = dp ) :: c2d_xiA
      integer na, nb

      c2d_xiA = (0.0_dp, 0.0_dp)
      do nb = 0, ubound(f,2)
         do na = 1, ubound(f,1)
            c2d_xiA = c2d_xiA + real(na,dp)*sqrt(real(na,dp)) * conjg(f(na-1,nb,i,j))*f(na,nb,i,j)
         enddo
      enddo

   end function

   function c2d_xiB( f, i, j )

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:,:,:)
      integer, intent(in) :: i, j
      complex( kind = dp ) :: c2d_xiB
      integer na, nb

      c2d_xiB = (0.0_dp, 0.0_dp)
      do nb = 1, ubound(f,2)
         do na = 0, ubound(f,1)
            c2d_xiB = c2d_xiB + real(nb,dp)*sqrt(real(nb,dp)) * conjg(f(na,nb-1,i,j))*f(na,nb,i,j)
         enddo
      enddo

   end function

   function c2d_etaA( f, i, j )

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:,:,:)
      integer, intent(in) :: i, j
      complex( kind = dp ) :: c2d_etaA
      integer na, nb

      c2d_etaA = (0.0_dp, 0.0_dp)
      do nb = 1, ubound(f,2)
         do na = 0, ubound(f,1)
            c2d_etaA = c2d_etaA + real(na,dp)*sqrt(real(nb,dp)) * conjg(f(na,nb-1,i,j))*f(na,nb,i,j)
         enddo
      enddo

   end function

   function c2d_etaB( f, i, j )

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:,:,:)
      integer, intent(in) :: i, j
      complex( kind = dp ) :: c2d_etaB
      integer na, nb

      c2d_etaB = (0.0_dp, 0.0_dp)
      do nb = 0, ubound(f,2)
         do na = 1, ubound(f,1)
            c2d_etaB = c2d_etaB + real(nb,dp)*sqrt(real(na,dp)) * conjg(f(na-1,nb,i,j))*f(na,nb,i,j)
         enddo
      enddo

   end function

! | --------------------------------------------------
   function c3d_xiA( f, i, j, k )

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:,:,:,:)
      integer, intent(in) :: i, j, k
      complex( kind = dp ) :: c3d_xiA
      integer na, nb

      c3d_xiA = (0.0_dp, 0.0_dp)
      do nb = 0, ubound(f,2)
         do na = 1, ubound(f,1)
            c3d_xiA = c3d_xiA + real(na,dp)*sqrt(real(na,dp)) * conjg(f(na-1,nb,i,j,k))*f(na,nb,i,j,k)
         enddo
      enddo

   end function

   function c3d_xiB( f, i, j, k )

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:,:,:,:)
      integer, intent(in) :: i, j, k
      complex( kind = dp ) :: c3d_xiB
      integer na, nb

      c3d_xiB = (0.0_dp, 0.0_dp)
      do nb = 1, ubound(f,2)
         do na = 0, ubound(f,1)
            c3d_xiB = c3d_xiB + real(nb,dp)*sqrt(real(nb,dp)) * conjg(f(na,nb-1,i,j,k))*f(na,nb,i,j,k)
         enddo
      enddo

   end function

   function c3d_etaA( f, i, j, k )

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:,:,:,:)
      integer, intent(in) :: i, j, k
      complex( kind = dp ) :: c3d_etaA
      integer na, nb

      c3d_etaA = (0.0_dp, 0.0_dp)
      do nb = 1, ubound(f,2)
         do na = 0, ubound(f,1)
            c3d_etaA = c3d_etaA + real(na,dp)*sqrt(real(nb,dp)) * conjg(f(na,nb-1,i,j,k))*f(na,nb,i,j,k)
         enddo
      enddo

   end function

   function c3d_etaB( f, i, j, k )

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:,:,:,:)
      integer, intent(in) :: i, j, k
      complex( kind = dp ) :: c3d_etaB
      integer na, nb

      c3d_etaB = (0.0_dp, 0.0_dp)
      do nb = 0, ubound(f,2)
         do na = 1, ubound(f,1)
            c3d_etaB = c3d_etaB + real(nb,dp)*sqrt(real(na,dp)) * conjg(f(na-1,nb,i,j,k))*f(na,nb,i,j,k)
         enddo
      enddo

   end function
! | --------------------------------------------------
   function c1d_ChemPot( f, ja, jb, ua, ub, uab )

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:,:)
      real( kind = dp ), intent(in) :: ua, ub, uab, ja, jb
      real( kind = dp ), dimension(2) :: c1d_ChemPot

      integer :: i ,na, nb
      real( kind = dp ) :: nna,nnb,na2,nb2,nab,na3,nb3,nanb2,na2nb
      complex( kind = dp ) :: ord,eta,xi,norda,nordb
      real( kind = dp ) :: coll_flucta,coll_fluctb,coll_fluctab
      real( kind = dp ) :: coeff_A, coeff_B, det
      complex( kind = dp ), allocatable :: ArrayOrderA(:), ArrayOrderB(:)
      complex( kind = dp ), allocatable :: ArrayXiA(:), ArrayXiB(:)
      complex( kind = dp ), allocatable :: ArrayEtaA(:), ArrayEtaB(:)

      allocate( ArrayOrderA( size(f,3) ), ArrayOrderB( size(f,3) ) )
      allocate( ArrayXiA( size(f,3) ), ArrayXiB( size(f,3) ) )
      allocate( ArrayEtaA( size(f,3) ), ArrayEtaB( size(f,3) ) )

      coll_flucta = 0.0_dp; coll_fluctb = 0.0_dp; coll_fluctab = 0.0_dp

      do i = 1, ubound(f,3)

         ord = c1d_OrderA( f, i )
         ArrayOrderA(i) = ord
         ord = c1d_OrderB( f, i )
         ArrayOrderB(i) = ord

         xi  = c1d_xiA( f, i )
         ArrayXiA(i) = xi
         xi  = c1d_xiB( f, i )
         ArrayXiB(i) = xi

         eta = c1d_etaA( f, i )
         ArrayEtaA(i) = eta
         eta = c1d_etaB( f, i )
         ArrayEtaB(i) = eta

      enddo

      do i = 1, ubound(f,3)

         ! nearest neighbour
         if( i == 1 ) then
            norda = ArrayOrderA(i+1) + ArrayOrderA( ubound(f,3) )
            nordb = ArrayOrderB(i+1) + ArrayOrderB( ubound(f,3) )
         elseif( i == ubound(f,3) ) then
            norda = ArrayOrderA(1) + ArrayOrderA( i-1 )
            nordb = ArrayOrderB(1) + ArrayOrderB( i-1 )
         else
            norda = ArrayOrderA(i+1) + ArrayOrderA(i-1)
            nordb = ArrayOrderB(i+1) + ArrayOrderB(i-1)
         endif

         nna = 0.0_dp; nnb = 0.0_dp; na2 = 0.0_dp; nb2 = 0.0_dp; nab = 0.0_dp
         na3 = 0.0_dp; nb3 = 0.0_dp; nanb2 = 0.0_dp; na2nb = 0.0_dp

         do nb = 0, ubound(f,2)
            do na = 0, ubound(f,1)

               ! density correlations
               nna  = nna + real(na,dp) * abs(f(na,nb,i))**2
               nnb  = nnb + real(nb,dp) * abs(f(na,nb,i))**2
               na2 = na2 + real(na,dp)**2 * abs(f(na,nb,i))**2
               nb2 = nb2 + real(nb,dp)**2 * abs(f(na,nb,i))**2
               nab = nab + real(na*nb,dp) * abs(f(na,nb,i))**2
               na3 = na3 + real(na,dp)**3 * abs(f(na,nb,i))**2
               nb3 = nb3 + real(nb,dp)**3 * abs(f(na,nb,i))**2
               nanb2 = nanb2 + real(na*nb*nb,dp) * abs(f(na,nb,i))**2
               na2nb = na2nb + real(na*na*nb,dp) * abs(f(na,nb,i))**2

             enddo
         enddo

         ! collect everything
         coll_flucta = coll_flucta + (na2 - nna**2)
         coll_fluctb = coll_fluctb + (nb2 - nnb**2)
         coll_fluctab = coll_fluctab + (nab - na*nb)

         coeff_A = coeff_A + 0.5_dp*ua*(na3 - na2*nna + nna**2 - na2) +&
                             0.5_dp*ub*(nanb2 - nna*nb2 + nna*nnb - nab) +&
                             uab*(na2nb - nna*nab) +&
                             2.0_dp*jb*real( conjg(nordb)*(nna*ArrayOrderB(i) - ArrayEtaA(i)) ) +&
                             ja*real( conjg(norda)*(nna*ArrayOrderA(i) + ArrayOrderA(i) - 2.0_dp*ArrayXiA(i)) )

         coeff_B = coeff_B + 0.5_dp*ua*(na2nb - nnb*na2 + nna*nnb - nab) +&
                             0.5_dp*ub*(nb3 - nb2*nnb + nnb**2 - nb2) +&
                             uab*(nanb2 - nnb*nab) +&
                             2.0_dp*ja*real( conjg(norda)*(nnb*ArrayOrderA(i) - ArrayEtaB(i)) ) +&
                             jb*real( conjg(nordb)*(nnb*ArrayOrderB(i) + ArrayOrderB(i) - 2.0_dp*ArrayXiB(i)) )
      enddo

      det = coll_flucta*coll_fluctb - coll_fluctab**2
      if (det .ne. 0.0_dp) then
         c1d_ChemPot(1) = (coeff_A * coll_fluctb - coeff_B * coll_fluctab)/det
         c1d_ChemPot(2) = (coeff_B * coll_flucta - coeff_A * coll_fluctab)/det
      else
         c1d_ChemPot(:) = 0.0_dp
         print *, "Cannot determine chemical potential!"
      endif

   end function

   function c2d_ChemPot( f, ja, jb, ua, ub, uab )

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:,:,:)
      real( kind = dp ), intent(in) :: ua, ub, uab, ja, jb
      real( kind = dp ), dimension(2) :: c2d_ChemPot

      integer :: i, j, na, nb
      real( kind = dp ) :: nna,nnb,na2,nb2,nab,na3,nb3,nanb2,na2nb
      complex( kind = dp ) :: ord,eta,xi,norda,nordb
      real( kind = dp ) :: coll_flucta,coll_fluctb,coll_fluctab
      real( kind = dp ) :: coeff_A, coeff_B, det
      complex( kind = dp ), allocatable :: ArrayOrderA(:,:), ArrayOrderB(:,:)
      complex( kind = dp ), allocatable :: ArrayXiA(:,:), ArrayXiB(:,:)
      complex( kind = dp ), allocatable :: ArrayEtaA(:,:), ArrayEtaB(:,:)

      allocate( ArrayOrderA( size(f,3), size(f,4) ) )
      allocate( ArrayOrderB( size(f,3), size(f,4) ) )
      allocate( ArrayXiA( size(f,3), size(f,4) ) )
      allocate( ArrayXiB( size(f,3), size(f,4) ) )
      allocate( ArrayEtaA( size(f,3), size(f,4) ) )
      allocate( ArrayEtaB( size(f,3), size(f,4) ) )

      coll_flucta = 0.0_dp; coll_fluctb = 0.0_dp; coll_fluctab = 0.0_dp

      do j = 1, ubound(f,4)
      do i = 1, ubound(f,3)

         ord = c2d_OrderA( f, i, j )
         ArrayOrderA(i,j) = ord
         ord = c2d_OrderB( f, i, j )
         ArrayOrderB(i,j) = ord

         xi  = c2d_xiA( f, i, j )
         ArrayXiA(i,j) = xi
         xi  = c2d_xiB( f, i, j )
         ArrayXiB(i,j) = xi

         eta = c2d_etaA( f, i, j )
         ArrayEtaA(i,j) = eta
         eta = c2d_etaB( f, i, j )
         ArrayEtaB(i,j) = eta

      enddo ! i
      enddo ! j

      do j = 1, ubound(f,4)
      do i = 1, ubound(f,3)

         ! nearest neighbour
         if( i == 1 ) then
            norda = ArrayOrderA(i+1,j) + ArrayOrderA( ubound(f,3),j )
            nordb = ArrayOrderB(i+1,j) + ArrayOrderB( ubound(f,3),j )
         elseif( i == ubound(f,3) ) then
            norda = ArrayOrderA(1,j) + ArrayOrderA( i-1,j )
            nordb = ArrayOrderB(1,j) + ArrayOrderB( i-1,j )
         else
            norda = ArrayOrderA(i+1,j) + ArrayOrderA(i-1,j)
            nordb = ArrayOrderB(i+1,j) + ArrayOrderB(i-1,j)
         endif

         if( j == 1 ) then
            norda = norda + ArrayOrderA(i,j+1) + ArrayOrderA( i,ubound(f,4) )
            nordb = nordb + ArrayOrderB(i,j+1) + ArrayOrderB( i,ubound(f,4) )
         elseif( j == ubound(f,4) ) then
            norda = norda + ArrayOrderA(i,1) + ArrayOrderA( i,j-1 )
            nordb = nordb + ArrayOrderB(i,1) + ArrayOrderB( i,j-1 )
         else
            norda = norda + ArrayOrderA(i,j+1) + ArrayOrderA(i,j-1)
            nordb = nordb + ArrayOrderB(i,j+1) + ArrayOrderB(i,j-1)
         endif

         nna = 0.0_dp; nnb = 0.0_dp; na2 = 0.0_dp; nb2 = 0.0_dp; nab = 0.0_dp
         na3 = 0.0_dp; nb3 = 0.0_dp; nanb2 = 0.0_dp; na2nb = 0.0_dp

         do nb = 0, ubound(f,2)
            do na = 0, ubound(f,1)

               ! density correlations
               nna  = nna + real(na,dp) * abs(f(na,nb,i,j))**2
               nnb  = nnb + real(nb,dp) * abs(f(na,nb,i,j))**2
               na2 = na2 + real(na,dp)**2 * abs(f(na,nb,i,j))**2
               nb2 = nb2 + real(nb,dp)**2 * abs(f(na,nb,i,j))**2
               nab = nab + real(na*nb,dp) * abs(f(na,nb,i,j))**2
               na3 = na3 + real(na,dp)**3 * abs(f(na,nb,i,j))**2
               nb3 = nb3 + real(nb,dp)**3 * abs(f(na,nb,i,j))**2
               nanb2 = nanb2 + real(na*nb*nb,dp) * abs(f(na,nb,i,j))**2
               na2nb = na2nb + real(na*na*nb,dp) * abs(f(na,nb,i,j))**2

             enddo
         enddo

         ! collect everything
         coll_flucta = coll_flucta + (na2 - nna**2)
         coll_fluctb = coll_fluctb + (nb2 - nnb**2)
         coll_fluctab = coll_fluctab + (nab - na*nb)

         coeff_A = coeff_A + 0.5_dp*ua*(na3 - na2*nna + nna**2 - na2) +&
                             0.5_dp*ub*(nanb2 - nna*nb2 + nna*nnb - nab) +&
                             uab*(na2nb - nna*nab) +&
                             2.0_dp*jb*real( conjg(nordb)*(nna*ArrayOrderB(i,j) - ArrayEtaA(i,j)) ) +&  
                             ja*real( conjg(norda)*(nna*ArrayOrderA(i,j) + ArrayOrderA(i,j) - 2.0_dp*ArrayXiA(i,j)) )

         coeff_B = coeff_B + 0.5_dp*ua*(na2nb - nnb*na2 + nna*nnb - nab) +&
                             0.5_dp*ub*(nb3 - nb2*nnb + nnb**2 - nb2) +&
                             uab*(nanb2 - nnb*nab) +&
                             2.0_dp*ja*real( conjg(norda)*(nnb*ArrayOrderA(i,j) - ArrayEtaB(i,j)) ) +&
                             jb*real( conjg(nordb)*(nnb*ArrayOrderB(i,j) + ArrayOrderB(i,j) - 2.0_dp*ArrayXiB(i,j)) )
      enddo ! i
      enddo ! j

      det = coll_flucta*coll_fluctb - coll_fluctab**2
      if (det .ne. 0.0_dp) then
         c2d_ChemPot(1) = (coeff_A * coll_fluctb - coeff_B * coll_fluctab)/det
         c2d_ChemPot(2) = (coeff_B * coll_flucta - coeff_A * coll_fluctab)/det
      else
         c2d_ChemPot(:) = 0.0_dp
         print *, "Cannot determine chemical potential!"
      endif

   end function

   function c3d_ChemPot( f, ja, jb, ua, ub, uab )

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:,:,:,:)
      real( kind = dp ), intent(in) :: ua, ub, uab, ja, jb
      real( kind = dp ), dimension(2) :: c3d_ChemPot

      integer :: i, j, k, na, nb
      real( kind = dp ) :: nna,nnb,na2,nb2,nab,na3,nb3,nanb2,na2nb
      complex( kind = dp ) :: ord,eta,xi,norda,nordb
      real( kind = dp ) :: coll_flucta,coll_fluctb,coll_fluctab
      real( kind = dp ) :: coeff_A, coeff_B, det
      complex( kind = dp ), allocatable :: ArrayOrderA(:,:,:), ArrayOrderB(:,:,:)
      complex( kind = dp ), allocatable :: ArrayXiA(:,:,:), ArrayXiB(:,:,:)
      complex( kind = dp ), allocatable :: ArrayEtaA(:,:,:), ArrayEtaB(:,:,:)

      allocate( ArrayOrderA( size(f,3), size(f,4), size(f,5) ) )
      allocate( ArrayOrderB( size(f,3), size(f,4), size(f,5) ) )
      allocate( ArrayXiA( size(f,3), size(f,4), size(f,5) ) )
      allocate( ArrayXiB( size(f,3), size(f,4), size(f,5) ) )
      allocate( ArrayEtaA( size(f,3), size(f,4), size(f,5) ) )
      allocate( ArrayEtaB( size(f,3), size(f,4), size(f,5) ) )

      coll_flucta = 0.0_dp; coll_fluctb = 0.0_dp; coll_fluctab = 0.0_dp

      do k = 1, ubound(f,5)
      do j = 1, ubound(f,4)
      do i = 1, ubound(f,3)

         ord = c3d_OrderA( f, i, j, k )
         ArrayOrderA(i,j,k) = ord
         ord = c3d_OrderB( f, i, j, k )
         ArrayOrderB(i,j,k) = ord

         xi  = c3d_xiA( f, i, j, k )
         ArrayXiA(i,j,k) = xi
         xi  = c3d_xiB( f, i, j, k )
         ArrayXiB(i,j,k) = xi

         eta = c3d_etaA( f, i, j, k )
         ArrayEtaA(i,j,k) = eta
         eta = c3d_etaB( f, i, j, k )
         ArrayEtaB(i,j,k) = eta

      enddo ! i
      enddo ! j
      enddo ! k

      do k = 1, ubound(f,5)
      do j = 1, ubound(f,4)
      do i = 1, ubound(f,3)

         ! nearest neighbour
         if( i == 1 ) then
            norda = ArrayOrderA(i+1,j,k) + ArrayOrderA( ubound(f,3),j,k )
            nordb = ArrayOrderB(i+1,j,k) + ArrayOrderB( ubound(f,3),j,k )
         elseif( i == ubound(f,3) ) then
            norda = ArrayOrderA(1,j,k) + ArrayOrderA( i-1,j,k )
            nordb = ArrayOrderB(1,j,k) + ArrayOrderB( i-1,j,k )
         else
            norda = ArrayOrderA(i+1,j,k) + ArrayOrderA(i-1,j,k)
            nordb = ArrayOrderB(i+1,j,k) + ArrayOrderB(i-1,j,k)
         endif

         if( j == 1 ) then
            norda = norda + ArrayOrderA(i,j+1,k) + ArrayOrderA( i,ubound(f,4),k )
            nordb = nordb + ArrayOrderB(i,j+1,k) + ArrayOrderB( i,ubound(f,4),k )
         elseif( j == ubound(f,4) ) then
            norda = norda + ArrayOrderA(i,1,k) + ArrayOrderA( i,j-1,k )
            nordb = nordb + ArrayOrderB(i,1,k) + ArrayOrderB( i,j-1,k )
         else
            norda = norda + ArrayOrderA(i,j+1,k) + ArrayOrderA(i,j-1,k)
            nordb = nordb + ArrayOrderB(i,j+1,k) + ArrayOrderB(i,j-1,k)
         endif

         if( k == 1 ) then
            norda = norda + ArrayOrderA(i,j,k+1) + ArrayOrderA( i,j,ubound(f,5) )
            nordb = nordb + ArrayOrderB(i,j,k+1) + ArrayOrderB( i,j,ubound(f,5) )
         elseif( k == ubound(f,5) ) then
            norda = norda + ArrayOrderA(i,j,1) + ArrayOrderA( i,j,k-1 )
            nordb = nordb + ArrayOrderB(i,j,1) + ArrayOrderB( i,j,k-1 )
         else
            norda = norda + ArrayOrderA(i,j,k+1) + ArrayOrderA(i,j,k-1)
            nordb = nordb + ArrayOrderB(i,j,k+1) + ArrayOrderB(i,j,k-1)
         endif

         nna = 0.0_dp; nnb = 0.0_dp; na2 = 0.0_dp; nb2 = 0.0_dp; nab = 0.0_dp
         na3 = 0.0_dp; nb3 = 0.0_dp; nanb2 = 0.0_dp; na2nb = 0.0_dp

         do nb = 0, ubound(f,2)
            do na = 0, ubound(f,1)

               ! density correlations
               nna  = nna + real(na,dp) * abs(f(na,nb,i,j,k))**2
               nnb  = nnb + real(nb,dp) * abs(f(na,nb,i,j,k))**2
               na2 = na2 + real(na,dp)**2 * abs(f(na,nb,i,j,k))**2
               nb2 = nb2 + real(nb,dp)**2 * abs(f(na,nb,i,j,k))**2
               nab = nab + real(na*nb,dp) * abs(f(na,nb,i,j,k))**2
               na3 = na3 + real(na,dp)**3 * abs(f(na,nb,i,j,k))**2
               nb3 = nb3 + real(nb,dp)**3 * abs(f(na,nb,i,j,k))**2
               nanb2 = nanb2 + real(na*nb*nb,dp) * abs(f(na,nb,i,j,k))**2
               na2nb = na2nb + real(na*na*nb,dp) * abs(f(na,nb,i,j,k))**2

             enddo
         enddo

         ! collect everything
         coll_flucta = coll_flucta + (na2 - nna**2)
         coll_fluctb = coll_fluctb + (nb2 - nnb**2)
         coll_fluctab = coll_fluctab + (nab - na*nb)

         coeff_A = coeff_A + 0.5_dp*ua*(na3 - na2*nna + nna**2 - na2) +&
                             0.5_dp*ub*(nanb2 - nna*nb2 + nna*nnb - nab) +&
                             uab*(na2nb - nna*nab) +&
                             2.0_dp*jb*real( conjg(nordb)*(nna*ArrayOrderB(i,j,k) - ArrayEtaA(i,j,k)) ) +&
                             ja*real( conjg(norda)*(nna*ArrayOrderA(i,j,k) + ArrayOrderA(i,j,k) - 2.0_dp*ArrayXiA(i,j,k)) )

         coeff_B = coeff_B + 0.5_dp*ua*(na2nb - nnb*na2 + nna*nnb - nab) +&
                             0.5_dp*ub*(nb3 - nb2*nnb + nnb**2 - nb2) +&
                             uab*(nanb2 - nnb*nab) +&
                             2.0_dp*ja*real( conjg(norda)*(nnb*ArrayOrderA(i,j,k) - ArrayEtaB(i,j,k)) ) +&
                             jb*real( conjg(nordb)*(nnb*ArrayOrderB(i,j,k) + ArrayOrderB(i,j,k) - 2.0_dp*ArrayXiB(i,j,k)) )
      enddo ! i
      enddo ! j
      enddo ! k

      det = coll_flucta*coll_fluctb - coll_fluctab**2
      if (det .ne. 0.0_dp) then
         c3d_ChemPot(1) = (coeff_A * coll_fluctb - coeff_B * coll_fluctab)/det
         c3d_ChemPot(2) = (coeff_B * coll_flucta - coeff_A * coll_fluctab)/det
      else
         c3d_ChemPot(:) = 0.0_dp
         print *, "Cannot determine chemical potential!"
      endif

   end function

end module
