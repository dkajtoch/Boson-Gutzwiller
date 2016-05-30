module measures

   use parameters

   interface MeanA
      module procedure c1d_MeanA
   end interface

   interface MeanB
      module procedure c1d_MeanB
   end interface

   interface MeanAB
      module procedure c1d_MeanAB
   end interface

   interface VarA
      module procedure c1d_VarA
   end interface 

   interface VarB
      module procedure c1d_VarB
   end interface

   interface norm
      module procedure c1d_norm
   end interface

   interface normalize
      module procedure c1d_normalize
   end interface

   interface OrderA
      module procedure c1d_OrderA
   end interface

   interface OrderB
      module procedure c1d_OrderB
   end interface

   interface TotEnergy
      module procedure c1d_TotEnergy
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
      do na = 0, ubound(f,2)
         do nb = 0, ubound(f,3)
            c1d_MeanA = c1d_MeanA + real( na,dp ) * abs(f(i,na,nb))**2
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
      do na = 0, ubound(f,2)
         do nb = 0, ubound(f,3)
            c1d_MeanB = c1d_MeanB + real( nb,dp ) * abs(f(i,na,nb))**2
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
      do na = 0, ubound(f,2)
         do nb = 0, ubound(f,3)
            c1d_MeanAB = c1d_MeanAB + real( na*nb,dp ) * abs(f(i,na,nb))**2
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

      do na = 0, ubound(f,2)
         do nb = 0, ubound(f,3)
            suma   = suma + real( na,dp ) * abs(f(i,na,nb))**2
            c1d_VarA = c1d_VarA + real( na*na,dp ) * abs(f(i,na,nb))**2
         enddo
      enddo

      c1d_VarA = c1d_VarA - suma**2

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

      do na = 0, ubound(f,2)
         do nb = 0, ubound(f,3)
            suma   = suma + real( nb,dp ) * abs(f(i,na,nb))**2
            c1d_VarB = c1d_VarB + real( nb*nb,dp ) * abs(f(i,na,nb))**2
         enddo
      enddo

      c1d_VarB = c1d_VarB - suma**2

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

      do na = 0, ubound(f,2)
         do nb = 0, ubound(f,3)
            c1d_norm = c1d_norm + abs(f(i,na,nb))**2
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

      do i = 1, ubound(f,1)
         nout = norm( f, i )
         f(i,:,:) = f(i,:,:)*1.0_dp/SQRT( nout )
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

      do na = 1, ubound(f,2)
         do nb = 0, ubound(f,3)
            c1d_orderA = c1d_orderA + SQRT( real( na,dp ) ) * conjg(f(i,na-1,nb)) * f(i,na,nb)
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

      do na = 0, ubound(f,2)
         do nb = 1, ubound(f,3)
            c1d_orderB = c1d_orderB + SQRT( real( nb,dp ) ) * conjg(f(i,na,nb-1)) * f(i,na,nb)
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

      allocate( ArrayOrderA( size(f,1) ), ArrayOrderB( size(f,1) ) )

      c1d_TotEnergy = 0.0_dp

      do i = 1, ubound(f,1)

         ! calculate order parameters
         ord = OrderA( f, i )
         ArrayOrderA(i) = ord
         ord = OrderB( f, i )
         ArrayOrderB(i) = ord

         do na = 0, ubound(f,2)
            do nb = 0, ubound(f,3)

               c1d_TotEnergy = c1d_TotEnergy + ua/2.0_dp * real( na*(na-1),dp ) * abs(f(i,na,nb))**2
               c1d_TotEnergy = c1d_TotEnergy + ub/2.0_dp * real( nb*(nb-1),dp ) * abs(f(i,na,nb))**2
               c1d_TotEnergy = c1d_TotEnergy + uab * real( na*nb,dp ) * abs(f(i,na,nb))**2

            enddo
         enddo

      enddo

      do i = 1, ubound(f,1)

         c1d_TotEnergy = c1d_TotEnergy - 2.0_dp * ja * real( conjg(ArrayOrderA(i)) * ArrayOrderA(i+1) )
         c1d_TotEnergy = c1d_TotEnergy - 2.0_dp * jb * real( conjg(ArrayOrderB(i)) * ArrayOrderB(i+1) )

      enddo

      c1d_TotEnergy = c1d_TotEnergy - 2.0_dp * ja * real( conjg(ArrayOrderA( ubound(f,1) )) * ArrayOrderA(1) )
      c1d_TotEnergy = c1d_TotEnergy - 2.0_dp * jb * real( conjg(ArrayOrderB( ubound(f,1) )) * ArrayOrderB(1) )

   end function

end module
