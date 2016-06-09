module solvers

   use parameters
   use measures, only: normalize
   use fixed_mean, only: const_mean

   private
   public GroundState, GroundStateNC, InitUniform, InitUniformNC

   interface GroundState
      module procedure c1d_GroundState
      module procedure c2d_GroundState
      module procedure c3d_GroundState
   end interface

   interface GroundStateNC
      module procedure c1d_GroundStateNC
      module procedure c2d_GroundStateNC
      module procedure c3d_GroundStateNC
   end interface

   interface InitUniform
      module procedure c1d_InitUniform
      module procedure c2d_InitUniform
      module procedure c3d_InitUniform
   end interface

   interface InitUniformNC
      module procedure c1d_InitUniformNC
      module procedure c2d_InitUniformNC
      module procedure c3d_InitUniformNC
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

      allocate( ArrayOrderA( size(f,3) ), ArrayOrderB( size(f,3) ) )

      ! collect order parameters
      do i = 1, ubound(f,3)
         
         ArrayOrderA(i) = 0.0_dp*re
         ArrayOrderB(i) = 0.0_dp*re

         do na = 0, ubound(f,1)
            do nb = 0, ubound(f,2)
               
               if( na > 0 ) then
                  ArrayOrderA(i) = ArrayOrderA(i) + SQRT( real(na,dp) ) * conjg( f(na-1,nb,i) ) * f(na,nb,i)
               endif
               if( nb > 0 ) then
                  ArrayOrderB(i) = ArrayOrderB(i) + SQRT( real(nb,dp) ) * conjg( f(na,nb-1,i) ) * f(na,nb,i)
               endif

            enddo
         enddo

      enddo

      ! solve coupled system of equations
      do i = 1, ubound(f,3)

         ! calculate nearest-neighbour order parameter (periodic boundary conditions)
         if( i .eq. 1 ) then
            ordA = ArrayOrderA( ubound(f,3) ) + ArrayOrderA(i+1)
            ordB = ArrayOrderB( ubound(f,3) ) + ArrayOrderB(i+1)
         elseif( i .eq. ubound(f,3) ) then
            ordA = ArrayOrderA(1) + ArrayOrderA(i-1)
            ordB = ArrayOrderB(1) + ArrayOrderB(i-1)
         else
            ordA = ArrayOrderA(i+1) + ArrayOrderA(i-1)
            ordB = ArrayOrderB(i+1) + ArrayOrderB(i-1)
         endif
         
         do na = 0, ubound(f,1)
            do nb = 0, ubound(f,2)

               ! chop
               if( abs( real( f(na,nb,i) ) ) < chopCutoff ) then
                  f(na,nb,i) = ( f(na,nb,i) - conjg(f(na,nb,i)) )/2.0_dp
               endif
               if( abs( aimag( f(na,nb,i) ) ) < chopCutoff ) then
                  f(na,nb,i) = ( f(na,nb,i) + conjg(f(na,nb,i)) )/2.0_dp
               endif


               k(na,nb,i) = mlt * ( real(na,dp) * ( mua - ua/2.0_dp * (real(na,dp) - 1.0_dp) ) + &
                                    real(nb,dp) * ( mub - ub/2.0_dp * (real(nb,dp) - 1.0_dp) ) - &
                                    uab * real(na,dp) * real(nb,dp) ) * f(na,nb,i)

               if( na < ubound(f,1) ) then
                  k(na,nb,i) = k(na,nb,i) + mlt * ja * ordA * SQRT( real(na+1,dp) ) * f(na+1,nb,i)
               endif
               if( na > 0 ) then
                  k(na,nb,i) = k(na,nb,i) + mlt * ja * ordA * SQRT( real(na,dp) ) * f(na-1,nb,i)
               endif

               if( nb < ubound(f,2) ) then
                  k(na,nb,i) = k(na,nb,i) + mlt * jb * ordB * SQRT( real(nb+1,dp) ) * f(na,nb+1,i)
               endif
               if( nb > 0 ) then
                  k(na,nb,i) = k(na,nb,i) + mlt * jb * ordB * SQRT( real(nb,dp) ) * f(na,nb-1,i)
               endif

            enddo
         enddo

      enddo

   end subroutine

   subroutine c2d_solve( f, k, ja, jb, ua, ub, uab, mua, mub, mode )

      implicit none
      complex( kind = dp), intent(inout), allocatable :: f(:,:,:,:)
      complex( kind = dp ), intent(inout), allocatable :: k(:,:,:,:)
      real( kind = dp ), intent(in) :: ja, jb, ua, ub, uab, mua, mub
      character(len=2), intent(in) :: mode

      ! local variables
      integer i, j, na, nb
      complex( kind = dp ), allocatable :: ArrayOrderA(:,:), ArrayOrderB(:,:)
      complex( kind = dp ) :: ordA, ordB, mlt

      if (mode .eq. 'EV') then
         mlt = im
      else if (mode .eq. 'GS') then
         mlt = re
      else
         mlt = (0.0_dp, 0.0_dp) 
         write(*,*)  "Wrong mode in SOLVE subroutine!"
      endif

      allocate( ArrayOrderA( size(f,3), size(f,4) ), ArrayOrderB( size(f,3), size(f,4) ) )

      ! collect order parameters
      do i = 1, ubound(f,3)
         do j = 1, ubound(f,4)
         
            ArrayOrderA(i,j) = 0.0_dp*re
            ArrayOrderB(i,j) = 0.0_dp*re

            do na = 0, ubound(f,1)
               do nb = 0, ubound(f,2)
               
                  if( na > 0 ) then
                     ArrayOrderA(i,j) = ArrayOrderA(i,j) + SQRT( real(na,dp) ) * conjg( f(na-1,nb,i,j) ) * f(na,nb,i,j)
                  endif
                  if( nb > 0 ) then
                     ArrayOrderB(i,j) = ArrayOrderB(i,j) + SQRT( real(nb,dp) ) * conjg( f(na,nb-1,i,j) ) * f(na,nb,i,j)
                  endif

               enddo
            enddo

         enddo
      enddo

      ! solve coupled system of equations
      do i = 1, ubound(f,3)
         do j = 1, ubound(f,4)

            ! calculate nearest-neighbour order parameter (periodic boundary conditions)
            if( i == 1 ) then

               if( j == 1 ) then

                  ordA = ArrayOrderA( ubound(f,3), j ) + ArrayOrderA(i+1, j) + &
                         ArrayOrderA( i, ubound(f,4) ) + ArrayOrderA(i, j+1)

                  ordB = ArrayOrderB( ubound(f,3), j ) + ArrayOrderB(i+1, j) + &
                         ArrayOrderB( i, ubound(f,4) ) + ArrayOrderB(i, j+1)
               
               elseif( j == ubound(f,4) ) then

                  ordA = ArrayOrderA( ubound(f,3), j ) + ArrayOrderA(i+1, j) + &
                         ArrayOrderA( i, 1 ) + ArrayOrderA(i, j-1)

                  ordB = ArrayOrderB( ubound(f,3), j ) + ArrayOrderB(i+1, j) + &
                         ArrayOrderB( i, 1 ) + ArrayOrderB(i, j-1)
               else

                  ordA = ArrayOrderA( ubound(f,3), j ) + ArrayOrderA(i+1, j) + &
                         ArrayOrderA( i, j+1 ) + ArrayOrderA(i, j-1)

                  ordB = ArrayOrderB( ubound(f,3), j ) + ArrayOrderB(i+1, j) + &
                         ArrayOrderB( i, j+1 ) + ArrayOrderB(i, j-1)

               endif

            elseif( i == ubound(f,3) ) then

               if( j == 1 ) then

                  ordA = ArrayOrderA( 1, j ) + ArrayOrderA(i-1, j) + &
                         ArrayOrderA( i, ubound(f,4) ) + ArrayOrderA(i, j+1)

                  ordB = ArrayOrderB( 1, j ) + ArrayOrderB(i-1, j) + &
                         ArrayOrderB( i, ubound(f,4) ) + ArrayOrderB(i, j+1)
               
               elseif( j == ubound(f,4) ) then

                  ordA = ArrayOrderA( 1, j ) + ArrayOrderA(i-1, j) + &
                         ArrayOrderA( i, 1 ) + ArrayOrderA(i, j-1)

                  ordB = ArrayOrderB( 1, j ) + ArrayOrderB(i-1, j) + &
                         ArrayOrderB( i, 1 ) + ArrayOrderB(i, j-1)
               else

                  ordA = ArrayOrderA( 1, j ) + ArrayOrderA(i-1, j) + &
                         ArrayOrderA( i, j+1 ) + ArrayOrderA(i, j-1)

                  ordB = ArrayOrderB( 1, j ) + ArrayOrderB(i-1, j) + &
                         ArrayOrderB( i, j+1 ) + ArrayOrderB(i, j-1)

               endif

            else

               if( j == 1 ) then

                  ordA = ArrayOrderA( i+1, j ) + ArrayOrderA(i-1, j) + &
                         ArrayOrderA( i, ubound(f,4) ) + ArrayOrderA(i, j+1)

                  ordB = ArrayOrderB( 1, j ) + ArrayOrderB(i-1, j) + &
                         ArrayOrderB( i, ubound(f,4) ) + ArrayOrderB(i, j+1)
               
               elseif( j == ubound(f,4) ) then

                  ordA = ArrayOrderA( i+1, j ) + ArrayOrderA(i-1, j) + &
                         ArrayOrderA( i, 1 ) + ArrayOrderA(i, j-1)

                  ordB = ArrayOrderB( i+1, j ) + ArrayOrderB(i-1, j) + &
                         ArrayOrderB( i, 1 ) + ArrayOrderB(i, j-1)
               else

                  ordA = ArrayOrderA( i+1, j ) + ArrayOrderA(i-1, j) + &
                         ArrayOrderA( i, j+1 ) + ArrayOrderA(i, j-1)

                  ordB = ArrayOrderB( i+1, j ) + ArrayOrderB(i-1, j) + &
                         ArrayOrderB( i, j+1 ) + ArrayOrderB(i, j-1)

               endif

            endif
         
            do na = 0, ubound(f,1)
               do nb = 0, ubound(f,2)

                  ! chop
                  if( abs( real( f(na,nb,i,j) ) ) < chopCutoff ) then
                     f(na,nb,i,j) = ( f(na,nb,i,j) - conjg(f(na,nb,i,j)) )/2.0_dp
                  endif
                  if( abs( aimag( f(na,nb,i,j) ) ) < chopCutoff ) then
                     f(na,nb,i,j) = ( f(na,nb,i,j) + conjg(f(na,nb,i,j)) )/2.0_dp
                  endif

                  k(na,nb,i,j) = mlt * ( real(na,dp) * ( mua - ua/2.0_dp * (real(na,dp) - 1.0_dp) ) + &
                                         real(nb,dp) * ( mub - ub/2.0_dp * (real(nb,dp) - 1.0_dp) ) - &
                                         uab * real(na,dp) * real(nb,dp) ) * f(na,nb,i,j)

                  if( na < ubound(f,1) ) then
                     k(na,nb,i,j) = k(na,nb,i,j) + mlt * ja * ordA * SQRT( real(na+1,dp) ) * f(na+1,nb,i,j)
                  endif
                  if( na > 0 ) then
                     k(na,nb,i,j) = k(na,nb,i,j) + mlt * ja * ordA * SQRT( real(na,dp) ) * f(na-1,nb,i,j)
                  endif

                  if( nb < ubound(f,2) ) then
                     k(na,nb,i,j) = k(na,nb,i,j) + mlt * jb * ordB * SQRT( real(nb+1,dp) ) * f(na,nb+1,i,j)
                  endif
                  if( nb > 0 ) then
                     k(na,nb,i,j) = k(na,nb,i,j) + mlt * jb * ordB * SQRT( real(nb,dp) ) * f(na,nb-1,i,j)
                  endif

               enddo
            enddo

         enddo
      enddo

   end subroutine

   subroutine c3d_solve( f, kv, ja, jb, ua, ub, uab, mua, mub, mode )

      implicit none
      complex( kind = dp), intent(inout), allocatable :: f(:,:,:,:,:)
      complex( kind = dp ), intent(inout), allocatable :: kv(:,:,:,:,:)
      real( kind = dp ), intent(in) :: ja, jb, ua, ub, uab, mua, mub
      character(len=2), intent(in) :: mode

      ! local variables
      integer i, j, k, na, nb
      complex( kind = dp ), allocatable :: ArrayOrderA(:,:,:), ArrayOrderB(:,:,:)
      complex( kind = dp ) :: ordA, ordB, mlt

      if (mode .eq. 'EV') then
         mlt = im
      else if (mode .eq. 'GS') then
         mlt = re
      else
         mlt = (0.0_dp, 0.0_dp) 
         write(*,*)  "Wrong mode in SOLVE subroutine!"
      endif

      allocate( ArrayOrderA( size(f,3), size(f,4), size(f,5) ), ArrayOrderB( size(f,3), size(f,4), size(f,5) ) )

      ! collect order parameters
      do i = 1, ubound(f,3)
         do j = 1, ubound(f,4)
            do k = 1, ubound(f,5)
         
               ArrayOrderA(i,j,k) = 0.0_dp*re
               ArrayOrderB(i,j,k) = 0.0_dp*re

               do na = 0, ubound(f,1)
                  do nb = 0, ubound(f,2)
               
                     if( na > 0 ) then
                        ArrayOrderA(i,j,k) = ArrayOrderA(i,j,k) + SQRT( real(na,dp) ) * conjg( f(na-1,nb,i,j,k) ) * f(na,nb,i,j,k)
                     endif
                     if( nb > 0 ) then
                        ArrayOrderB(i,j,k) = ArrayOrderB(i,j,k) + SQRT( real(nb,dp) ) * conjg( f(na,nb-1,i,j,k) ) * f(na,nb,i,j,k)
                     endif

                  enddo
               enddo

            enddo
         enddo
      enddo

      ! solve coupled system of equations
      do i = 1, ubound(f,3)
         do j = 1, ubound(f,4)
            do k = 1, ubound(f,5)

               ! calculate nearest-neighbour order parameter (periodic boundary conditions)
               if( i == 1 ) then

                  if( j == 1 ) then

                     if( k == 1 ) then

                        ordA = ArrayOrderA( ubound(f,3), j, k ) + ArrayOrderA(i+1, j, k) + &
                               ArrayOrderA( i, ubound(f,4), k ) + ArrayOrderA(i, j+1, k) + &
                               ArrayOrderA( i, j, ubound(f,5) ) + ArrayOrderA(i, j, k+1)

                        ordB = ArrayOrderB( ubound(f,3), j, k ) + ArrayOrderB(i+1, j, k) + &
                               ArrayOrderB( i, ubound(f,4), k ) + ArrayOrderB(i, j+1, k) + &
                               ArrayOrderB( i, j, ubound(f,5) ) + ArrayOrderB(i, j, k+1)

                     elseif( k == ubound(f,5) ) then
                            
                        ordA = ArrayOrderA( ubound(f,3), j, k ) + ArrayOrderA(i+1, j, k) + &
                               ArrayOrderA( i, ubound(f,4), k ) + ArrayOrderA(i, j+1, k) + &
                               ArrayOrderA( i, j, 1 ) + ArrayOrderA(i, j, k-1)

                        ordB = ArrayOrderB( ubound(f,3), j, k ) + ArrayOrderB(i+1, j, k) + &
                               ArrayOrderB( i, ubound(f,4), k ) + ArrayOrderB(i, j+1, k) + &
                               ArrayOrderB( i, j, 1 ) + ArrayOrderB(i, j, k-1)

                     else

                        ordA = ArrayOrderA( ubound(f,3), j, k ) + ArrayOrderA(i+1, j, k) + &
                               ArrayOrderA( i, ubound(f,4), k ) + ArrayOrderA(i, j+1, k) + &
                               ArrayOrderA( i, j, k+1 ) + ArrayOrderA(i, j, k-1)

                        ordB = ArrayOrderB( ubound(f,3), j, k ) + ArrayOrderB(i+1, j, k) + &
                               ArrayOrderB( i, ubound(f,4), k ) + ArrayOrderB(i, j+1, k) + &
                               ArrayOrderB( i, j, k+1 ) + ArrayOrderB(i, j, k-1)

                     endif

                  elseif( j == ubound(f,4) ) then

                     if( k == 1 ) then

                        ordA = ArrayOrderA( ubound(f,3), j, k ) + ArrayOrderA(i+1, j, k) + &
                               ArrayOrderA( i, 1, k ) + ArrayOrderA(i, j-1, k) + &
                               ArrayOrderA( i, j, ubound(f,5) ) + ArrayOrderA(i, j, k+1)

                        ordB = ArrayOrderB( ubound(f,3), j, k ) + ArrayOrderB(i+1, j, k) + &
                               ArrayOrderB( i, 1, k ) + ArrayOrderB(i, j-1, k) + &
                               ArrayOrderB( i, j, ubound(f,5) ) + ArrayOrderB(i, j, k+1)

                     elseif( k == ubound(f,5) ) then
                            
                        ordA = ArrayOrderA( ubound(f,3), j, k ) + ArrayOrderA(i+1, j, k) + &
                               ArrayOrderA( i, 1, k ) + ArrayOrderA(i, j-1, k) + &
                               ArrayOrderA( i, j, 1 ) + ArrayOrderA(i, j, k-1)

                        ordB = ArrayOrderB( ubound(f,3), j, k ) + ArrayOrderB(i+1, j, k) + &
                               ArrayOrderB( i, 1, k ) + ArrayOrderB(i, j-1, k) + &
                               ArrayOrderB( i, j, 1 ) + ArrayOrderB(i, j, k-1)

                     else

                        ordA = ArrayOrderA( ubound(f,3), j, k ) + ArrayOrderA(i+1, j, k) + &
                               ArrayOrderA( i, 1, k ) + ArrayOrderA(i, j-1, k) + &
                               ArrayOrderA( i, j, k+1 ) + ArrayOrderA(i, j, k-1)

                        ordB = ArrayOrderB( ubound(f,3), j, k ) + ArrayOrderB(i+1, j, k) + &
                               ArrayOrderB( i, 1, k ) + ArrayOrderB(i, j-1, k) + &
                               ArrayOrderB( i, j, k+1 ) + ArrayOrderB(i, j, k-1)

                     endif

                  else

                     if( k == 1 ) then

                        ordA = ArrayOrderA( ubound(f,3), j, k ) + ArrayOrderA(i+1, j, k) + &
                               ArrayOrderA( i, j-1, k ) + ArrayOrderA(i, j+1, k) + &
                               ArrayOrderA( i, j, ubound(f,5) ) + ArrayOrderA(i, j, k+1)

                        ordB = ArrayOrderB( ubound(f,3), j, k ) + ArrayOrderB(i+1, j, k) + &
                               ArrayOrderB( i, j-1, k ) + ArrayOrderB(i, j+1, k) + &
                               ArrayOrderB( i, j, ubound(f,5) ) + ArrayOrderB(i, j, k+1)

                     elseif( k == ubound(f,5) ) then
                            
                        ordA = ArrayOrderA( ubound(f,3), j, k ) + ArrayOrderA(i+1, j, k) + &
                               ArrayOrderA( i, j-1, k ) + ArrayOrderA(i, j+1, k) + &
                               ArrayOrderA( i, j, 1 ) + ArrayOrderA(i, j, k-1)

                        ordB = ArrayOrderB( ubound(f,3), j, k ) + ArrayOrderB(i+1, j, k) + &
                               ArrayOrderB( i, j-1, k ) + ArrayOrderB(i, j+1, k) + &
                               ArrayOrderB( i, j, 1 ) + ArrayOrderB(i, j, k-1)

                     else

                        ordA = ArrayOrderA( ubound(f,3), j, k ) + ArrayOrderA(i+1, j, k) + &
                               ArrayOrderA( i, j-1, k ) + ArrayOrderA(i, j+1, k) + &
                               ArrayOrderA( i, j, k+1 ) + ArrayOrderA(i, j, k-1)

                        ordB = ArrayOrderB( ubound(f,3), j, k ) + ArrayOrderB(i+1, j, k) + &
                               ArrayOrderB( i, j-1, k ) + ArrayOrderB(i, j+1, k) + &
                               ArrayOrderB( i, j, k+1 ) + ArrayOrderB(i, j, k-1)

                     endif

                  endif

               elseif( i == ubound(f,3) ) then

                  if( j == 1 ) then

                     if( k == 1 ) then

                        ordA = ArrayOrderA( 1, j, k ) + ArrayOrderA(i-1, j, k) + &
                               ArrayOrderA( i, ubound(f,4), k ) + ArrayOrderA(i, j+1, k) + &
                               ArrayOrderA( i, j, ubound(f,5) ) + ArrayOrderA(i, j, k+1)

                        ordB = ArrayOrderB( 1, j, k ) + ArrayOrderB(i-1, j, k) + &
                               ArrayOrderB( i, ubound(f,4), k ) + ArrayOrderB(i, j+1, k) + &
                               ArrayOrderB( i, j, ubound(f,5) ) + ArrayOrderB(i, j, k+1)

                     elseif( k == ubound(f,5) ) then
                            
                        ordA = ArrayOrderA( 1, j, k ) + ArrayOrderA(i-1, j, k) + &
                               ArrayOrderA( i, ubound(f,4), k ) + ArrayOrderA(i, j+1, k) + &
                               ArrayOrderA( i, j, 1 ) + ArrayOrderA(i, j, k-1)

                        ordB = ArrayOrderB( 1, j, k ) + ArrayOrderB(i-1, j, k) + &
                               ArrayOrderB( i, ubound(f,4), k ) + ArrayOrderB(i, j+1, k) + &
                               ArrayOrderB( i, j, 1 ) + ArrayOrderB(i, j, k-1)

                     else

                        ordA = ArrayOrderA( 1, j, k ) + ArrayOrderA(i-1, j, k) + &
                               ArrayOrderA( i, ubound(f,4), k ) + ArrayOrderA(i, j+1, k) + &
                               ArrayOrderA( i, j, k+1 ) + ArrayOrderA(i, j, k-1)

                        ordB = ArrayOrderB( 1, j, k ) + ArrayOrderB(i-1, j, k) + &
                               ArrayOrderB( i, ubound(f,4), k ) + ArrayOrderB(i, j+1, k) + &
                               ArrayOrderB( i, j, k+1 ) + ArrayOrderB(i, j, k-1)

                     endif

                  elseif( j == ubound(f,4) ) then

                     if( k == 1 ) then

                        ordA = ArrayOrderA( 1, j, k ) + ArrayOrderA(i-1, j, k) + &
                               ArrayOrderA( i, 1, k ) + ArrayOrderA(i, j-1, k) + &
                               ArrayOrderA( i, j, ubound(f,5) ) + ArrayOrderA(i, j, k+1)

                        ordB = ArrayOrderB( 1, j, k ) + ArrayOrderB(i-1, j, k) + &
                               ArrayOrderB( i, 1, k ) + ArrayOrderB(i, j-1, k) + &
                               ArrayOrderB( i, j, ubound(f,5) ) + ArrayOrderB(i, j, k+1)

                     elseif( k == ubound(f,5) ) then
                            
                        ordA = ArrayOrderA( 1, j, k ) + ArrayOrderA(i-1, j, k) + &
                               ArrayOrderA( i, 1, k ) + ArrayOrderA(i, j-1, k) + &
                               ArrayOrderA( i, j, 1 ) + ArrayOrderA(i, j, k-1)

                        ordB = ArrayOrderB( 1, j, k ) + ArrayOrderB(i-1, j, k) + &
                               ArrayOrderB( i, 1, k ) + ArrayOrderB(i, j-1, k) + &
                               ArrayOrderB( i, j, 1 ) + ArrayOrderB(i, j, k-1)

                     else

                        ordA = ArrayOrderA( 1, j, k ) + ArrayOrderA(i-1, j, k) + &
                               ArrayOrderA( i, 1, k ) + ArrayOrderA(i, j-1, k) + &
                               ArrayOrderA( i, j, k+1 ) + ArrayOrderA(i, j, k-1)

                        ordB = ArrayOrderB( 1, j, k ) + ArrayOrderB(i-1, j, k) + &
                               ArrayOrderB( i, 1, k ) + ArrayOrderB(i, j-1, k) + &
                               ArrayOrderB( i, j, k+1 ) + ArrayOrderB(i, j, k-1)

                     endif

                  else

                     if( k == 1 ) then

                        ordA = ArrayOrderA( 1, j, k ) + ArrayOrderA(i-1, j, k) + &
                               ArrayOrderA( i, j-1, k ) + ArrayOrderA(i, j+1, k) + &
                               ArrayOrderA( i, j, ubound(f,5) ) + ArrayOrderA(i, j, k+1)

                        ordB = ArrayOrderB( 1, j, k ) + ArrayOrderB(i-1, j, k) + &
                               ArrayOrderB( i, j-1, k ) + ArrayOrderB(i, j+1, k) + &
                               ArrayOrderB( i, j, ubound(f,5) ) + ArrayOrderB(i, j, k+1)

                     elseif( k == ubound(f,5) ) then
                            
                        ordA = ArrayOrderA( 1, j, k ) + ArrayOrderA(i-1, j, k) + &
                               ArrayOrderA( i, j-1, k ) + ArrayOrderA(i, j+1, k) + &
                               ArrayOrderA( i, j, 1 ) + ArrayOrderA(i, j, k-1)

                        ordB = ArrayOrderB( 1, j, k ) + ArrayOrderB(i-1, j, k) + &
                               ArrayOrderB( i, j-1, k ) + ArrayOrderB(i, j+1, k) + &
                               ArrayOrderB( i, j, 1 ) + ArrayOrderB(i, j, k-1)

                     else

                        ordA = ArrayOrderA( 1, j, k ) + ArrayOrderA(i-1, j, k) + &
                               ArrayOrderA( i, j-1, k ) + ArrayOrderA(i, j+1, k) + &
                               ArrayOrderA( i, j, k+1 ) + ArrayOrderA(i, j, k-1)

                        ordB = ArrayOrderB( 1, j, k ) + ArrayOrderB(i-1, j, k) + &
                               ArrayOrderB( i, j-1, k ) + ArrayOrderB(i, j+1, k) + &
                               ArrayOrderB( i, j, k+1 ) + ArrayOrderB(i, j, k-1)

                     endif

                  endif

               else

                  if( j == 1 ) then

                     if( k == 1 ) then

                        ordA = ArrayOrderA( i-1, j, k ) + ArrayOrderA(i+1, j, k) + &
                               ArrayOrderA( i, ubound(f,4), k ) + ArrayOrderA(i, j+1, k) + &
                               ArrayOrderA( i, j, ubound(f,5) ) + ArrayOrderA(i, j, k+1)

                        ordB = ArrayOrderB( i-1, j, k ) + ArrayOrderB(i+1, j, k) + &
                               ArrayOrderB( i, ubound(f,4), k ) + ArrayOrderB(i, j+1, k) + &
                               ArrayOrderB( i, j, ubound(f,5) ) + ArrayOrderB(i, j, k+1)

                     elseif( k == ubound(f,5) ) then
                            
                        ordA = ArrayOrderA( i-1, j, k ) + ArrayOrderA(i+1, j, k) + &
                               ArrayOrderA( i, ubound(f,4), k ) + ArrayOrderA(i, j+1, k) + &
                               ArrayOrderA( i, j, 1 ) + ArrayOrderA(i, j, k-1)

                        ordB = ArrayOrderB( i-1, j, k ) + ArrayOrderB(i+1, j, k) + &
                               ArrayOrderB( i, ubound(f,4), k ) + ArrayOrderB(i, j+1, k) + &
                               ArrayOrderB( i, j, 1 ) + ArrayOrderB(i, j, k-1)

                     else

                        ordA = ArrayOrderA( i-1, j, k ) + ArrayOrderA(i+1, j, k) + &
                               ArrayOrderA( i, ubound(f,4), k ) + ArrayOrderA(i, j+1, k) + &
                               ArrayOrderA( i, j, k+1 ) + ArrayOrderA(i, j, k-1)

                        ordB = ArrayOrderB( i-1, j, k ) + ArrayOrderB(i+1, j, k) + &
                               ArrayOrderB( i, ubound(f,4), k ) + ArrayOrderB(i, j+1, k) + &
                               ArrayOrderB( i, j, k+1 ) + ArrayOrderB(i, j, k-1)

                     endif

                  elseif( j == ubound(f,4) ) then

                     if( k == 1 ) then

                        ordA = ArrayOrderA( i-1, j, k ) + ArrayOrderA(i+1, j, k) + &
                               ArrayOrderA( i, 1, k ) + ArrayOrderA(i, j-1, k) + &
                               ArrayOrderA( i, j, ubound(f,5) ) + ArrayOrderA(i, j, k+1)

                        ordB = ArrayOrderB( i-1, j, k ) + ArrayOrderB(i+1, j, k) + &
                               ArrayOrderB( i, 1, k ) + ArrayOrderB(i, j-1, k) + &
                               ArrayOrderB( i, j, ubound(f,5) ) + ArrayOrderB(i, j, k+1)

                     elseif( k == ubound(f,5) ) then
                            
                        ordA = ArrayOrderA( i-1, j, k ) + ArrayOrderA(i+1, j, k) + &
                               ArrayOrderA( i, 1, k ) + ArrayOrderA(i, j-1, k) + &
                               ArrayOrderA( i, j, 1 ) + ArrayOrderA(i, j, k-1)

                        ordB = ArrayOrderB( i-1, j, k ) + ArrayOrderB(i+1, j, k) + &
                               ArrayOrderB( i, 1, k ) + ArrayOrderB(i, j-1, k) + &
                               ArrayOrderB( i, j, 1 ) + ArrayOrderB(i, j, k-1)

                     else

                        ordA = ArrayOrderA( i-1, j, k ) + ArrayOrderA(i+1, j, k) + &
                               ArrayOrderA( i, 1, k ) + ArrayOrderA(i, j-1, k) + &
                               ArrayOrderA( i, j, k+1 ) + ArrayOrderA(i, j, k-1)

                        ordB = ArrayOrderB( i-1, j, k ) + ArrayOrderB(i+1, j, k) + &
                               ArrayOrderB( i, 1, k ) + ArrayOrderB(i, j-1, k) + &
                               ArrayOrderB( i, j, k+1 ) + ArrayOrderB(i, j, k-1)

                     endif

                  else

                     if( k == 1 ) then

                        ordA = ArrayOrderA( i-1, j, k ) + ArrayOrderA(i+1, j, k) + &
                               ArrayOrderA( i, j-1, k ) + ArrayOrderA(i, j+1, k) + &
                               ArrayOrderA( i, j, ubound(f,5) ) + ArrayOrderA(i, j, k+1)

                        ordB = ArrayOrderB( i-1, j, k ) + ArrayOrderB(i+1, j, k) + &
                               ArrayOrderB( i, j-1, k ) + ArrayOrderB(i, j+1, k) + &
                               ArrayOrderB( i, j, ubound(f,5) ) + ArrayOrderB(i, j, k+1)

                     elseif( k == ubound(f,5) ) then
                            
                        ordA = ArrayOrderA( i-1, j, k ) + ArrayOrderA(i+1, j, k) + &
                               ArrayOrderA( i, j-1, k ) + ArrayOrderA(i, j+1, k) + &
                               ArrayOrderA( i, j, 1 ) + ArrayOrderA(i, j, k-1)

                        ordB = ArrayOrderB( i-1, j, k ) + ArrayOrderB(i+1, j, k) + &
                               ArrayOrderB( i, j-1, k ) + ArrayOrderB(i, j+1, k) + &
                               ArrayOrderB( i, j, 1 ) + ArrayOrderB(i, j, k-1)

                     else

                        ordA = ArrayOrderA( i-1, j, k ) + ArrayOrderA(i+1, j, k) + &
                               ArrayOrderA( i, j-1, k ) + ArrayOrderA(i, j+1, k) + &
                               ArrayOrderA( i, j, k+1 ) + ArrayOrderA(i, j, k-1)

                        ordB = ArrayOrderB( i-1, j, k ) + ArrayOrderB(i+1, j, k) + &
                               ArrayOrderB( i, j-1, k ) + ArrayOrderB(i, j+1, k) + &
                               ArrayOrderB( i, j, k+1 ) + ArrayOrderB(i, j, k-1)

                     endif

                  endif

               endif

                       
               do na = 0, ubound(f,1)
                  do nb = 0, ubound(f,2)

                     ! chop
                     if( abs( real( f(na,nb,i,j,k) ) ) < chopCutoff ) then
                        f(na,nb,i,j,k) = ( f(na,nb,i,j,k) - conjg(f(na,nb,i,j,k)) )/2.0_dp
                     endif
                     if( abs( aimag( f(na,nb,i,j,k) ) ) < chopCutoff ) then
                        f(na,nb,i,j,k) = ( f(na,nb,i,j,k) + conjg(f(na,nb,i,j,k)) )/2.0_dp
                     endif

                     kv(na,nb,i,j,k) = mlt * ( real(na,dp) * ( mua - ua/2.0_dp * (real(na,dp) - 1.0_dp) ) + &
                                               real(nb,dp) * ( mub - ub/2.0_dp * (real(nb,dp) - 1.0_dp) ) - &
                                               uab * real(na,dp) * real(nb,dp) ) * f(na,nb,i,j,k)

                     if( na < ubound(f,1) ) then
                        kv(na,nb,i,j,k) = kv(na,nb,i,j,k) + mlt * ja * ordA * SQRT( real(na+1,dp) ) * f(na+1,nb,i,j,k)
                     endif
                     if( na > 0 ) then
                        kv(na,nb,i,j,k) = kv(na,nb,i,j,k) + mlt * ja * ordA * SQRT( real(na,dp) ) * f(na-1,nb,i,j,k)
                     endif

                     if( nb < ubound(f,2) ) then
                        kv(na,nb,i,j,k) = kv(na,nb,i,j,k) + mlt * jb * ordB * SQRT( real(nb+1,dp) ) * f(na,nb+1,i,j,k)
                     endif
                     if( nb > 0 ) then
                        kv(na,nb,i,j,k) = kv(na,nb,i,j,k) + mlt * jb * ordB * SQRT( real(nb,dp) ) * f(na,nb-1,i,j,k)
                     endif

                  enddo
               enddo

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
         call c1d_solve( f, k1, ja, jb, ua, ub, uab, mua, mub, 'GS' )
         ! update vector to the 2nd step
         tmp(:,:,:) = f(:,:,:) + dt/2.0_dp * k1(:,:,:)
         ! 2nd Runge-Kutta step
         call c1d_solve( tmp, k2, ja, jb, ua, ub, uab, mua, mub, 'GS' )
         ! update vector for the 3rd step
         tmp(:,:,:) = f(:,:,:) + dt/2.0_dp * k2(:,:,:)
         ! 3rd Runge-Kutta step
         call c1d_solve( tmp, k3, ja, jb, ua, ub, uab, mua, mub, 'GS' )
         ! update vector for the 4th step
         tmp(:,:,:) = f(:,:,:) + dt * k3(:,:,:)
         ! 4th Runge-Kutta step
         call c1d_solve( tmp, k4, ja, jb, ua, ua, uab, mua, mub, 'GS' )
         ! update vector
         f(:,:,:) = f(:,:,:) + dt/6.0_dp * ( k1(:,:,:) + 2.0_dp *&
                    k2(:,:,:) + 2.0_dp * k3(:,:,:) + k4(:,:,:) )

         ! normalize vector
         call const_mean( f, mna, mnb, info )
         
         if( info /= 0 ) then
            ! suppress stdout if quietON
            if( .not. quietON ) print *, 'const_mean info = ', info
            exit
         endif

         ! calculate error
         if( modulo(cnt, stepsForJudge) .eq. 0 ) then

            print *, cnt
            error = 0.0_dp
            do i = 1, ubound(f,3)
               do na = 0, ubound(f,1)
                  do nb = 0, ubound(f,2)

                     ! absolute error
                     ierr = abs( ( f(nA,nB,i) - fpom(nA,nB,i) ) )
                     if (ierr .gt. error) then
                        error = ierr
                     endif
                  enddo
               enddo
            enddo

            if( error < convCriterion ) then
               exit
            endif

            if( .not. quietON ) print *, 'GroundStateNC: error = ', error
            fpom(:,:,:) = f(:,:,:)

         endif

      end do

   end subroutine

   subroutine c2d_GroundStateNC( f, ja, jb, ua, ub, uab, mna, mnb )

      implicit none
      complex( kind = dp ), intent(inout), allocatable :: f(:,:,:,:)
      real( kind = dp ), intent(in) :: ja, jb, ua, ub, uab, mna, mnb

      ! | -------------------- |
      ! | For Runge-Kutta step |
      ! | -------------------- |
      complex( kind = dp ), allocatable :: k1(:,:,:,:)
      complex( kind = dp ), allocatable :: k2(:,:,:,:)
      complex( kind = dp ), allocatable :: k3(:,:,:,:)
      complex( kind = dp ), allocatable :: k4(:,:,:,:)
      complex( kind = dp ), allocatable :: tmp(:,:,:,:)
      complex( kind = dp ), allocatable :: fpom(:,:,:,:)
      
      real( kind = dp ) :: error, ierr, mua, mub
      integer i, j, na, nb, cnt, info

      allocate( k1(lbound(f,1):ubound(f,1), lbound(f,2):ubound(f,2), lbound(f,3):ubound(f,3), lbound(f,4):ubound(f,4)) )
      allocate( k2(lbound(f,1):ubound(f,1), lbound(f,2):ubound(f,2), lbound(f,3):ubound(f,3), lbound(f,4):ubound(f,4)) )
      allocate( k3(lbound(f,1):ubound(f,1), lbound(f,2):ubound(f,2), lbound(f,3):ubound(f,3), lbound(f,4):ubound(f,4)) )
      allocate( k4(lbound(f,1):ubound(f,1), lbound(f,2):ubound(f,2), lbound(f,3):ubound(f,3), lbound(f,4):ubound(f,4)) )
      allocate( tmp(lbound(f,1):ubound(f,1), lbound(f,2):ubound(f,2), lbound(f,3):ubound(f,3), lbound(f,4):ubound(f,4)) )
      allocate( fpom(lbound(f,1):ubound(f,1), lbound(f,2):ubound(f,2), lbound(f,3):ubound(f,3), lbound(f,4):ubound(f,4)) )

      fpom(:,:,:,:) = f(:,:,:,:)
      cnt = 0
      mua = 0.0_dp; mub = 0.0_dp

      do while( .true. )

         cnt = cnt + 1 
 
         ! 1st Runge-Kutta step
         call c2d_solve( f, k1, ja, jb, ua, ub, uab, mua, mub, 'GS' )
         ! update vector to the 2nd step
         tmp(:,:,:,:) = f(:,:,:,:) + dt/2.0_dp * k1(:,:,:,:)
         ! 2nd Runge-Kutta step
         call c2d_solve( tmp, k2, ja, jb, ua, ub, uab, mua, mub, 'GS' )
         ! update vector for the 3rd step
         tmp(:,:,:,:) = f(:,:,:,:) + dt/2.0_dp * k2(:,:,:,:)
         ! 3rd Runge-Kutta step
         call c2d_solve( tmp, k3, ja, jb, ua, ub, uab, mua, mub, 'GS' )
         ! update vector for the 4th step
         tmp(:,:,:,:) = f(:,:,:,:) + dt * k3(:,:,:,:)
         ! 4th Runge-Kutta step
         call c2d_solve( tmp, k4, ja, jb, ua, ua, uab, mua, mub, 'GS' )
         ! update vector
         f(:,:,:,:) = f(:,:,:,:) + dt/6.0_dp * ( k1(:,:,:,:) + 2.0_dp *&
                      k2(:,:,:,:) + 2.0_dp * k3(:,:,:,:) + k4(:,:,:,:) )

         ! normalize vector
         call const_mean( f, mna, mnb, info )
         
         if( info /= 0 ) then
            ! suppress stdout if quietON
            if( .not. quietON ) print *, 'const_mean info = ', info
            exit
         endif

         ! calculate error
         if( modulo(cnt, stepsForJudge) .eq. 0 ) then

            print *, cnt
            error = 0.0_dp
            do i = 1, ubound(f,3)
               do j = 1, ubound(f,4)
                  do na = 0, ubound(f,1)
                     do nb = 0, ubound(f,2)
                        ierr = abs( ( f(nA,nB,i,j) - fpom(nA,nB,i,j) ) )
                        if (ierr .gt. error) then
                           error = ierr
                        endif
                     enddo
                  enddo
               enddo
            enddo

            if( error < convCriterion ) then
               exit
            endif

            if( .not. quietON ) print *, 'GroundStateNC: error = ', error
            fpom(:,:,:,:) = f(:,:,:,:)

         endif

      end do

   end subroutine

   subroutine c3d_GroundStateNC( f, ja, jb, ua, ub, uab, mna, mnb )

      implicit none
      complex( kind = dp ), intent(inout), allocatable :: f(:,:,:,:,:)
      real( kind = dp ), intent(in) :: ja, jb, ua, ub, uab, mna, mnb

      ! | -------------------- |
      ! | For Runge-Kutta step |
      ! | -------------------- |
      complex( kind = dp ), allocatable :: k1(:,:,:,:,:)
      complex( kind = dp ), allocatable :: k2(:,:,:,:,:)
      complex( kind = dp ), allocatable :: k3(:,:,:,:,:)
      complex( kind = dp ), allocatable :: k4(:,:,:,:,:)
      complex( kind = dp ), allocatable :: tmp(:,:,:,:,:)
      complex( kind = dp ), allocatable :: fpom(:,:,:,:,:)
      
      real( kind = dp ) :: error, ierr, mua, mub
      integer i, j, k, na, nb, cnt, info

      allocate( k1(lbound(f,1):ubound(f,1), lbound(f,2):ubound(f,2), lbound(f,3):ubound(f,3), &
                lbound(f,4):ubound(f,4), lbound(f,5):ubound(f,5)) )
      allocate( k2(lbound(f,1):ubound(f,1), lbound(f,2):ubound(f,2), lbound(f,3):ubound(f,3), &
                lbound(f,4):ubound(f,4), lbound(f,5):ubound(f,5)) )
      allocate( k3(lbound(f,1):ubound(f,1), lbound(f,2):ubound(f,2), lbound(f,3):ubound(f,3), &
                lbound(f,4):ubound(f,4), lbound(f,5):ubound(f,5)) )
      allocate( k4(lbound(f,1):ubound(f,1), lbound(f,2):ubound(f,2), lbound(f,3):ubound(f,3), &
                lbound(f,4):ubound(f,4), lbound(f,5):ubound(f,5)) )
      allocate( tmp(lbound(f,1):ubound(f,1), lbound(f,2):ubound(f,2), lbound(f,3):ubound(f,3), &
                lbound(f,4):ubound(f,4), lbound(f,5):ubound(f,5)) )
      allocate( fpom(lbound(f,1):ubound(f,1), lbound(f,2):ubound(f,2), lbound(f,3):ubound(f,3), &
                lbound(f,4):ubound(f,4), lbound(f,5):ubound(f,5)) )

      fpom(:,:,:,:,:) = f(:,:,:,:,:)
      cnt = 0
      mua = 0.0_dp; mub = 0.0_dp

      do while( .true. )

         cnt = cnt + 1 
 
         ! 1st Runge-Kutta step
         call c3d_solve( f, k1, ja, jb, ua, ub, uab, mua, mub, 'GS' )
         ! update vector to the 2nd step
         tmp(:,:,:,:,:) = f(:,:,:,:,:) + dt/2.0_dp * k1(:,:,:,:,:)
         ! 2nd Runge-Kutta step
         call c3d_solve( tmp, k2, ja, jb, ua, ub, uab, mua, mub, 'GS' )
         ! update vector for the 3rd step
         tmp(:,:,:,:,:) = f(:,:,:,:,:) + dt/2.0_dp * k2(:,:,:,:,:)
         ! 3rd Runge-Kutta step
         call c3d_solve( tmp, k3, ja, jb, ua, ub, uab, mua, mub, 'GS' )
         ! update vector for the 4th step
         tmp(:,:,:,:,:) = f(:,:,:,:,:) + dt * k3(:,:,:,:,:)
         ! 4th Runge-Kutta step
         call c3d_solve( tmp, k4, ja, jb, ua, ua, uab, mua, mub, 'GS' )
         ! update vector
         f(:,:,:,:,:) = f(:,:,:,:,:) + dt/6.0_dp * ( k1(:,:,:,:,:) + 2.0_dp *&
                        k2(:,:,:,:,:) + 2.0_dp * k3(:,:,:,:,:) + k4(:,:,:,:,:) )

         ! normalize vector
         call const_mean( f, mna, mnb, info )
         
         if( info /= 0 ) then
            ! suppress stdout if quietON
            if( .not. quietON ) print *, 'const_mean info = ', info
            exit
         endif


         ! calculate error
         if( modulo(cnt, stepsForJudge) .eq. 0 ) then

            print *, cnt
            error = 0.0_dp
            do i = 1, ubound(f,3)
               do j = 1, ubound(f,4)
                  do k = 1, ubound(f,5)
                     do na = 0, ubound(f,1)
                        do nb = 0, ubound(f,2)
                           ierr = abs( ( f(nA,nB,i,j,k) - fpom(nA,nB,i,j,k) ) )
                           if (ierr .gt. error) then
                              error = ierr
                           endif
                        enddo
                     enddo
                  enddo
               enddo
            enddo

            if( error < convCriterion ) then
               exit
            endif

            if( .not. quietON ) print *, 'GroundStateNC: error = ', error
            fpom(:,:,:,:,:) = f(:,:,:,:,:)

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
         call c1d_solve( f, k1, ja, jb, ua, ub, uab, mua, mub, 'GS' )
         ! update vector to the 2nd step
         tmp(:,:,:) = f(:,:,:) + dt/2.0_dp * k1(:,:,:)
         ! 2nd Runge-Kutta step
         call c1d_solve( tmp, k2, ja, jb, ua, ub, uab, mua, mub, 'GS' )
         ! update vector for the 3rd step
         tmp(:,:,:) = f(:,:,:) + dt/2.0_dp * k2(:,:,:)
         ! 3rd Runge-Kutta step
         call c1d_solve( tmp, k3, ja, jb, ua, ub, uab, mua, mub, 'GS' )
         ! update vector for the 4th step
         tmp(:,:,:) = f(:,:,:) + dt * k3(:,:,:)
         ! 4th Runge-Kutta step
         call c1d_solve( tmp, k4, ja, jb, ua, ua, uab, mua, mub, 'GS' )
         ! update vector
         f(:,:,:) = f(:,:,:) + dt/6.0_dp * ( k1(:,:,:) + 2.0_dp *&
                    k2(:,:,:) + 2.0_dp * k3(:,:,:) + k4(:,:,:) )

         ! normalize vector
         call normalize( f )

         ! calculate error
         if( modulo(cnt, stepsForJudge) .eq. 0 ) then

            print *, cnt
            error = 0.0_dp
            do i = 1, ubound(f,3)
               do na = 0, ubound(f,1)
                  do nb = 0, ubound(f,2)
                     ierr = abs( ( f(nA,nB,i) - fpom(nA,nB,i) ) )
                     if (ierr .gt. error) then
                        error = ierr
                     endif
                  enddo
               enddo
            enddo

            if( error .lt. convCriterion ) then
               exit
            endif

            if( .not. quietON ) print *, 'GroundStateNC: error = ', error
            fpom(:,:,:) = f(:,:,:)

         endif

      end do

   end subroutine

   subroutine c2d_GroundState( f, ja, jb, ua, ub, uab, mua, mub )

      implicit none
      complex( kind = dp ), intent(inout), allocatable :: f(:,:,:,:)
      real( kind = dp ), intent(in) :: ja, jb, ua, ub, uab, mua, mub

      ! | -------------------- |
      ! | For Runge-Kutta step |
      ! | -------------------- |
      complex( kind = dp ), allocatable :: k1(:,:,:,:)
      complex( kind = dp ), allocatable :: k2(:,:,:,:)
      complex( kind = dp ), allocatable :: k3(:,:,:,:)
      complex( kind = dp ), allocatable :: k4(:,:,:,:)
      complex( kind = dp ), allocatable :: tmp(:,:,:,:)
      complex( kind = dp ), allocatable :: fpom(:,:,:,:)
      
      real( kind = dp ) :: error, ierr
      integer i, j, na, nb, cnt

      allocate( k1(lbound(f,1):ubound(f,1), lbound(f,2):ubound(f,2), lbound(f,3):ubound(f,3), &
                lbound(f,4):ubound(f,4)) )
      allocate( k2(lbound(f,1):ubound(f,1), lbound(f,2):ubound(f,2), lbound(f,3):ubound(f,3), &
                lbound(f,4):ubound(f,4)) )
      allocate( k3(lbound(f,1):ubound(f,1), lbound(f,2):ubound(f,2), lbound(f,3):ubound(f,3), &
                lbound(f,4):ubound(f,4)) )
      allocate( k4(lbound(f,1):ubound(f,1), lbound(f,2):ubound(f,2), lbound(f,3):ubound(f,3), &
                lbound(f,4):ubound(f,4)) )
      allocate( tmp(lbound(f,1):ubound(f,1), lbound(f,2):ubound(f,2), lbound(f,3):ubound(f,3), &
                lbound(f,4):ubound(f,4)) )
      allocate( fpom(lbound(f,1):ubound(f,1), lbound(f,2):ubound(f,2), lbound(f,3):ubound(f,3), &
                lbound(f,4):ubound(f,4)) )

      fpom(:,:,:,:) = f(:,:,:,:)
      cnt = 0

      do while( .true. )

         cnt = cnt + 1 
 
         ! 1st Runge-Kutta step
         call c2d_solve( f, k1, ja, jb, ua, ub, uab, mua, mub, 'GS' )
         ! update vector to the 2nd step
         tmp(:,:,:,:) = f(:,:,:,:) + dt/2.0_dp * k1(:,:,:,:)
         ! 2nd Runge-Kutta step
         call c2d_solve( tmp, k2, ja, jb, ua, ub, uab, mua, mub, 'GS' )
         ! update vector for the 3rd step
         tmp(:,:,:,:) = f(:,:,:,:) + dt/2.0_dp * k2(:,:,:,:)
         ! 3rd Runge-Kutta step
         call c2d_solve( tmp, k3, ja, jb, ua, ub, uab, mua, mub, 'GS' )
         ! update vector for the 4th step
         tmp(:,:,:,:) = f(:,:,:,:) + dt * k3(:,:,:,:)
         ! 4th Runge-Kutta step
         call c2d_solve( tmp, k4, ja, jb, ua, ua, uab, mua, mub, 'GS' )
         ! update vector
         f(:,:,:,:) = f(:,:,:,:) + dt/6.0_dp * ( k1(:,:,:,:) + 2.0_dp *&
                      k2(:,:,:,:) + 2.0_dp * k3(:,:,:,:) + k4(:,:,:,:) )

         ! normalize vector
         call normalize( f )
         
         ! calculate error
         if( modulo(cnt, stepsForJudge) .eq. 0 ) then

            print *, cnt
            error = 0.0_dp
            do i = 1, ubound(f,3)
               do j = 1, ubound(f,4)
                  do na = 0, ubound(f,1)
                     do nb = 0, ubound(f,2)
                        ierr = abs( ( f(nA,nB,i,j) - fpom(nA,nB,i,j) ) )
                        if (ierr .gt. error) then
                           error = ierr
                        endif
                     enddo
                  enddo
               enddo
            enddo

            if( error < convCriterion ) then
               exit
            endif

            if( .not. quietON ) print *, 'GroundStateNC: error = ', error
            fpom(:,:,:,:) = f(:,:,:,:)

         endif

      end do

   end subroutine

   subroutine c3d_GroundState( f, ja, jb, ua, ub, uab, mua, mub )

      implicit none
      complex( kind = dp ), intent(inout), allocatable :: f(:,:,:,:,:)
      real( kind = dp ), intent(in) :: ja, jb, ua, ub, uab, mua, mub

      ! | -------------------- |
      ! | For Runge-Kutta step |
      ! | -------------------- |
      complex( kind = dp ), allocatable :: k1(:,:,:,:,:)
      complex( kind = dp ), allocatable :: k2(:,:,:,:,:)
      complex( kind = dp ), allocatable :: k3(:,:,:,:,:)
      complex( kind = dp ), allocatable :: k4(:,:,:,:,:)
      complex( kind = dp ), allocatable :: tmp(:,:,:,:,:)
      complex( kind = dp ), allocatable :: fpom(:,:,:,:,:)
      
      real( kind = dp ) :: error, ierr
      integer i, j, k, na, nb, cnt

      allocate( k1(lbound(f,1):ubound(f,1), lbound(f,2):ubound(f,2), lbound(f,3):ubound(f,3), &
                lbound(f,4):ubound(f,4), lbound(f,5):ubound(f,5)) )
      allocate( k2(lbound(f,1):ubound(f,1), lbound(f,2):ubound(f,2), lbound(f,3):ubound(f,3), &
                lbound(f,4):ubound(f,4), lbound(f,5):ubound(f,5)) )
      allocate( k3(lbound(f,1):ubound(f,1), lbound(f,2):ubound(f,2), lbound(f,3):ubound(f,3), &
                lbound(f,4):ubound(f,4), lbound(f,5):ubound(f,5)) )
      allocate( k4(lbound(f,1):ubound(f,1), lbound(f,2):ubound(f,2), lbound(f,3):ubound(f,3), &
                lbound(f,4):ubound(f,4), lbound(f,5):ubound(f,5)) )
      allocate( tmp(lbound(f,1):ubound(f,1), lbound(f,2):ubound(f,2), lbound(f,3):ubound(f,3), &
                lbound(f,4):ubound(f,4), lbound(f,5):ubound(f,5)) )
      allocate( fpom(lbound(f,1):ubound(f,1), lbound(f,2):ubound(f,2), lbound(f,3):ubound(f,3), &
                lbound(f,4):ubound(f,4), lbound(f,5):ubound(f,5)) )

      fpom(:,:,:,:,:) = f(:,:,:,:,:)
      cnt = 0

      do while( .true. )

         cnt = cnt + 1 
 
         ! 1st Runge-Kutta step
         call c3d_solve( f, k1, ja, jb, ua, ub, uab, mua, mub, 'GS' )
         ! update vector to the 2nd step
         tmp(:,:,:,:,:) = f(:,:,:,:,:) + dt/2.0_dp * k1(:,:,:,:,:)
         ! 2nd Runge-Kutta step
         call c3d_solve( tmp, k2, ja, jb, ua, ub, uab, mua, mub, 'GS' )
         ! update vector for the 3rd step
         tmp(:,:,:,:,:) = f(:,:,:,:,:) + dt/2.0_dp * k2(:,:,:,:,:)
         ! 3rd Runge-Kutta step
         call c3d_solve( tmp, k3, ja, jb, ua, ub, uab, mua, mub, 'GS' )
         ! update vector for the 4th step
         tmp(:,:,:,:,:) = f(:,:,:,:,:) + dt * k3(:,:,:,:,:)
         ! 4th Runge-Kutta step
         call c3d_solve( tmp, k4, ja, jb, ua, ua, uab, mua, mub, 'GS' )
         ! update vector
         f(:,:,:,:,:) = f(:,:,:,:,:) + dt/6.0_dp * ( k1(:,:,:,:,:) + 2.0_dp *&
                        k2(:,:,:,:,:) + 2.0_dp * k3(:,:,:,:,:) + k4(:,:,:,:,:) )

         ! normalize vector
         call normalize( f )
         
         ! calculate error
         if( modulo(cnt, stepsForJudge) .eq. 0 ) then

            print *, cnt
            error = 0.0_dp
            do i = 1, ubound(f,3)
               do j = 1, ubound(f,4)
                  do k = 1, ubound(f,5)
                     do na = 0, ubound(f,1)
                        do nb = 0, ubound(f,2)
                           ierr = abs( ( f(nA,nB,i,j,k) - fpom(nA,nB,i,j,k) ) )
                           if (ierr .gt. error) then
                              error = ierr
                           endif
                        enddo
                     enddo
                  enddo
               enddo
            enddo

            if( error < convCriterion ) then
               exit
            endif

            if( .not. quietON ) print *, 'GroundStateNC: error = ', error
            fpom(:,:,:,:,:) = f(:,:,:,:,:)

         endif

      end do

   end subroutine

!  | ---------------------------------------------- |
!  | Initialize uniformly with unit normalization   | 
!  | ---------------------------------------------- |
   subroutine c1d_InitUniform( f )
 
      implicit none
      complex( kind = dp ), intent(inout), allocatable :: f(:,:,:)

      f(:,:,:) = (1.0_dp, 1.0_dp)
      call normalize( f )

   end subroutine

   subroutine c2d_InitUniform( f )
 
      implicit none
      complex( kind = dp ), intent(inout), allocatable :: f(:,:,:,:)

      f(:,:,:,:) = (1.0_dp, 1.0_dp)
      call normalize( f )

   end subroutine

   subroutine c3d_InitUniform( f )
 
      implicit none
      complex( kind = dp ), intent(inout), allocatable :: f(:,:,:,:,:)

      f(:,:,:,:,:) = (1.0_dp, 1.0_dp)
      call normalize( f )

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

   subroutine c2d_InitUniformNC( f, mna, mnb )

      implicit none
      complex( kind = dp ), intent(inout), allocatable :: f(:,:,:,:)
      real( kind = dp ), intent(in) :: mna, mnb
      integer info

      f(:,:,:,:) = (1.0_dp, 1.0_dp)
      call const_mean( f, mna, mnb, info )

   end subroutine

   subroutine c3d_InitUniformNC( f, mna, mnb )

      implicit none
      complex( kind = dp ), intent(inout), allocatable :: f(:,:,:,:,:)
      real( kind = dp ), intent(in) :: mna, mnb
      integer info

      f(:,:,:,:,:) = (1.0_dp, 1.0_dp)
      call const_mean( f, mna, mnb, info )

   end subroutine

end module
