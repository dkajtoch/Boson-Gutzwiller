module solvers

   use parameters
   use measures, only: normalize
   use fixed_mean, only: const_mean

   private
   public GroundState, GroundStateNC, InitUniform, InitUniformNC, &
          InitRandom, InitRandomNC

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

   interface InitRandom
      module procedure c1d_InitRandom
      module procedure c2d_InitRandom
      module procedure c3d_InitRandom
   end interface

   interface InitUniformNC
      module procedure c1d_InitUniformNC
      module procedure c2d_InitUniformNC
      module procedure c3d_InitUniformNC
   end interface

   interface InitRandomNC
      module procedure c1d_InitRandomNC
      module procedure c2d_InitRandomNC
      module procedure c3d_InitRandomNC
   end interface

   contains

!  | ---------------------------------------------- |
!  | Update Gutzwiller coeffcients using right-hand | 
!  | side of Gutzwiller eqution                    |
!  | ---------------------------------------------- |
   subroutine c1d_solve( f, k, ArrayOrder, ja, u, mu, mode )

      implicit none
      complex( kind = dp), intent(inout), allocatable :: f(:,:)
      complex( kind = dp ), intent(inout), allocatable :: k(:,:)
      complex( kind = dp ), intent(inout), allocatable :: ArrayOrder(:)
      real( kind = dp ), intent(in) :: ja, u, mu
      character(len=2), intent(in) :: mode

      ! local variables
      integer i, na
      complex( kind = dp ) :: ord, mlt

      if (mode .eq. 'EV') then
         mlt = im
      else if (mode .eq. 'GS') then
         mlt = re
      else
         mlt = (0.0_dp, 0.0_dp) 
         write(*,*)  "Wrong mode in SOLVE subroutine!"
      endif

      ! collect order parameters
      do i = 1, ubound(f,2)
         
         ArrayOrder(i) = 0.0_dp*re

         do na = 0, ubound(f,1)
               
            if( na > 0 ) then
               ArrayOrder(i) = ArrayOrder(i) + SQRT( real(na,dp) ) * conjg( f(na-1,i) ) * f(na,i)
            endif

         enddo

      enddo

      ! solve coupled system of equtions
      do i = 1, ubound(f,2)

         ! calculate nearest-neighbour order parameter (periodic boundary conditions)

         if( i == 1 ) then
            ord = ArrayOrder( 2 ) + ArrayOrder( ubound(f,2) )
         elseif( i == ubound(f,2) ) then
            ord = ArrayOrder( 1 ) + ArrayOrder( ubound(f,2)-1 )
         else
            ord = ArrayOrder( i-1 ) + ArrayOrder( i+1 )
         endif

         do na = 0, ubound(f,1)

            ! chop
            if( abs( real( f(na,i) ) ) < chopCutoff ) then
                  f(na,i) = ( f(na,i) - conjg(f(na,i)) )/2.0_dp
            endif
            if( abs( aimag( f(na,i) ) ) < chopCutoff ) then
               f(na,i) = ( f(na,i) + conjg(f(na,i)) )/2.0_dp
            endif


            k(na,i) = mlt * ( real(na,dp) * ( mu - u/2.0_dp * (real(na,dp) - 1.0_dp) ) * f(na,i)

            if( na < ubound(f,1) ) then
               k(na,i) = k(na,i) + mlt * ja * ord * SQRT( real(na+1,dp) ) * f(na+1,i)
            endif
            if( na > 0 ) then
               k(na,i) = k(na,i) + mlt * ja * ord * SQRT( real(na,dp) ) * f(na-1,i)
            endif

         enddo ! na

      enddo ! i

   end subroutine

   subroutine c2d_solve( f, k, ArrayOrder, ja, u, mu, mode )

      implicit none
      complex( kind = dp), intent(inout), allocatable :: f(:,:,:)
      complex( kind = dp ), intent(inout), allocatable :: k(:,:,:)
      complex( kind = dp ), allocatable :: ArrayOrder(:,:)
      real( kind = dp ), intent(in) :: ja, u, mu
      character(len=2), intent(in) :: mode

      ! local variables
      integer i, j, na
      complex( kind = dp ) :: ord, mlt

      if (mode .eq. 'EV') then
         mlt = im
      else if (mode .eq. 'GS') then
         mlt = re
      else
         mlt = (0.0_dp, 0.0_dp) 
         write(*,*)  "Wrong mode in SOLVE subroutine!"
      endif

      ! collect order parameters
      do i = 1, ubound(f,2)
         do j = 1, ubound(f,3)
         
            ArrayOrder(i,j) = 0.0_dp*re

            do na = 0, ubound(f,1)
               
               if( na > 0 ) then
                  ArrayOrder(i,j) = ArrayOrder(i,j) + SQRT( real(na,dp) ) * conjg( f(na-1,i,j) ) * f(na,i,j)
               endif

            enddo

         enddo
      enddo

      ! solve coupled system of equtions
      do i = 1, ubound(f,2)
         do j = 1, ubound(f,3)

            ! calculate nearest-neighbour order parameter (periodic boundary conditions)
            ord = (0.0_dp, 0.0_dp)

            if( i == 1 ) then
               ord = ord + ArrayOrder( 2, j ) + ArrayOrder( ubound(f,2), j )
            elseif( i == ubound(f,2) ) then
               ord = ord + ArrayOrder( 1, j ) + ArrayOrder( ubound(f,2)-1, j )
            else
               ord = ord + ArrayOrder( i+1, j ) + ArrayOrder( i-1, j )
            endif

            if( j == 1 ) then
               ord = ord + ArrayOrder( i, 2 ) + ArrayOrder( i, ubound(f,3) )
            elseif( j == ubound(f,3) ) then
               ord = ord + ArrayOrder( i, 1 ) + ArrayOrder( i, ubound(f,3)-1 )
            else
               ord = ord + ArrayOrder( i, j+1 ) + ArrayOrder( i, j-1 )
            endif

            
            do na = 0, ubound(f,1)

               ! chop
               if( abs( real( f(na,i,j) ) ) < chopCutoff ) then
                  f(na,i,j) = ( f(na,i,j) - conjg(f(na,i,j)) )/2.0_dp
               endif
               if( abs( aimag( f(na,i,j) ) ) < chopCutoff ) then
                  f(na,i,j) = ( f(na,i,j) + conjg(f(na,i,j)) )/2.0_dp
               endif

               k(na,i,j) = mlt * ( real(na,dp) * ( mu - u/2.0_dp * (real(na,dp) - 1.0_dp) ) * f(na,i,j)

               if( na < ubound(f,1) ) then
                  k(na,i,j) = k(na,i,j) + mlt * ja * ord * SQRT( real(na+1,dp) ) * f(na+1,i,j)
               endif
               if( na > 0 ) then
                  k(na,i,j) = k(na,i,j) + mlt * ja * ord * SQRT( real(na,dp) ) * f(na-1,i,j)
               endif

            enddo ! na

         enddo ! i
      enddo ! j

   end subroutine

   subroutine c3d_solve( f, kv, ArrayOrder, ja, u, mu, mode )

      implicit none
      complex( kind = dp ), intent(inout), allocatable :: f(:,:,:,:)
      complex( kind = dp ), intent(inout), allocatable :: kv(:,:,:,:)
      complex( kind = dp ), intent(inout), allocatable :: ArrayOrder(:,:,:)
      real( kind = dp ), intent(in) :: ja, u, mu
      character(len=2), intent(in) :: mode

      ! local variables
      integer i, j, k, na
      complex( kind = dp ) :: ord, mlt

      if (mode .eq. 'EV') then
         mlt = im
      else if (mode .eq. 'GS') then
         mlt = re
      else
         mlt = (0.0_dp, 0.0_dp) 
         write(*,*)  "Wrong mode in SOLVE subroutine!"
      endif

      ! collect order parameters
      do i = 1, ubound(f,2)
         do j = 1, ubound(f,3)
            do k = 1, ubound(f,4)
         
               ArrayOrder(i,j,k) = 0.0_dp*re

               do na = 0, ubound(f,1)
               
                  if( na > 0 ) then
                     ArrayOrder(i,j,k) = ArrayOrder(i,j,k) + SQRT( real(na,dp) ) * conjg( f(na-1,i,j,k) ) * f(na,i,j,k)
                  endif

               enddo

            enddo
         enddo
      enddo

      ! solve coupled system of equtions
      do i = 1, ubound(f,2)
         do j = 1, ubound(f,3)
            do k = 1, ubound(f,4)


               ! calculate nearest-neighbour order parameter (periodic boundary conditions)
               ord = (0.0_dp, 0.0_dp)

               if( i == 1 ) then
                  ord = ord + ArrayOrder( 2, j, k ) + ArrayOrder( ubound(f,2), j, k )
               elseif( i == ubound(f,2) ) then
                  ord = ord + ArrayOrder( 1, j, k ) + ArrayOrder( ubound(f,2)-1, j, k )
               else
                  ord = ord + ArrayOrder( i+1, j, k ) + ArrayOrder( i-1, j, k )
               endif

               if( j == 1 ) then
                  ord = ord + ArrayOrder( i, 2, k ) + ArrayOrder( i, ubound(f,3), k )
               elseif( j == ubound(f,3) ) then
                  ord = ord + ArrayOrder( i, 1, k ) + ArrayOrder( i, ubound(f,3)-1, k )
               else
                  ord = ord + ArrayOrder( i, j+1, k ) + ArrayOrder( i, j-1, k )
               endif

               if( k == 1 ) then
                  ord = ord + ArrayOrder( i, j, 2 ) + ArrayOrder( i, j, ubound(f,4) )
               elseif( k == ubound(f,4) ) then
                  ord = ord + ArrayOrder( i, j, 1 ) + ArrayOrder( i, j, ubound(f,4)-1 )
               else
                  ord = ord + ArrayOrder( i, j, k+1 ) + ArrayOrder( i, j, k-1 )
               endif
                                      
               do na = 0, ubound(f,1)

                  ! chop
                  if( abs( real( f(na,i,j,k) ) ) < chopCutoff ) then
                     f(na,i,j,k) = ( f(na,i,j,k) - conjg(f(na,i,j,k)) )/2.0_dp
                  endif
                  if( abs( aimag( f(na,i,j,k) ) ) < chopCutoff ) then
                     f(na,i,j,k) = ( f(na,i,j,k) + conjg(f(na,i,j,k)) )/2.0_dp
                  endif

                  kv(na,i,j,k) = mlt * ( real(na,dp) * ( mu - u/2.0_dp * (real(na,dp) - 1.0_dp) ) * f(na,i,j,k)

                  if( na < ubound(f,1) ) then
                     kv(na,i,j,k) = kv(na,i,j,k) + mlt * ja * ord * SQRT( real(na+1,dp) ) * f(na+1,i,j,k)
                  endif
                  if( na > 0 ) then
                     kv(na,i,j,k) = kv(na,i,j,k) + mlt * ja * ord * SQRT( real(na,dp) ) * f(na-1,i,j,k)
                  endif

               enddo ! na

            enddo ! i
         enddo ! j
      enddo ! k

   end subroutine

!  | ---------------------------------------------- |
!  | Propagate in imaginary time to find the ground | 
!  | state configuration                            |
!  | ---------------------------------------------- |
   subroutine c1d_GroundStateNC( f, ja, u, mn )

      implicit none
      complex( kind = dp ), intent(inout), allocatable :: f(:,:)
      real( kind = dp ), intent(in) :: ja, u, mn

      ! | -------------------- |
      ! | For Runge-Kutta step |
      ! | -------------------- |
      complex( kind = dp ), allocatable :: k1(:,:)
      complex( kind = dp ), allocatable :: k2(:,:)
      complex( kind = dp ), allocatable :: k3(:,:)
      complex( kind = dp ), allocatable :: k4(:,:)
      complex( kind = dp ), allocatable :: tmp(:,:)
      complex( kind = dp ), allocatable :: fpom(:,:)
      complex( kind = dp ), allocatable :: ArrayOrder(:)
      
      real( kind = dp ) :: error, ierr, mu
      integer i, na, cnt, info

      allocate( k1(lbound(f,1):ubound(f,1), lbound(f,2):ubound(f,2)) )
      allocate( k2(lbound(f,1):ubound(f,1), lbound(f,2):ubound(f,2)) )
      allocate( k3(lbound(f,1):ubound(f,1), lbound(f,2):ubound(f,2)) )
      allocate( k4(lbound(f,1):ubound(f,1), lbound(f,2):ubound(f,2)) )
      allocate( tmp(lbound(f,1):ubound(f,1), lbound(f,2):ubound(f,2)) )
      allocate( fpom(lbound(f,1):ubound(f,1), lbound(f,2):ubound(f,2)) )
      allocate( ArrayOrder( size(f,2) ) )

      fpom(:,:) = f(:,:)
      cnt = 0
      mu = 0.0_dp

      do while( .true. )

         cnt = cnt + 1 
 
         ! 1st Runge-Kutta step
         call c1d_solve( f, k1, ArrayOrder, ja, u, mu, 'GS' )
         ! update vector to the 2nd step
         tmp(:,:) = f(:,:) + dt/2.0_dp * k1(:,:)
         ! 2nd Runge-Kutta step
         call c1d_solve( tmp, k2, ArrayOrder, ja, u, mu, 'GS' )
         ! update vector for the 3rd step
         tmp(:,:) = f(:,:) + dt/2.0_dp * k2(:,:)
         ! 3rd Runge-Kutta step
         call c1d_solve( tmp, k3, ArrayOrder, ja, u, mu, 'GS' )
         ! update vector for the 4th step
         tmp(:,:) = f(:,:) + dt * k3(:,:)
         ! 4th Runge-Kutta step
         call c1d_solve( tmp, k4, ArrayOrder, ja, u, mu, 'GS' )
         ! update vector
         f(:,:) = f(:,:) + dt/6.0_dp * ( k1(:,:) + 2.0_dp *&
                    k2(:,:) + 2.0_dp * k3(:,:) + k4(:,:) )

         ! normalize vector
         call const_mean( f, mn, info )
         
         if( info /= 0 ) then
            ! suppress stdout if quietON
            if( .not. quietON ) print *, 'const_mean info = ', info
            exit
         endif

         ! calculate error
         if( cnt .eq. stepsForJudge ) then

            cnt = 0
            error = 0.0_dp
            do i = 1, ubound(f,2)
               do na = 0, ubound(f,1)

                  ! absolute error
                  ierr = abs( ( f(na,i) - fpom(na,i) ) )
                  if (ierr .gt. error) then
                     error = ierr
                  endif

               enddo
            enddo

            if( error < convCriterion ) then
               exit
            endif

            if( .not. quietON ) print *, 'GroundStateNC: error = ', error
            fpom(:,:) = f(:,:)

         endif

      enddo

   end subroutine

   subroutine c2d_GroundStateNC( f, ja, u, mn )

      implicit none
      complex( kind = dp ), intent(inout), allocatable :: f(:,:,:)
      real( kind = dp ), intent(in) :: ja, u, mn

      ! | -------------------- |
      ! | For Runge-Kutta step |
      ! | -------------------- |
      complex( kind = dp ), allocatable :: k1(:,:,:)
      complex( kind = dp ), allocatable :: k2(:,:,:)
      complex( kind = dp ), allocatable :: k3(:,:,:)
      complex( kind = dp ), allocatable :: k4(:,:,:)
      complex( kind = dp ), allocatable :: tmp(:,:,:)
      complex( kind = dp ), allocatable :: fpom(:,:,:)
      complex( kind = dp ), allocatable :: ArrayOrder(:,:)
      
      real( kind = dp ) :: error, ierr, mu
      integer i, j, na, cnt, info

      allocate( k1(lbound(f,1):ubound(f,1), lbound(f,2):ubound(f,2), lbound(f,3):ubound(f,3)) )
      allocate( k2(lbound(f,1):ubound(f,1), lbound(f,2):ubound(f,2), lbound(f,3):ubound(f,3)) )
      allocate( k3(lbound(f,1):ubound(f,1), lbound(f,2):ubound(f,2), lbound(f,3):ubound(f,3)) )
      allocate( k4(lbound(f,1):ubound(f,1), lbound(f,2):ubound(f,2), lbound(f,3):ubound(f,3)) )
      allocate( tmp(lbound(f,1):ubound(f,1), lbound(f,2):ubound(f,2), lbound(f,3):ubound(f,3)) )
      allocate( fpom(lbound(f,1):ubound(f,1), lbound(f,2):ubound(f,2), lbound(f,3):ubound(f,3)) )
      allocate( ArrayOrder( size(f,2), size(f,3) ) )

      fpom(:,:,:) = f(:,:,:)
      cnt = 0
      mu = 0.0_dp

      do while( .true. )

         cnt = cnt + 1 
 
         ! 1st Runge-Kutta step
         call c2d_solve( f, k1, ArrayOrder, ja, u, mu, 'GS' )
         ! update vector to the 2nd step
         tmp(:,:,:) = f(:,:,:) + dt/2.0_dp * k1(:,:,:)
         ! 2nd Runge-Kutta step
         call c2d_solve( tmp, k2, ArrayOrder, ja, u, mu, 'GS' )
         ! update vector for the 3rd step
         tmp(:,:,:) = f(:,:,:) + dt/2.0_dp * k2(:,:,:)
         ! 3rd Runge-Kutta step
         call c2d_solve( tmp, k3, ArrayOrder, ja, u, mu, 'GS' )
         ! update vector for the 4th step
         tmp(:,:,:) = f(:,:,:) + dt * k3(:,:,:)
         ! 4th Runge-Kutta step
         call c2d_solve( tmp, k4, ArrayOrder, ja, u, mu, 'GS' )
         ! update vector
         f(:,:,:) = f(:,:,:) + dt/6.0_dp * ( k1(:,:,:) + 2.0_dp *&
                      k2(:,:,:) + 2.0_dp * k3(:,:,:) + k4(:,:,:) )

         ! normalize vector
         call const_mean( f, mn, info )
         
         if( info /= 0 ) then
            ! suppress stdout if quietON
            if( .not. quietON ) print *, 'const_mean info = ', info
            exit
         endif

         ! calculate error
         if( cnt .eq. stepsForJudge ) then

            cnt = 0
            error = 0.0_dp
            do i = 1, ubound(f,2)
               do j = 1, ubound(f,3)
                  do na = 0, ubound(f,1)
                     ierr = abs( ( f(na,i,j) - fpom(na,i,j) ) )
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

      enddo

   end subroutine

   subroutine c3d_GroundStateNC( f, ja, u, mn )

      implicit none
      complex( kind = dp ), intent(inout), allocatable :: f(:,:,:,:)
      real( kind = dp ), intent(in) :: ja, u, mn

      ! | -------------------- |
      ! | For Runge-Kutta step |
      ! | -------------------- |
      complex( kind = dp ), allocatable :: k1(:,:,:,:)
      complex( kind = dp ), allocatable :: k2(:,:,:,:)
      complex( kind = dp ), allocatable :: k3(:,:,:,:)
      complex( kind = dp ), allocatable :: k4(:,:,:,:)
      complex( kind = dp ), allocatable :: tmp(:,:,:,:)
      complex( kind = dp ), allocatable :: fpom(:,:,:,:)
      complex( kind = dp ), allocatable :: ArrayOrder(:,:,:)
      
      real( kind = dp ) :: error, ierr, mu
      integer i, j, k, na, cnt, info

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
      allocate( ArrayOrder( size(f,2), size(f,3), size(f,4) ) )
      

      fpom(:,:,:,:) = f(:,:,:,:)
      cnt = 0
      mu = 0.0_dp

      do while( .true. )

         cnt = cnt + 1 
 
         ! 1st Runge-Kutta step
         call c3d_solve( f, k1, ArrayOrder, ja, u, mu, 'GS' )
         ! update vector to the 2nd step
         tmp(:,:,:,:) = f(:,:,:,:) + dt/2.0_dp * k1(:,:,:,:)
         ! 2nd Runge-Kutta step
         call c3d_solve( tmp, k2, ArrayOrder, ja, u, mu, 'GS' )
         ! update vector for the 3rd step
         tmp(:,:,:,:) = f(:,:,:,:) + dt/2.0_dp * k2(:,:,:,:)
         ! 3rd Runge-Kutta step
         call c3d_solve( tmp, k3, ArrayOrder, ja, u, mu, 'GS' )
         ! update vector for the 4th step
         tmp(:,:,:,:) = f(:,:,:,:) + dt * k3(:,:,:,:)
         ! 4th Runge-Kutta step
         call c3d_solve( tmp, k4, ArrayOrder, ja, u, mu, 'GS' )
         ! update vector
         f(:,:,:,:) = f(:,:,:,:) + dt/6.0_dp * ( k1(:,:,:,:) + 2.0_dp *&
                        k2(:,:,:,:) + 2.0_dp * k3(:,:,:,:) + k4(:,:,:,:) )

         ! normalize vector
         call const_mean( f, mn, info )
         
         if( info /= 0 ) then
            ! suppress stdout if quietON
            if( .not. quietON ) print *, 'const_mean info = ', info
            exit
         endif


         ! calculate error
         if( cnt .eq. stepsForJudge ) then

            cnt = 0
            error = 0.0_dp
            do i = 1, ubound(f,2)
               do j = 1, ubound(f,3)
                  do k = 1, ubound(f,4)
                     do na = 0, ubound(f,1)
                        ierr = abs( ( f(na,i,j,k) - fpom(na,i,j,k) ) )
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

!  | ---------------------------------------------- |
!  | Propagate in imaginary time to find the ground | 
!  | state configuration                            |
!  | ---------------------------------------------- |
   subroutine c1d_GroundState( f, ja, u, mu )

      implicit none
      complex( kind = dp ), intent(inout), allocatable :: f(:,:)
      real( kind = dp ), intent(in) :: ja, u, mu

      ! | -------------------- |
      ! | For Runge-Kutta step |
      ! | -------------------- |
      complex( kind = dp ), allocatable :: k1(:,:)
      complex( kind = dp ), allocatable :: k2(:,:)
      complex( kind = dp ), allocatable :: k3(:,:)
      complex( kind = dp ), allocatable :: k4(:,:)
      complex( kind = dp ), allocatable :: tmp(:,:)
      complex( kind = dp ), allocatable :: fpom(:,:)
      complex( kind = dp ), allocatable :: ArrayOrder(:)
      
      real( kind = dp ) :: error, ierr
      integer i, na, cnt, info

      allocate( k1(lbound(f,1):ubound(f,1), lbound(f,2):ubound(f,2)) )
      allocate( k2(lbound(f,1):ubound(f,1), lbound(f,2):ubound(f,2)) )
      allocate( k3(lbound(f,1):ubound(f,1), lbound(f,2):ubound(f,2)) )
      allocate( k4(lbound(f,1):ubound(f,1), lbound(f,2):ubound(f,2)) )
      allocate( tmp(lbound(f,1):ubound(f,1), lbound(f,2):ubound(f,2)) )
      allocate( fpom(lbound(f,1):ubound(f,1), lbound(f,2):ubound(f,2)) )
      allocate( ArrayOrder( size(f,2) ) )

      fpom(:,:) = f(:,:)
      cnt = 0

      do while( .true. )

         cnt = cnt + 1 
 
         ! 1st Runge-Kutta step
         call c1d_solve( f, k1, ArrayOrder, ja, u, mu, 'GS' )
         ! update vector to the 2nd step
         tmp(:,:) = f(:,:) + dt/2.0_dp * k1(:,:)
         ! 2nd Runge-Kutta step
         call c1d_solve( tmp, k2, ArrayOrder, ja, u, mu, 'GS' )
         ! update vector for the 3rd step
         tmp(:,:) = f(:,:) + dt/2.0_dp * k2(:,:)
         ! 3rd Runge-Kutta step
         call c1d_solve( tmp, k3, ArrayOrder, ja, u, mu, 'GS' )
         ! update vector for the 4th step
         tmp(:,:) = f(:,:) + dt * k3(:,:)
         ! 4th Runge-Kutta step
         call c1d_solve( tmp, k4, ArrayOrder, ja, u, mu, 'GS' )
         ! update vector
         f(:,:) = f(:,:) + dt/6.0_dp * ( k1(:,:) + 2.0_dp *&
                    k2(:,:) + 2.0_dp * k3(:,:) + k4(:,:) )

         ! normalize vector
         call normalize( f )
         
         if( info /= 0 ) then
            ! suppress stdout if quietON
            if( .not. quietON ) print *, 'const_mean info = ', info
            exit
         endif

         ! calculate error
         if( cnt .eq. stepsForJudge ) then

            cnt = 0
            error = 0.0_dp
            do i = 1, ubound(f,2)
               do na = 0, ubound(f,1)

                  ! absolute error
                  ierr = abs( ( f(na,i) - fpom(na,i) ) )
                  if (ierr .gt. error) then
                     error = ierr
                  endif

               enddo
            enddo

            if( error < convCriterion ) then
               exit
            endif

            if( .not. quietON ) print *, 'GroundStateNC: error = ', error
            fpom(:,:) = f(:,:)

         endif

      enddo

   end subroutine

   subroutine c2d_GroundState( f, ja, u, mu )

      implicit none
      complex( kind = dp ), intent(inout), allocatable :: f(:,:,:)
      real( kind = dp ), intent(in) :: ja, u, mu

      ! | -------------------- |
      ! | For Runge-Kutta step |
      ! | -------------------- |
      complex( kind = dp ), allocatable :: k1(:,:,:)
      complex( kind = dp ), allocatable :: k2(:,:,:)
      complex( kind = dp ), allocatable :: k3(:,:,:)
      complex( kind = dp ), allocatable :: k4(:,:,:)
      complex( kind = dp ), allocatable :: tmp(:,:,:)
      complex( kind = dp ), allocatable :: fpom(:,:,:)
      complex( kind = dp ), allocatable :: ArrayOrder(:,:)
      
      real( kind = dp ) :: error, ierr
      integer i, j, na, cnt, info

      allocate( k1(lbound(f,1):ubound(f,1), lbound(f,2):ubound(f,2), lbound(f,3):ubound(f,3)) )
      allocate( k2(lbound(f,1):ubound(f,1), lbound(f,2):ubound(f,2), lbound(f,3):ubound(f,3)) )
      allocate( k3(lbound(f,1):ubound(f,1), lbound(f,2):ubound(f,2), lbound(f,3):ubound(f,3)) )
      allocate( k4(lbound(f,1):ubound(f,1), lbound(f,2):ubound(f,2), lbound(f,3):ubound(f,3)) )
      allocate( tmp(lbound(f,1):ubound(f,1), lbound(f,2):ubound(f,2), lbound(f,3):ubound(f,3)) )
      allocate( fpom(lbound(f,1):ubound(f,1), lbound(f,2):ubound(f,2), lbound(f,3):ubound(f,3)) )
      allocate( ArrayOrder( size(f,2), size(f,3) ) )

      fpom(:,:,:) = f(:,:,:)
      cnt = 0

      do while( .true. )

         cnt = cnt + 1 
 
         ! 1st Runge-Kutta step
         call c2d_solve( f, k1, ArrayOrder, ja, u, mu, 'GS' )
         ! update vector to the 2nd step
         tmp(:,:,:) = f(:,:,:) + dt/2.0_dp * k1(:,:,:)
         ! 2nd Runge-Kutta step
         call c2d_solve( tmp, k2, ArrayOrder, ja, u, mu, 'GS' )
         ! update vector for the 3rd step
         tmp(:,:,:) = f(:,:,:) + dt/2.0_dp * k2(:,:,:)
         ! 3rd Runge-Kutta step
         call c2d_solve( tmp, k3, ArrayOrder, ja, u, mu, 'GS' )
         ! update vector for the 4th step
         tmp(:,:,:) = f(:,:,:) + dt * k3(:,:,:)
         ! 4th Runge-Kutta step
         call c2d_solve( tmp, k4, ArrayOrder, ja, u, mu, 'GS' )
         ! update vector
         f(:,:,:) = f(:,:,:) + dt/6.0_dp * ( k1(:,:,:) + 2.0_dp *&
                      k2(:,:,:) + 2.0_dp * k3(:,:,:) + k4(:,:,:) )

         ! normalize vector
         call normalize( f )
         
         if( info /= 0 ) then
            ! suppress stdout if quietON
            if( .not. quietON ) print *, 'const_mean info = ', info
            exit
         endif

         ! calculate error
         if( cnt .eq. stepsForJudge ) then

            cnt = 0
            error = 0.0_dp
            do i = 1, ubound(f,2)
               do j = 1, ubound(f,3)
                  do na = 0, ubound(f,1)
                     ierr = abs( ( f(na,i,j) - fpom(na,i,j) ) )
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

      enddo

   end subroutine

   subroutine c3d_GroundState( f, ja, u, mu )

      implicit none
      complex( kind = dp ), intent(inout), allocatable :: f(:,:,:,:)
      real( kind = dp ), intent(in) :: ja, u, mu

      ! | -------------------- |
      ! | For Runge-Kutta step |
      ! | -------------------- |
      complex( kind = dp ), allocatable :: k1(:,:,:,:)
      complex( kind = dp ), allocatable :: k2(:,:,:,:)
      complex( kind = dp ), allocatable :: k3(:,:,:,:)
      complex( kind = dp ), allocatable :: k4(:,:,:,:)
      complex( kind = dp ), allocatable :: tmp(:,:,:,:)
      complex( kind = dp ), allocatable :: fpom(:,:,:,:)
      complex( kind = dp ), allocatable :: ArrayOrder(:,:,:)
      
      real( kind = dp ) :: error, ierr, mu
      integer i, j, k, na, cnt, info

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
      allocate( ArrayOrder( size(f,2), size(f,3), size(f,4) ) )
      

      fpom(:,:,:,:) = f(:,:,:,:)
      cnt = 0

      do while( .true. )

         cnt = cnt + 1 
 
         ! 1st Runge-Kutta step
         call c3d_solve( f, k1, ArrayOrder, ja, u, mu, 'GS' )
         ! update vector to the 2nd step
         tmp(:,:,:,:) = f(:,:,:,:) + dt/2.0_dp * k1(:,:,:,:)
         ! 2nd Runge-Kutta step
         call c3d_solve( tmp, k2, ArrayOrder, ja, u, mu, 'GS' )
         ! update vector for the 3rd step
         tmp(:,:,:,:) = f(:,:,:,:) + dt/2.0_dp * k2(:,:,:,:)
         ! 3rd Runge-Kutta step
         call c3d_solve( tmp, k3, ArrayOrder, ja, u, mu, 'GS' )
         ! update vector for the 4th step
         tmp(:,:,:,:) = f(:,:,:,:) + dt * k3(:,:,:,:)
         ! 4th Runge-Kutta step
         call c3d_solve( tmp, k4, ArrayOrder, ja, u, mu, 'GS' )
         ! update vector
         f(:,:,:,:) = f(:,:,:,:) + dt/6.0_dp * ( k1(:,:,:,:) + 2.0_dp *&
                        k2(:,:,:,:) + 2.0_dp * k3(:,:,:,:) + k4(:,:,:,:) )

         ! normalize vector
         call normalize( f )
         
         if( info /= 0 ) then
            ! suppress stdout if quietON
            if( .not. quietON ) print *, 'const_mean info = ', info
            exit
         endif


         ! calculate error
         if( cnt .eq. stepsForJudge ) then

            cnt = 0
            error = 0.0_dp
            do i = 1, ubound(f,2)
               do j = 1, ubound(f,3)
                  do k = 1, ubound(f,4)
                     do na = 0, ubound(f,1)
                        ierr = abs( ( f(na,i,j,k) - fpom(na,i,j,k) ) )
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

!  | ---------------------------------------------- |
!  | Initialize uniformly with unit normalization   | 
!  | ---------------------------------------------- |
   subroutine c1d_InitUniform( f )
 
      implicit none
      complex( kind = dp ), intent(inout), allocatable :: f(:,:)

      f(:,:) = (1.0_dp, 1.0_dp)
      call normalize( f )

   end subroutine

   subroutine c2d_InitUniform( f )
 
      implicit none
      complex( kind = dp ), intent(inout), allocatable :: f(:,:,:)

      f(:,:,:) = (1.0_dp, 1.0_dp)
      call normalize( f )

   end subroutine

   subroutine c3d_InitUniform( f )
 
      implicit none
      complex( kind = dp ), intent(inout), allocatable :: f(:,:,:,:)

      f(:,:,:,:) = (1.0_dp, 1.0_dp)
      call normalize( f )

   end subroutine

!  | ---------------------------------------------- |
!  | Initialize randomly  with unit normalization   | 
!  | ---------------------------------------------- |
   subroutine c1d_InitRandom( f )
 
      implicit none
      complex( kind = dp ), intent(inout), allocatable :: f(:,:)

      integer :: i, na
      real :: ran(2)

      call init_random_seed()

      do na = 0, ubound(f,1)
         do i = 1, ubound(f,2)
            call random_number( ran )
            f(na,i) = ran(1)*re + im*ran(2) 
         enddo
      enddo

      call normalize( f )

   end subroutine

   subroutine c2d_InitRandom( f )
 
      implicit none
      complex( kind = dp ), intent(inout), allocatable :: f(:,:,:)

      integer :: i, j, na
      real :: ran(2)

      do na = 0, ubound(f,1)
         do i = 1, ubound(f,2)
            do j = 1, ubound(f,3)
               call random_number( ran )
               f(na,i,j) = ran(1)*re + im*ran(2)
            enddo
         enddo
      enddo

      call normalize( f )

   end subroutine

   subroutine c3d_InitRandom( f )
 
      implicit none
      complex( kind = dp ), intent(inout), allocatable :: f(:,:,:,:)

      integer :: i, j, k, na
      real :: ran(2)

      do na = 0, ubound(f,1)
         do i = 1, ubound(f,2)
            do j = 1, ubound(f,3)
               do k = 1, ubound(f,4)
                  call random_number( ran )
                  f(na,i,j,k) = ran(1)*re + im*ran(2)
               enddo
            enddo
         enddo
      enddo

      call normalize( f )

   end subroutine
!  | ---------------------------------------------- |
!  | Initialize uniformly with unit normalization   |
!  | and fixed mean number of particles in each     |
!  | ---------------------------------------------- |
   subroutine c1d_InitUniformNC( f, mn )

      implicit none
      complex( kind = dp ), intent(inout), allocatable :: f(:,:)
      real( kind = dp ), intent(in) :: mn
      integer info

      f(:,:) = (1.0_dp, 1.0_dp)
      call const_mean( f, mn, info )

   end subroutine

   subroutine c2d_InitUniformNC( f, mn )

      implicit none
      complex( kind = dp ), intent(inout), allocatable :: f(:,:,:)
      real( kind = dp ), intent(in) :: mn
      integer info

      f(:,:,:) = (1.0_dp, 1.0_dp)
      call const_mean( f, mn, info )

   end subroutine

   subroutine c3d_InitUniformNC( f, mn )

      implicit none
      complex( kind = dp ), intent(inout), allocatable :: f(:,:,:,:)
      real( kind = dp ), intent(in) :: mn
      integer info

      f(:,:,:,:) = (1.0_dp, 1.0_dp)
      call const_mean( f, mn, info )

   end subroutine

!  | ---------------------------------------------- |
!  | Initialize randomly with unit normalization    |
!  | and fixed mean number of particles in each     |
!  | ---------------------------------------------- |
   subroutine c1d_InitRandomNC( f, mn )

      implicit none
      complex( kind = dp ), intent(inout), allocatable :: f(:,:)
      real( kind = dp ), intent(in) :: mn
      integer info

      integer :: i, na
      real :: ran(2)

      call init_random_seed()

      do na = 0, ubound(f,1)
         do i = 1, ubound(f,2)
            call random_number( ran )
            f(na,i) = ran(1)*re + im*ran(2) 
         enddo
      enddo
      
      call const_mean( f, mn, info )

   end subroutine

   subroutine c2d_InitRandomNC( f, mn )

      implicit none
      complex( kind = dp ), intent(inout), allocatable :: f(:,:,:)
      real( kind = dp ), intent(in) :: mn
      integer info

      integer :: i, j, na
      real :: ran(2)

      do na = 0, ubound(f,1)
         do i = 1, ubound(f,2)
            do j = 1, ubound(f,3)
               call random_number( ran )
               f(na,i,j) = ran(1)*re + im*ran(2)
            enddo
         enddo
      enddo

      call const_mean( f, mn, info )

   end subroutine

   subroutine c3d_InitRandomNC( f, mn )

      implicit none
      complex( kind = dp ), intent(inout), allocatable :: f(:,:,:,:)
      real( kind = dp ), intent(in) :: mn
      integer info

      integer :: i, j, k, na
      real :: ran(2)

      do na = 0, ubound(f,1)
         do i = 1, ubound(f,2)
            do j = 1, ubound(f,3)
               do k = 1, ubound(f,4)
                  call random_number( ran )
                  f(na,i,j,k) = ran(1)*re + im*ran(2)
               enddo
            enddo
         enddo
      enddo

      call const_mean( f, mn, info )

   end subroutine

! | ------------------------------------------------ |
! | simple cpu time seeding for random number gen    |
! | ------------------------------------------------ |
   subroutine init_random_seed()
      implicit none
      integer, parameter :: int64 = selected_int_kind(16)
      integer, allocatable :: seed(:)
      integer :: i, n, dt(8)
      integer(int64) :: t

      call random_seed( size = n )
      allocate( seed(n) )

      call date_and_time( values = dt )
      t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
        + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
        + dt(3) * 24_int64 * 60 * 60 * 1000 &
        + dt(5) * 60 * 60 * 1000 &
        + dt(6) * 60 * 1000 + dt(7) * 1000 &
        + dt(8)

      do i = 1, n
         seed(i) = lcg(t)
      enddo

      call random_seed( put = seed )
      contains

      function lcg(s)
         integer :: lcg
         integer(int64) :: s
         if( s == 0 ) then
            s = 104729
         else
            s = mod( s, 4294967296_int64 )
         endif
         s = mod( s * 279470273_int64, 4294967291_int64 )
         lcg = int( mod( s, int(huge(0), int64 ) ), kind(0) )
      end function lcg
   end subroutine

end module
