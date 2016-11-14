module squeezing

   use parameters

   private
   public SpinMatrix, SpinSqueezing

   interface SpinMatrix
      module procedure c1d_SpinMatrix
      module procedure c2d_SpinMatrix
      module procedure c3d_SpinMatrix
   end interface

   interface SpinSqueezing
      module procedure c1d_SpinSqueezing
      module procedure c2d_SpinSqueezing
      module procedure c3d_SpinSqueezing
   end interface

   contains

   ! -------------------------------------------------------------------------------------------
   ! Spin correlation matrix 
   ! -------------------------------------------------------------------------------------------

   function c1d_SpinMatrix( f ) result mat

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:,:)
      complex( kind = dp ) mat(4,3)

      ! local variables
      integer i, na, nb
      complex( kind = dp ) sumab, sumba, sumaa, sumbb
      complex( kind = dp ) sum_aabb, sum_bbaa, sum_abab, &
                           sum_aaab, sum_abaa, sum_abbb, &
                           sum_bbba, sum_aaaa, sum_bbbb, &
                           sx, sy, sz
      real( kind = dp ) Msite

      ! number of lattice sites
      Msite = real( size(f,3), dp )

      sumab = 0.0_dp
      sumba = 0.0_dp
      sumaa = 0.0_dp
      sumbb = 0.0_dp
      sum_aabb = 0.0_dp
      sum_bbaa = 0.0_dp
      sum_abab = 0.0_dp
      sum_aaab = 0.0_dp
      sum_abaa = 0.0_dp
      sum_abbb = 0.0_dp
      sum_bbba = 0.0_dp
      sum_aaaa = 0.0_dp
      sum_bbbb = 0.0_dp

      mat(:,:) = 0.0_dp

      ! summation over lattice sites
      do i = 1, ubound(f,3)

         do nb = 0, ubound(f,2)
            do na = 0, ubound(f,1)

               ! -----------------------------------------------------------------------------------
               ! spin expectation value
               ! -----------------------------------------------------------------------------------
               sumaa    = sumaa + real(n,dp) * conjg( f(na,nb,i) ) * f(na,nb,i)
               sumbb    = sumbb + real(m,dp) * conjg( f(na,nb,i) ) * f(na,nb,i)
               sum_abab = sum_abab + 2.0_dp * real(n,dp) * real(m,dp) * conjg( f(na,nb,i) ) * f(na,nb,i)
               sum_aaaa = sum_aaaa + real(na*(na-1),dp) * conjg( f(na,nb,i) ) * f(na,nb,i)
               sum_bbbb = sum_bbbb + real(nb*(nb-1),dp) * conjg( f(na,nb,i) ) * f(na,nb,i)

               if( na > 0 ) then
                  sumab    = sumab + sqrt( real(na*(nb+1),dp) ) * 
                                     conjg( f(na,nb,i) ) * f(na-1,nb+1,i)
                  sum_aaab = sum_aaab + sqrt( real(na*(na-1)*(na-1)*(nb+1),dp) ) * &
                                        conjg( f(na,nb,i) ) * f(na-1,nb+1,i)
                  sum_abbb = sum_abbb + sqrt( real(na*nb*nb*(nb+1),dp) ) * &
                                        conjg( f(na,nb,i) ) * f(na-1,nb+1,i)
               endif

               if( nb > 0 ) then
                  sumba    = sumba + sqrt( real(nb*(na+1),dp) ) * 
                                     conjg( f(na,nb,i) ) * f(na+1,nb-1,i)
                  sum_abaa = sum_abaa + sqrt( real(na*na*(na+1)*nb,dp) ) * &
                                        conjg( f(na,nb,i) ) * f(na+1,nb-1,i)
                  sum_bbba = sum_bbba + sqrt( real((na+1)*nb*(nb-1)*(nb-1),dp) ) * &
                                        conjg( f(na,nb,i) ) * f(na+1,nb-1,i)
               endif

               if( na > 1 ) then
                  sum_aabb = sum_aabb + sqrt( real(na*(na-1)*(nb+1)*(nb+2),dp) ) * &
                                        conjg( f(na,nb,i) ) * f(na-2,nb+2,i)
               endif

               if( nb > 1 ) then
                  sum_bbaa = sum_bbaa + sqrt( real(nb*(nb-1)*(na+1)*(na+2),dp) ) * &
                                        conjg( f(na,nb,i) ) * f(na+2,nb-2,i)
               endif

            enddo
         enddo

         sx = 0.5_dp * (sumab + sumba)
         sy = 0.5_dp * im * (sumba - sumab)
         sz = 0.5_dp * (sumaa - sumbb)

         mat(1,1) = mat(1,1) + sx
         mat(1,2) = mat(1,2) + sy
         mat(1,3) = mat(1,3) + sz

         ! Jx*Jx
         mat(2,1) = mat(2,1) + 0.25_dp * (sumaa + sumbb + sum_aabb + sum_bbaa + sum_abab)
         mat(2,1) = mat(2,1) - sx * sx
         ! Jx*Jy
         mat(2,2) = mat(2,2) + 0.25_dp * im * (sumaa - sumbb - sum_aabb + sum_bbaa)
         mat(2,2) = mat(2,2) - sx * sy
         ! Jx*Jz
         mat(2,3) = mat(2,3) + 0.25_dp * (sumba - sumab + sum_aaab + sum_abaa - sum_abbbb - sum_bbba)
         mat(2,3) = mat(2,3) - sx * sz
         ! Jy*Jx
         mat(3,1) = mat(3,1) + 0.25_dp * im * (sumbb - sumaa + sum_bbaa - sum_aabb)
         mat(3,1) = mat(3,1) - sy * sz
         ! Jy*Jy
         mat(3,2) = mat(3,2) + 0.25_dp * (sumaa + sumbb - sum_aabb - sum_bbaa + sum_abab)
         mat(3,2) = mat(3,2) - sy * sy
         ! Jy*Jz 
         mat(3,3) = mat(3,3) + 0.25_dp * im * (sumab + sumba - sum_aaab + sum_abaa + sum_abbb - sum_bbba)
         mat(3,3) = mat(3,3) - sy * sz
         ! Jz*Jx
         mat(4,1) = mat(4,1) + 0.25_dp * (sumab - sumba + sum_aaab + sum_abaa - sum_abbb - sum_bbba)
         mat(4,1) = mat(4,1) - sz * sx
         ! Jz*Jy
         mat(4,2) = mat(4,2) + 0.25_dp * im * (sum_abaa + sum_abbb - sum_bbba - sum_aaab - sumab - sumba)
         mat(4,2) = mat(4,2) - sz * sy
         ! Jz*Jz
         mat(4,3) = mat(4,3) + 0.25_dp * (sumaa + sumbb + sum_aaaa + sum_bbbb - sum_abab)
         mat(4,3) = mat(4,3) - sz * sz

      enddo

      ! Jx*Jx
      mat(2,1) = mat(2,1) + mat(1,1) * mat(1,1)
      ! Jx*Jy
      mat(2,2) = mat(2,2) + mat(1,1) * mat(1,2)
      ! Jx*Jz
      mat(2,3) = mat(2,3) + mat(1,1) * mat(1,3)
      ! Jy*Jx
      mat(3,1) = mat(3,1) + mat(1,2) * mat(1,1)
      ! Jy*Jy
      mat(3,2) = mat(3,2) + mat(1,2) * mat(1,2)
      ! Jy*Jz 
      mat(3,3) = mat(3,3) + mat(1,2) * mat(1,3)
      ! Jz*Jx
      mat(4,1) = mat(4,1) + mat(1,3) * mat(1,1)
      ! Jz*Jy
      mat(4,2) = mat(4,2) + mat(1,3) * mat(1,2)
      ! Jz*Jz
      mat(4,3) = mat(4,3) + mat(1,3) * mat(1,3)


      mat(1,:) = mat(1,:) / Msite

      mat(2:4,:) = mat(2:4,:) / Msite**2

   end function

   function c2d_SpinMatrix( f ) result mat

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:,:,:)
      complex( kind = dp ) mat(4,3)

      ! local variables
      integer i, j, na, nb
      complex( kind = dp ) sumab, sumba, sumaa, sumbb
      complex( kind = dp ) sum_aabb, sum_bbaa, sum_abab, &
                           sum_aaab, sum_abaa, sum_abbb, &
                           sum_bbba, sum_aaaa, sum_bbbb, &
                           sx, sy, sz
      real( kind = dp ) Msite

      ! number of lattice sites
      Msite = real( size(f,3)*size(f,4), dp )

      sumab = 0.0_dp
      sumba = 0.0_dp
      sumaa = 0.0_dp
      sumbb = 0.0_dp
      sum_aabb = 0.0_dp
      sum_bbaa = 0.0_dp
      sum_abab = 0.0_dp
      sum_aaab = 0.0_dp
      sum_abaa = 0.0_dp
      sum_abbb = 0.0_dp
      sum_bbba = 0.0_dp
      sum_aaaa = 0.0_dp
      sum_bbbb = 0.0_dp

      mat(:,:) = 0.0_dp

      ! summation over lattice sites
      do j = 1, ubound(f,4)
      do i = 1, ubound(f,3)

         do nb = 0, ubound(f,2)
            do na = 0, ubound(f,1)

               ! -----------------------------------------------------------------------------------
               ! spin expectation value
               ! -----------------------------------------------------------------------------------
               sumaa    = sumaa + real(n,dp) * conjg( f(na,nb,i,j) ) * f(na,nb,i,j)
               sumbb    = sumbb + real(m,dp) * conjg( f(na,nb,i,j) ) * f(na,nb,i,j)
               sum_abab = sum_abab + 2.0_dp * real(n,dp) * real(m,dp) * conjg( f(na,nb,i,j) ) * f(na,nb,i,j)
               sum_aaaa = sum_aaaa + real(na*(na-1),dp) * conjg( f(na,nb,i,j) ) * f(na,nb,i,j)
               sum_bbbb = sum_bbbb + real(nb*(nb-1),dp) * conjg( f(na,nb,i,j) ) * f(na,nb,i,j)

               if( na > 0 ) then
                  sumab    = sumab + sqrt( real(na*(nb+1),dp) ) * 
                                     conjg( f(na,nb,i,j) ) * f(na-1,nb+1,i,j)
                  sum_aaab = sum_aaab + sqrt( real(na*(na-1)*(na-1)*(nb+1),dp) ) * &
                                        conjg( f(na,nb,i,j) ) * f(na-1,nb+1,i,j)
                  sum_abbb = sum_abbb + sqrt( real(na*nb*nb*(nb+1),dp) ) * &
                                        conjg( f(na,nb,i,j) ) * f(na-1,nb+1,i,j)
               endif

               if( nb > 0 ) then
                  sumba    = sumba + sqrt( real(nb*(na+1),dp) ) * 
                                     conjg( f(na,nb,i,j) ) * f(na+1,nb-1,i,j)
                  sum_abaa = sum_abaa + sqrt( real(na*na*(na+1)*nb,dp) ) * &
                                        conjg( f(na,nb,i,j) ) * f(na+1,nb-1,i,j)
                  sum_bbba = sum_bbba + sqrt( real((na+1)*nb*(nb-1)*(nb-1),dp) ) * &
                                        conjg( f(na,nb,i,j) ) * f(na+1,nb-1,i,j)
               endif

               if( na > 1 ) then
                  sum_aabb = sum_aabb + sqrt( real(na*(na-1)*(nb+1)*(nb+2),dp) ) * &
                                        conjg( f(na,nb,i,j) ) * f(na-2,nb+2,i,j)
               endif

               if( nb > 1 ) then
                  sum_bbaa = sum_bbaa + sqrt( real(nb*(nb-1)*(na+1)*(na+2),dp) ) * &
                                        conjg( f(na,nb,i,j) ) * f(na+2,nb-2,i,j)
               endif

            enddo
         enddo

         sx = 0.5_dp * (sumab + sumba)
         sy = 0.5_dp * im * (sumba - sumab)
         sz = 0.5_dp * (sumaa - sumbb)

         mat(1,1) = mat(1,1) + sx
         mat(1,2) = mat(1,2) + sy
         mat(1,3) = mat(1,3) + sz

         ! Jx*Jx
         mat(2,1) = mat(2,1) + 0.25_dp * (sumaa + sumbb + sum_aabb + sum_bbaa + sum_abab)
         mat(2,1) = mat(2,1) - sx * sx
         ! Jx*Jy
         mat(2,2) = mat(2,2) + 0.25_dp * im * (sumaa - sumbb - sum_aabb + sum_bbaa)
         mat(2,2) = mat(2,2) - sx * sy
         ! Jx*Jz
         mat(2,3) = mat(2,3) + 0.25_dp * (sumba - sumab + sum_aaab + sum_abaa - sum_abbbb - sum_bbba)
         mat(2,3) = mat(2,3) - sx * sz
         ! Jy*Jx
         mat(3,1) = mat(3,1) + 0.25_dp * im * (sumbb - sumaa + sum_bbaa - sum_aabb)
         mat(3,1) = mat(3,1) - sy * sz
         ! Jy*Jy
         mat(3,2) = mat(3,2) + 0.25_dp * (sumaa + sumbb - sum_aabb - sum_bbaa + sum_abab)
         mat(3,2) = mat(3,2) - sy * sy
         ! Jy*Jz 
         mat(3,3) = mat(3,3) + 0.25_dp * im * (sumab + sumba - sum_aaab + sum_abaa + sum_abbb - sum_bbba)
         mat(3,3) = mat(3,3) - sy * sz
         ! Jz*Jx
         mat(4,1) = mat(4,1) + 0.25_dp * (sumab - sumba + sum_aaab + sum_abaa - sum_abbb - sum_bbba)
         mat(4,1) = mat(4,1) - sz * sx
         ! Jz*Jy
         mat(4,2) = mat(4,2) + 0.25_dp * im * (sum_abaa + sum_abbb - sum_bbba - sum_aaab - sumab - sumba)
         mat(4,2) = mat(4,2) - sz * sy
         ! Jz*Jz
         mat(4,3) = mat(4,3) + 0.25_dp * (sumaa + sumbb + sum_aaaa + sum_bbbb - sum_abab)
         mat(4,3) = mat(4,3) - sz * sz

      enddo !i
      enddo !j

      ! Jx*Jx
      mat(2,1) = mat(2,1) + mat(1,1) * mat(1,1)
      ! Jx*Jy
      mat(2,2) = mat(2,2) + mat(1,1) * mat(1,2)
      ! Jx*Jz
      mat(2,3) = mat(2,3) + mat(1,1) * mat(1,3)
      ! Jy*Jx
      mat(3,1) = mat(3,1) + mat(1,2) * mat(1,1)
      ! Jy*Jy
      mat(3,2) = mat(3,2) + mat(1,2) * mat(1,2)
      ! Jy*Jz 
      mat(3,3) = mat(3,3) + mat(1,2) * mat(1,3)
      ! Jz*Jx
      mat(4,1) = mat(4,1) + mat(1,3) * mat(1,1)
      ! Jz*Jy
      mat(4,2) = mat(4,2) + mat(1,3) * mat(1,2)
      ! Jz*Jz
      mat(4,3) = mat(4,3) + mat(1,3) * mat(1,3)


      mat(1,:) = mat(1,:) / Msite

      mat(2:4,:) = mat(2:4,:) / Msite**2

   end function

   function c3d_SpinMatrix( f ) result mat

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:,:,:)
      complex( kind = dp ) mat(4,3)

      ! local variables
      integer i, j, k, na, nb
      complex( kind = dp ) sumab, sumba, sumaa, sumbb
      complex( kind = dp ) sum_aabb, sum_bbaa, sum_abab, &
                           sum_aaab, sum_abaa, sum_abbb, &
                           sum_bbba, sum_aaaa, sum_bbbb, &
                           sx, sy, sz
      real( kind = dp ) Msite

      ! number of lattice sites
      Msite = real( size(f,3)*size(f,4)*size(f,5), dp )

      sumab = 0.0_dp
      sumba = 0.0_dp
      sumaa = 0.0_dp
      sumbb = 0.0_dp
      sum_aabb = 0.0_dp
      sum_bbaa = 0.0_dp
      sum_abab = 0.0_dp
      sum_aaab = 0.0_dp
      sum_abaa = 0.0_dp
      sum_abbb = 0.0_dp
      sum_bbba = 0.0_dp
      sum_aaaa = 0.0_dp
      sum_bbbb = 0.0_dp

      mat(:,:) = 0.0_dp

      ! summation over lattice sites
      do k = 1, ubound(f,5)
      do j = 1, ubound(f,4)
      do i = 1, ubound(f,3)

         do nb = 0, ubound(f,2)
            do na = 0, ubound(f,1)

               ! -----------------------------------------------------------------------------------
               ! spin expectation value
               ! -----------------------------------------------------------------------------------
               sumaa    = sumaa + real(n,dp) * conjg( f(na,nb,i,j,k) ) * f(na,nb,i,j,k)
               sumbb    = sumbb + real(m,dp) * conjg( f(na,nb,i,j,k) ) * f(na,nb,i,j,k)
               sum_abab = sum_abab + 2.0_dp * real(n,dp) * real(m,dp) * conjg( f(na,nb,i,j,k) ) * f(na,nb,i,j,k)
               sum_aaaa = sum_aaaa + real(na*(na-1),dp) * conjg( f(na,nb,i,j,k) ) * f(na,nb,i,j,k)
               sum_bbbb = sum_bbbb + real(nb*(nb-1),dp) * conjg( f(na,nb,i,j,k) ) * f(na,nb,i,j,k)

               if( na > 0 ) then
                  sumab    = sumab + sqrt( real(na*(nb+1),dp) ) * 
                                     conjg( f(na,nb,i,j,k) ) * f(na-1,nb+1,i,j,k)
                  sum_aaab = sum_aaab + sqrt( real(na*(na-1)*(na-1)*(nb+1),dp) ) * &
                                        conjg( f(na,nb,i,j,k) ) * f(na-1,nb+1,i,j,k)
                  sum_abbb = sum_abbb + sqrt( real(na*nb*nb*(nb+1),dp) ) * &
                                        conjg( f(na,nb,i,j,k) ) * f(na-1,nb+1,i,j,k)
               endif

               if( nb > 0 ) then
                  sumba    = sumba + sqrt( real(nb*(na+1),dp) ) * 
                                     conjg( f(na,nb,i,j,k) ) * f(na+1,nb-1,i,j,k)
                  sum_abaa = sum_abaa + sqrt( real(na*na*(na+1)*nb,dp) ) * &
                                        conjg( f(na,nb,i,j,k) ) * f(na+1,nb-1,i,j,k)
                  sum_bbba = sum_bbba + sqrt( real((na+1)*nb*(nb-1)*(nb-1),dp) ) * &
                                        conjg( f(na,nb,i,j,k) ) * f(na+1,nb-1,i,j,k)
               endif

               if( na > 1 ) then
                  sum_aabb = sum_aabb + sqrt( real(na*(na-1)*(nb+1)*(nb+2),dp) ) * &
                                        conjg( f(na,nb,i,j,k) ) * f(na-2,nb+2,i,j,k)
               endif

               if( nb > 1 ) then
                  sum_bbaa = sum_bbaa + sqrt( real(nb*(nb-1)*(na+1)*(na+2),dp) ) * &
                                        conjg( f(na,nb,i,j,k) ) * f(na+2,nb-2,i,j,k)
               endif

            enddo
         enddo

         sx = 0.5_dp * (sumab + sumba)
         sy = 0.5_dp * im * (sumba - sumab)
         sz = 0.5_dp * (sumaa - sumbb)

         mat(1,1) = mat(1,1) + sx
         mat(1,2) = mat(1,2) + sy
         mat(1,3) = mat(1,3) + sz

         ! Jx*Jx
         mat(2,1) = mat(2,1) + 0.25_dp * (sumaa + sumbb + sum_aabb + sum_bbaa + sum_abab)
         mat(2,1) = mat(2,1) - sx * sx
         ! Jx*Jy
         mat(2,2) = mat(2,2) + 0.25_dp * im * (sumaa - sumbb - sum_aabb + sum_bbaa)
         mat(2,2) = mat(2,2) - sx * sy
         ! Jx*Jz
         mat(2,3) = mat(2,3) + 0.25_dp * (sumba - sumab + sum_aaab + sum_abaa - sum_abbbb - sum_bbba)
         mat(2,3) = mat(2,3) - sx * sz
         ! Jy*Jx
         mat(3,1) = mat(3,1) + 0.25_dp * im * (sumbb - sumaa + sum_bbaa - sum_aabb)
         mat(3,1) = mat(3,1) - sy * sz
         ! Jy*Jy
         mat(3,2) = mat(3,2) + 0.25_dp * (sumaa + sumbb - sum_aabb - sum_bbaa + sum_abab)
         mat(3,2) = mat(3,2) - sy * sy
         ! Jy*Jz 
         mat(3,3) = mat(3,3) + 0.25_dp * im * (sumab + sumba - sum_aaab + sum_abaa + sum_abbb - sum_bbba)
         mat(3,3) = mat(3,3) - sy * sz
         ! Jz*Jx
         mat(4,1) = mat(4,1) + 0.25_dp * (sumab - sumba + sum_aaab + sum_abaa - sum_abbb - sum_bbba)
         mat(4,1) = mat(4,1) - sz * sx
         ! Jz*Jy
         mat(4,2) = mat(4,2) + 0.25_dp * im * (sum_abaa + sum_abbb - sum_bbba - sum_aaab - sumab - sumba)
         mat(4,2) = mat(4,2) - sz * sy
         ! Jz*Jz
         mat(4,3) = mat(4,3) + 0.25_dp * (sumaa + sumbb + sum_aaaa + sum_bbbb - sum_abab)
         mat(4,3) = mat(4,3) - sz * sz

      enddo !i
      enddo !j
      enddo !k

      ! Jx*Jx
      mat(2,1) = mat(2,1) + mat(1,1) * mat(1,1)
      ! Jx*Jy
      mat(2,2) = mat(2,2) + mat(1,1) * mat(1,2)
      ! Jx*Jz
      mat(2,3) = mat(2,3) + mat(1,1) * mat(1,3)
      ! Jy*Jx
      mat(3,1) = mat(3,1) + mat(1,2) * mat(1,1)
      ! Jy*Jy
      mat(3,2) = mat(3,2) + mat(1,2) * mat(1,2)
      ! Jy*Jz 
      mat(3,3) = mat(3,3) + mat(1,2) * mat(1,3)
      ! Jz*Jx
      mat(4,1) = mat(4,1) + mat(1,3) * mat(1,1)
      ! Jz*Jy
      mat(4,2) = mat(4,2) + mat(1,3) * mat(1,2)
      ! Jz*Jz
      mat(4,3) = mat(4,3) + mat(1,3) * mat(1,3)


      mat(1,:) = mat(1,:) / Msite

      mat(2:4,:) = mat(2:4,:) / Msite**2

   end function

   ! ------------------------------------------------------------------------------------
   ! Spin squeezing functions
   ! ------------------------------------------------------------------------------------

   function c1d_SpinSqueezing( f ) result vout
      ! Not normalized to coherent state value

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:,:)
      real( kind = dp ) :: vout

      complex( kind = dp ), dimension(4,3) :: c_spin
      real( kind = dp ) :: MSx, MSy, MSz, MSxSy, MSxSz, MSySz, VSx, VSy, VSz, sina, &
                           sinb, cosa, cosb, A, B

      c_spin = c1d_SpinMatrix( f )

      MSx = real( c_spin(1,1), dp )
      MSy = real( c_spin(1,2), dp )
      MSz = real( c_spin(1,3), dp )
      VSx = real( c_spin(2,1), dp ) - MSx**2
      VSy = real( c_spin(3,2), dp ) - MSy**2
      VSz = real( c_spin(4,3), dp ) - MSz**2
      MSxSy = real( c_spin(2,2) + c_spin(3,1), dp )
      MSxSz = real( c_spin(2,3) + c_spin(4,1), dp )
      MSySz = real( c_spin(3,3) + c_spin(4,2), dp )

      sina = 0.0_dp; sinb = 0.0_dp; cosa = 1.0_dp; cosb = 1.0_dp
      if (sqrt(MSx**2 + MSy**2) .gt. 1.0d-10) then
         sina = sqrt( (MSx**2 + MSy**2)/(MSx**2 + MSy**2 + MSz**2) )
         sinb = sqrt( MSy**2/(MSx**2 + MSy**2) )
         cosa = sqrt( MSz**2/(MSx**2 + MSy**2 + MSz**2) )
         cosb = sqrt( MSx**2/(MSx**2 + MSy**2) )
      endif

      A = (sinb**2 - cosa**2 * cosb**2) * VSx + (cosb**2 - cosa**2 * sinb**2) * VSy - sina**2 * VSz &
          - (1.0_dp + cosa**2) * cosb * sinb * (MSxSy - 2.0_dp * MSx*MSy) &
          + sina * cosa * cosb * (MSxSz - 2.0_dp * MSx * MSz) &
          + sina * cosa * sinb * (MSySz - 2.0_dp * Msy * MSz)

      B = 2.0_dp * cosa * sinb * cosb * (VSx - VSy) - & 
          cosa * (cosb**2 - sinb**2) * (MSxSy - 2.0_dp * MSx * MSy) - &
          sina * sinb * (MSxSz - 2.0_dp * MSx * MSz) + sina * cosb * (MSySz - 2.0_dp * MSy * MSz)

      vout = 0.5_dp * (cosa**2 * cosb**2 + sinb**2) * VSx + 0.5_dp * (cosa**2 * sinb**2 + cosb**2) * VSy &
             + 0.5_dp * sina**2 * VSz - 0.5_dp * sina**2 * sinb * cosb * (MSxSy - 2.0_dp * MSx * MSy) & 
             - 0.5_dp * sina * cosa * cosb * (MSxSz - 2.0_dp * MSx * MSz) &
             - 0.5_dp * sina * cosa * sinb * (MSySz - 2.0_dp * MSy * MSz) - 0.5_dp * sqrt(A**2 + B**2)
      
      if ( sqrt(MSx**2 + MSy**2 + MSz**2) .eq. 0.0_dp) then
         vout = 0.0_dp
         print *, "Mean spin equal to 0!"
      else
         vout = vout/(MSx**2 + MSy**2 + MSz**2)
      endif

   end function

   function c2d_SpinSqueezing( f ) result vout
      ! Not normalized to coherent state value

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:,:,:)
      real( kind = dp ) :: vout

      complex( kind = dp ), dimension(4,3) :: c_spin
      real( kind = dp ) :: MSx, MSy, MSz, MSxSy, MSxSz, MSySz, VSx, VSy, VSz, sina, &
                           sinb, cosa, cosb, A, B

      c_spin = c2d_SpinMatrix( f )

      MSx = real( c_spin(1,1), dp )
      MSy = real( c_spin(1,2), dp )
      MSz = real( c_spin(1,3), dp )
      VSx = real( c_spin(2,1), dp ) - MSx**2
      VSy = real( c_spin(3,2), dp ) - MSy**2
      VSz = real( c_spin(4,3), dp ) - MSz**2
      MSxSy = real( c_spin(2,2) + c_spin(3,1), dp )
      MSxSz = real( c_spin(2,3) + c_spin(4,1), dp )
      MSySz = real( c_spin(3,3) + c_spin(4,2), dp )

      sina = 0.0_dp; sinb = 0.0_dp; cosa = 1.0_dp; cosb = 1.0_dp
      if (sqrt(MSx**2 + MSy**2) .gt. 1.0d-10) then
         sina = sqrt( (MSx**2 + MSy**2)/(MSx**2 + MSy**2 + MSz**2) )
         sinb = sqrt( MSy**2/(MSx**2 + MSy**2) )
         cosa = sqrt( MSz**2/(MSx**2 + MSy**2 + MSz**2) )
         cosb = sqrt( MSx**2/(MSx**2 + MSy**2) )
      endif

      A = (sinb**2 - cosa**2 * cosb**2) * VSx + (cosb**2 - cosa**2 * sinb**2) * VSy - sina**2 * VSz &
          - (1.0_dp + cosa**2) * cosb * sinb * (MSxSy - 2.0_dp * MSx*MSy) &
          + sina * cosa * cosb * (MSxSz - 2.0_dp * MSx * MSz) &
          + sina * cosa * sinb * (MSySz - 2.0_dp * Msy * MSz)

      B = 2.0_dp * cosa * sinb * cosb * (VSx - VSy) - & 
          cosa * (cosb**2 - sinb**2) * (MSxSy - 2.0_dp * MSx * MSy) - &
          sina * sinb * (MSxSz - 2.0_dp * MSx * MSz) + sina * cosb * (MSySz - 2.0_dp * MSy * MSz)

      vout = 0.5_dp * (cosa**2 * cosb**2 + sinb**2) * VSx + 0.5_dp * (cosa**2 * sinb**2 + cosb**2) * VSy &
             + 0.5_dp * sina**2 * VSz - 0.5_dp * sina**2 * sinb * cosb * (MSxSy - 2.0_dp * MSx * MSy) & 
             - 0.5_dp * sina * cosa * cosb * (MSxSz - 2.0_dp * MSx * MSz) &
             - 0.5_dp * sina * cosa * sinb * (MSySz - 2.0_dp * MSy * MSz) - 0.5_dp * sqrt(A**2 + B**2)
      
      if ( sqrt(MSx**2 + MSy**2 + MSz**2) .eq. 0.0_dp) then
         vout = 0.0_dp
         print *, "Mean spin equal to 0!"
      else
         vout = vout/(MSx**2 + MSy**2 + MSz**2)
      endif

   end function

   function c3d_SpinSqueezing( f ) result vout
      ! Not normalized to coherent state value

      implicit none
      complex( kind = dp ), intent(in), allocatable :: f(:,:,:,:,:)
      real( kind = dp ) :: vout

      complex( kind = dp ), dimension(4,3) :: c_spin
      real( kind = dp ) :: MSx, MSy, MSz, MSxSy, MSxSz, MSySz, VSx, VSy, VSz, sina, &
                           sinb, cosa, cosb, A, B

      c_spin = c3d_SpinMatrix( f )

      MSx = real( c_spin(1,1), dp )
      MSy = real( c_spin(1,2), dp )
      MSz = real( c_spin(1,3), dp )
      VSx = real( c_spin(2,1), dp ) - MSx**2
      VSy = real( c_spin(3,2), dp ) - MSy**2
      VSz = real( c_spin(4,3), dp ) - MSz**2
      MSxSy = real( c_spin(2,2) + c_spin(3,1), dp )
      MSxSz = real( c_spin(2,3) + c_spin(4,1), dp )
      MSySz = real( c_spin(3,3) + c_spin(4,2), dp )

      sina = 0.0_dp; sinb = 0.0_dp; cosa = 1.0_dp; cosb = 1.0_dp
      if (sqrt(MSx**2 + MSy**2) .gt. 1.0d-10) then
         sina = sqrt( (MSx**2 + MSy**2)/(MSx**2 + MSy**2 + MSz**2) )
         sinb = sqrt( MSy**2/(MSx**2 + MSy**2) )
         cosa = sqrt( MSz**2/(MSx**2 + MSy**2 + MSz**2) )
         cosb = sqrt( MSx**2/(MSx**2 + MSy**2) )
      endif

      A = (sinb**2 - cosa**2 * cosb**2) * VSx + (cosb**2 - cosa**2 * sinb**2) * VSy - sina**2 * VSz &
          - (1.0_dp + cosa**2) * cosb * sinb * (MSxSy - 2.0_dp * MSx*MSy) &
          + sina * cosa * cosb * (MSxSz - 2.0_dp * MSx * MSz) &
          + sina * cosa * sinb * (MSySz - 2.0_dp * Msy * MSz)

      B = 2.0_dp * cosa * sinb * cosb * (VSx - VSy) - & 
          cosa * (cosb**2 - sinb**2) * (MSxSy - 2.0_dp * MSx * MSy) - &
          sina * sinb * (MSxSz - 2.0_dp * MSx * MSz) + sina * cosb * (MSySz - 2.0_dp * MSy * MSz)

      vout = 0.5_dp * (cosa**2 * cosb**2 + sinb**2) * VSx + 0.5_dp * (cosa**2 * sinb**2 + cosb**2) * VSy &
             + 0.5_dp * sina**2 * VSz - 0.5_dp * sina**2 * sinb * cosb * (MSxSy - 2.0_dp * MSx * MSy) & 
             - 0.5_dp * sina * cosa * cosb * (MSxSz - 2.0_dp * MSx * MSz) &
             - 0.5_dp * sina * cosa * sinb * (MSySz - 2.0_dp * MSy * MSz) - 0.5_dp * sqrt(A**2 + B**2)
      
      if ( sqrt(MSx**2 + MSy**2 + MSz**2) .eq. 0.0_dp) then
         vout = 0.0_dp
         print *, "Mean spin equal to 0!"
      else
         vout = vout/(MSx**2 + MSy**2 + MSz**2)
      endif

   end function

end module
