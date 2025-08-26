c     Generalized inverse Courrieu (2005)
c     The GINV of mxn X is X'LMML'G', with m>n
c     where M = inv(L'L) and L'L is full rank rxr
c     The subroutine only returns L and r
c     when m=n, L is non singular and the solution is
c     inv(L')inv(L)G'

      subroutine pinv(x, m, n, r, l)
      integer n, m, k, r
      double precision x(m,n), ix(n,m), a(n,n), tol
      double precision l(n,n), tmp(n,1)

      a = matmul(transpose(x), x)
      tol = 0.0d0
      do i=1,n
         if (a(i,i) .gt. 0) then
            if (tol .eq. 0.0d0) then
               tol = a(i,i)
            else
               if (a(i,i) .lt. tol) then
                  tol = a(i,i)
               end if
            endif
         end if
      end do
      tol = tol*1.0d-9
      l = 0.0d0
      r = 0
      do k=1,n
         r = r+1
         if (r .eq. 1) then
            l(k:n,r) = a(k:n,k)
         else
            l(k:n,r) = a(k:n,k)-matmul(l(k:n,1:(r-1)), l(k,1:(r-1)))
         end if
         if (l(k,r) .gt. tol) then
            l(k,r) = sqrt(l(k,r))
            if (k .lt. n) then
               l((k+1):n,r) = l((k+1):n,r)/l(k,r)
            end if
         else
            r = r-1
         end if
      end do
      end
      
