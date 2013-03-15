c...5...10...14...15....

      subroutine svd(m, n, a, v, w, xxx)

cccc  input, output
      implicit none
      integer n,m,i, k, j
      double precision a(m,n), v(n,n),w(n)
cccc  internal
c      parameter (NMAX = 500)
      integer l, its, jj,  nm
      double precision anorm, g, scale, rv1(500), c, f, h, s, x, y, 
     +     z, pythag, xxx(m,n)

c      common/svd_par/ n, m, a(m,n), w(n), v(n,n), xxx(m,n)


cccc  save A
      do i=1, n
         do j=1,m 
            xxx(j,i) = a(j,i)
         enddo
      enddo

cccc  initialize
      g = 0.0
      scale = 0.0
      anorm = 0.0
cccc  1st of four loops

      do i = 1, n
         l = i + 1
         rv1(i) = scale * g
         g = 0.0
         s = 0.0 
         scale = 0.0
         if (i.le.m) then
            do k = i, m
               scale = scale + abs(a(k,i))
            enddo
            if (scale.ne.0) then
               do k = i,m
                  a(k,i) = a(k,i) / scale              
                  s = s + a(k,i) * a (k,i)
               enddo
               f = a(i,i)
               g = -sign(sqrt(s),f)
               h = f*g - s
               a(i,i) = f-g
               do j = l, n
                  s = 0.0
                  do k = i,m
                     s = s +  a(k,i) * a(k,j)
                  enddo
                  f = s/h
                  do k = i, m
                     a(k,j) = a(k,j) + f*a(k,i)
                  enddo
               enddo
               do k = i,m
                  a(k,i) = scale * a(k,i)
               enddo
            endif
         endif     

         w(i) = scale * g
         g = 0.0
         s = 0.0
         scale = 0.0
         if( (i.le.m) .and. (i.ne.n)) then 
            do k = l, n
               scale = scale + abs(a(i,k))
            enddo
            if (scale.ne.0.0) then 
               do k = l , n
                  a(i,k) = a(i,k) /scale
                  s = s + a(i,k) * a(i,k)
               enddo
               f = a(i,l)
               g = -sign(sqrt(s),f)
               h = f*g - s
               a(i,l) = f - g
               do k = l, n
                  rv1(k) = a(i,k)/h
               enddo
               do j = l, m
                  s = 0.0
                  do k = l, n
                     s = s + a(j,k ) * a(i,k)
                  enddo
                  do k= l,n
                     a(j,k) = a(j,k) + s* rv1(k)
                  enddo
               enddo
               do k = l, n
                  a(i,k) = scale*a(i,k)
                enddo
            endif
         endif
         anorm = max ( anorm, (abs(w(i)) + abs(rv1(i))))
      enddo
cccc  2nd of four loops
      do i = n, 1, -1
         if (i.lt.n) then
            if(g.ne.0.0) then
               do j = l,n
                  v(j,i) = (a(i,j)/a(i,l))/g
               enddo
               do j = l, n
                  s = 0.0
                  do k=l,n
                     s = s + a(i,k)* v(k,j)
                  enddo
                  do k = l, n
                     v(k,j) = v(k,j) + s*v(k,i)
                  enddo
               enddo
            endif
            do j = l,n 
               v (i,j) = 0.0
               v(j,i) = 0.0
            enddo
         endif
         v(i,i) = 1.0
         g = rv1(i)
         l=i
      enddo
         
cccc 3rd of four loops
      do i = min(m,n), 1, -1 
         l = i+1
         g = w(i)
         do j = l, n
            a(i,j) = dble(0.0)
         enddo
         if (g.ne.0.0) then
            g = dble(1.0)/g
            do j = l, n
               s = dble(0.0)
               do k = l,m
                  s = s + a(k,i) * a(k,j)
               enddo
               f = (s/a(i,i)) * g
               do k = i, m
                  a(k,j) = a(k,j) + f*a(k,i)
               enddo
            enddo
            do j = i, m
               a(j,i) = a(j,i) * g
            enddo
         else
            do j = i, m
               a(j,i) = dble(0.0)
            enddo
         endif   
         a(i,i) = a(i,i) + dble(1.0)
      enddo

cccc 4th of four loops
      do k = n, 1, -1 
         do its = 1, 30 
            do l = k, 1, -1 
               nm = l - 1
               if ((abs(rv1(l)) + anorm).eq.anorm) goto 2
               if ((abs(w(nm)) + anorm).eq.anorm) goto 1
            enddo
 1          c = dble(0.0)
            s = dble(1.0)
            do i = l,k 
               f = s*rv1(i)
               rv1(i) = c * rv1(i)
               if ((abs(f) + anorm).eq. anorm) goto 2
               g = w(i)
               h = pythag(f,g)
               w(i) = h
               h = dble(1.0)/h
               c = (g*h)
               s = -(f*h)
               do j = 1, m
                  y = a(j, nm)
                  z = a(j,i)
                  a(j, nm) = (y *c) + (z*s)
                  a(j,i) = - (y*s) + (z*c)
              enddo
            enddo
 2          z = w(k)
            if(l.eq.k) then 
               if (z.lt.0.0) then
                  w(k) = -z
                  do j = 1, n
                     v(j,k) = -v(j,k)
                  enddo
               endif
               goto 3
            endif
            if (its.eq.30) pause 'no convergence in svd'
            x = w(l)
            nm = k - 1
            y = w(nm)
            g = rv1(nm)
            h = rv1(k)
            f = ((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y)
            g = pythag(f,dble(1.0))
            f = ((x-z)*(x+z)+h*((y/(f+sign(g,f))) - h)) / x
            c = dble(1.0)
            s = dble(1.0)
            do j = l, nm
               i = j+1
               g = rv1(i)
               y = w(i)
               h = s * g
               g = c * g
               z = pythag(f,h)
               rv1(j) = z
               c = f/z
               s = h/z
               f = (x*c) + (g*s)
               g = -(x*s) + (g * c)
               h = y*s
               y = y*c
               do jj = 1, n
                  x = v(jj,j)
                  z = v(jj,i)
                  v(jj,j) = (x*c) + (z*s)
                  v(jj,i) = - (x*s) + (z*c)
               enddo
               z = pythag(f,h)
               w(j) = z
               if (z.ne.0.0) then
                  z = dble(1.0)/z
                  c = f*z
                  s = h*z
               endif
               f = (c*g) + (s*y)
               x = -(s*g) + (c*y)
               do jj = 1, m
                  y = a(jj,j)
                  z = a(jj,i)
                  a(jj,j) = (y*c) + (z*s)
                  a(jj,i) = - (y*s) + (z*c)
c     here
               enddo
            enddo
            rv1(l) = 0.0
            rv1(k) = f
            w(k) = x
         enddo
 3       continue
      enddo


c     
c     transfer to old physical dimensions for return to main program
c


c      do j = 1, n
c         ww(j) = w(j)
c         do i = 1, m
c            xx(i, j) = xxx (i,j)
c            aa(i, j) =  a(i, j)
c            vv(i, j) = v (i,j)
c         enddo
c      enddo

      return
      end




      function pythag (a, b) 
      
      double precision a, b, pythag
      double precision absa, absb
      
      absa = abs(a)
      absb = abs(b)

      if (absa.gt.absb) then 
         pythag = absa * sqrt (1. + (absb/absa)**2)
      else
         if (absb.eq.0) then
            pythag = 0.
         else
            pythag = absb * sqrt (1. + (absa/absb)**2)
         endif
      endif

      return
      end




      subroutine singvd (ma, nobs, a, b, ax, x, aix, invs)


c     a = the matrix (replaced by inverse)
c     b = the rhs (da) 
c     if inv = 1 then invert, not solve
c     ax = original matrix a
c     aix = inverse??????
c     returns x = original matrix matrix solution or the inverse
c     x = the solution

      implicit none
      integer ma, nobs, invs
      double precision a(nobs, ma), b(nobs), x(ma), ax(nobs,ma)
      
      double precision v(ma, ma), w(ma), s, tmp(ma), tu(ma, nobs)
      double precision u(ma, nobs), aix(ma, nobs)
      integer j, i, jj

c     save a
      do i= 1, ma
         do j = 1, nobs
            ax(j,i) = a(j,i)
         enddo
      enddo

c     call svd to get the things... a is replaced by u
      call svd(nobs, ma, a, v, w, aix)

c     change low values to 0
      do j = 1, ma 
         if (w(j).lt.1) then 
            w(j) = 0
         endif 
      enddo
 
c     solve system x = (v##(ww##(transpose(u)##b)))
      do j = 1, ma
         s = 0.
         if (w(j).ne.0) then
            do i = 1, nobs
               s = s+a(i,j)*b(i)
            enddo
            s = s/w(j)
         endif
         tmp(j) = s
      enddo
      do j =1, ma
         s = 0.
         do jj = 1, ma
            s = s + v(j, jj)*tmp(jj)
         enddo
         x(j) = s
      enddo

      if (invs.eq.1) then 
c     calculate transpose(u)

         do i = 1, ma
            do j = 1, nobs 
               tu(i,j) = a(j,i)
            enddo
         enddo



c     now calculate w##transpose(u)  = aw

         do i = 1, ma 
            do j = 1, nobs
               u (i, j) =  tu( i,j)/w(i)
            enddo
         enddo

c     now calculate v##aw
         do i = 1, nobs
            do j = 1, ma
               s= 0.
               do jj = 1, ma 
                  s = s + v(j,jj) * u(jj,i)
               enddo
               aix(j,i) = s
            enddo
         enddo

       endif

       
      end

