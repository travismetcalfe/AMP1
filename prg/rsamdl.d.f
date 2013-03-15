      subroutine rsamdl(x,aa,data,n1,n2,n3,nn,ia,w,iw)
c
c  resets adiabatic model.
c  **********************
c
c  when n3 .gt. 0 and .lt. nn sets rho derivative numerically for
c  n .ge. n3.
c
c  when n1, n2 .gt. 0 and n1. lt. n2 resets rho and derivative for
c  n between n1 and n2 such that both rho ann derivative are
c  continuous at n1 and n2.
c
c  Modified 16/3/89, passing work array w through argument list.
c     Note that dimension iw must be at least 5.
c
c  Double precision version.
c  ++++++++++++++++++++++++
c
c  Dated 11/3/90
c
      implicit double precision (a-h,o-z)
      dimension x(1),aa(ia,1),data(1),w(iw,1)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      save
c
      if(istdpr.gt.0) write(istdpr,100) n1,n2,n3,nn
c
c  set rho hat and its derivative
c
      do 10 n=2,nn
      w(1,n)=aa(1,n)*aa(5,n)
   10 w(2,n)=-(aa(2,n)+aa(4,n))*w(1,n)/x(n)
      w(1,1)=aa(1,1)*aa(5,1)
      w(2,1)=0
c
c  test for resetting of derivative by numerical differentiation
c
      if(n3.le.0.or.n3.ge.nn) go to 30
      nn3=nn-n3+1
      call derive(x(n3),w(1,n3),w(4,n3),nn3,iw,iw,1,1)
c  diagnostic output
      if(istdpr.gt.0) write(istdpr,110) 
     *  (n,x(n),w(1,n),w(2,n),w(4,n),n=n3,nn)
      do 20 n=n3,nn
   20 w(2,n)=w(4,n)
c
c  reset rho hat between n1 and n2
c
   30 if(n1.le.0.or.n2.le.0.or.n1.ge.n2) go to 50
      init=1
      do 40 n=n1,n2
      w(3,n)=cubfit(x(n),x(n1),x(n2),w(1,n1),w(1,n2),w(2,n1),
     *  w(2,n2),w(4,n),init)
   40 init=0
c
      if(istdpr.gt.0) write(istdpr,120) (n,x(n),(w(i,n),i=1,4),n=n1,n2)
c  reset vg and u by scaling to new rho
      do 45 n=n1,n2
      w(2,n)=w(4,n)
      w(4,n)=w(3,n)/w(1,n)
      w(1,n)=w(3,n)
      aa(2,n)=w(4,n)*aa(2,n)
   45 aa(5,n)=w(4,n)*aa(5,n)
c
c  reset aa(4,n)
c
   50 nr1=n1
      if(n3.gt.0.and.n3.lt.n1) nr1=n3
      nr2=n2
      if(nn3.gt.0.and.n3.lt.nn) nr2=nn
      do 60 n=nr1,nr2
   60 aa(4,n)=-x(n)*w(2,n)/w(1,n)-aa(2,n)
      if(istdpr.gt.0) write(istdpr,130) 
     *  (n,x(n),(aa(i,n),i=1,5),n=nr1,nr2)
      return
  100 format(///' reset adiabatic model. n1, n2, n3, nn =',
     *  4i5)
  110 format(//' numerical differentiation of rho. n, x, rho hat,',
     *  ' old derivative, new derivative:'//(i4,0pf12.6,1p3e15.7))
  120 format(//' continuous resetting of rho hat derivative.',
     *  ' n, x, old rho hat and derivative, new rho hat and',
     *  ' derivative:'//(i4,0pf12.6,1p4e15.7))
  130 format(//' reset model. n, x, aa(1-5):'//
     *  (i4,0pf12.6,1p5e15.7))
      end
      double precision function cubfit(x,x1,x2,y1,y2,dy1,dy2,dy,init)
c
c  calculates third degree polynomia  at x, with values and
c  derivatives yi and dyi at xi, i = 1,2.
c
c  derivative of polynomial is returned in dy.
c
c  init must be 1 in first call and can be reset to .ne. 1
c  for subsequent calls.
c
c  if x1 = x2 a diagnostics is written, cubfit and dy returned as
c  zero and init reset to -1
c
c  Double precision version.
c  ++++++++++++++++++++++++
c
c  Dated 11/3/90
c
      implicit double precision (a-h,o-z)
      dimension c(4)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      dx=x2-x1
      if(dx.eq.0) go to 90
c
      if(init.ne.1) go to 20
c
      c(1)=y1
      c(2)=dy1*dx
      c(3)=3*(y2-y1)-(2*dy1+dy2)*dx
      c(4)=2*(y1-y2)+(dy1+dy2)*dx
c
   20 t=(x-x1)/dx
      cubfit=c(1)+t*(c(2)+t*(c(3)+t*c(4)))
      dy=(c(2)+t*(2*c(3)+3*t*c(4)))/dx
      return
c  diagnostics for zero range
   90 if(istdpr.gt.0) write(istdpr,100) x1,x2
      cubfit=0
      dy=0
      init=-1
      return
  100 format(//1x,10(1h*),' zero range in cubfit. x1, x2 =',
     *  1p2e15.7)
      end
