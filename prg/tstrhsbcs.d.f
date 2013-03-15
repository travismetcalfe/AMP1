      subroutine tstrhs(rhs,x,y,zk,ap,aq,ii,ii1,kk,nn,iy,epsd,nstp)
c
c  test right hand subroutine for tnrkt
c
c  test derivatives in dzdy and fd and, as of 30/12/1984,
c  also that solution satisfies equations.
c
c  note that time derivatives are not included. thus testing
c  is only effective for ii2 = ii3 = 0.
c
c  Double precision version.
c  ++++++++++++++++++++++++
c
c  Dated 11/3/90
c
      implicit double precision (a-h,o-z)
      dimension x(1),y(iy,1),zk(1),ap(1),aq(1)
      dimension z(10,3),dzdy(10,10,3),alam(1),alamd(1),h(1),
     *  hd(1),y1(10),f(10,3),fd(10,10,3)
      common/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      external rhs
c
      save
c
c  step through model
c
      do 70 n=1,nn,nstp
      if(istdpr.gt.0) write(istdpr,100) n,x(n),epsd,(y(i,n),i=1,ii)
c
      do 50 i=1,ii
c
      call store(y(1,n),y1,ii1)
      y1(i)=y1(i)-epsd
      do 20 k=1,3
      call rhs(x(n),y1,zk,z(1,k),dzdy(1,1,k),ap,aq,f(1,k),
     *  fd(1,1,k),lm,almd,h,hd,10,10,1,1,n,1)
   20 y1(i)=y1(i)+epsd
c  set and print derivatives
      if(istdpr.gt.0) write(istdpr,110) i
      do 25 j=1,ii
      zd1=(z(j,3)-z(j,1))/(2*epsd)
      dzd=(zd1-dzdy(j,i,2))/(abs(dzdy(j,i,2))+1.e-10)
      if(istdpr.gt.0) write(istdpr,120) j,zd1,dzdy(j,i,2),dzd
   25 continue
      if(istdpr.gt.0) write(istdpr,130) i
      do 30 j=1,ii
      fd1=(f(j,3)-f(j,1))/(2*epsd)
      dfd=(fd1-fd(j,i,2))/(abs(fd(j,i,2))+1.e-10)
      if(istdpr.gt.0) write(istdpr,120) j,fd1,fd(j,i,2),dfd
   30 continue
   50 continue
c
c  for n .lt. nn, test that equations are satisfied
c
   52 if(n.eq.nn) go to 70
c
      do 55 k=1,2
      n1=n+k-1
   55 call rhs(x(n1),y(1,n1),zk,z(1,k),dzdy(1,1,k),ap,aq,f(1,k),
     *  fd(1,1,k),lm,almd,h,hd,10,10,1,1,n1,1)
c
      dx=x(n+1)-x(n)
      if(istdpr.gt.0) write(istdpr,140) n,n1,dx
c
      do 60 i=1,ii
      dz1=(z(i,2)-z(i,1))/dx
      dz2=0.5*(f(i,1)+f(i,2))
      ddz=(dz1-dz2)/(abs(dz2)+1.e-30)
      if(istdpr.gt.0) write(istdpr,145) 
     *  i,(z(i,k),k=1,2),(f(i,k),k=1,2),dz1,dz2,ddz
   60 continue
   70 continue
      return
  100 format(///' test s/r rhs at n =',i3,'   x =',1pe13.5,
     *  '  with epsd =',e13.5/
     *  '  y:', 10e13.5)
  110 format(//' dzdy(j,',i2,'):'/)
  120 format(i4,1p2e13.5,5x,e11.3)
  130 format(//' fd(j,',i2,'):'/)
  140 format(/' test that equations are satisfied, between n =',
     *  i5,'  and n =',i5,' .  dx =',1pe15.7/
     *  ' i, z(i,1), z(i,2), f(i,1), f(i,2), dzi/dx, f aver., error'/)
  145 format(i4,1p2e15.7,3x,2e15.7,3x,2e15.7,e11.3)
      end
      subroutine tstbcs(bcs,x1,x2,y1,y2,zk,ap,aq,ii,ii1,kk,ka,kb,
     *  nn,epsd)
c  test boundary condition subroutine for nrk
      implicit double precision (a-h,o-z)
      dimension y1(1),y2(1),zk(1),ap(1),aq(1)
      dimension g(10,3),gd(10,10,3),ys1(10),ys2(10)
      common/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng
      common/noiter/ iter
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      external bcs
c
      save
c
      if(istdpr.gt.0) write(istdpr,100) nn,x1,x2,epsd
c
      iterp=iter
      idgbcp=idgbcs
      iter=1
c
      do 50 i=1,ii
      call store(y1,ys1,ii1)
      call store(y2,ys2,ii1)
      ys1(i)=ys1(i)-epsd
      ys2(i)=ys2(i)-epsd
      do 20 k=1,3
      idgbcs=0
      if(i.eq.1.and.k.eq.2) idgbcs=1
      call bcs(x1,x2,ys1,ys2,zk,ap,aq,g(1,k),gd(1,1,k),10,10,nn)
      ys1(i)=ys1(i)+epsd
   20 ys2(i)=ys2(i)+epsd
c  set and print derivatives
      if(istdpr.gt.0) write(istdpr,110) i
      do 25 j=1,ii
      gd1=(g(j,3)-g(j,1))/(2*epsd)
      dgd=(gd1-gd(j,i,2))/(abs(gd(j,i,2))+1.e-10)
      if(istdpr.gt.0) write(istdpr,120) j,gd1,gd(j,i,2),dgd
   25 continue
   50 continue
c
      iter=iterp
      idgbcs=idgbcp
c
      return
  100 format(///' test s/r bcs with nn =',i3,'   x1 =',1pe13.5,
     *  '  x2 =',e13.5,'  with epsd =',e13.5)
  110 format(//' gd(j,',i2,'):'/)
  120 format(i4,1p2e13.5,5x,e11.3)
      end
