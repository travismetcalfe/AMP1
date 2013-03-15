      subroutine stdifv(x,y,am,idiffus,nn,ii3,iy)
c
c  sets diffusion velocity variable into proper place in y,
c  assuming (for the moment, at least) uniform composition
c  On input, am the mass in solar units
c
c  Original version 16/8/92
c
c  Modified 5/8/95, to allow diffusion of several species.
c
c  Modified 20/5/95, moving zh from argument list to common/heavy/
c
c  Corrected 5/8/97, by including engenr.nnz.d.incl, to set nspdmx etc
c
      implicit double precision (a-h, o-z)
c
c  include statement setting storage indices for reaction rate
c  and energy generation parameters
c
      include 'engenr.nnz.d.incl'
      logical noder 
      dimension x(1), y(iy,1)
      dimension grad(6),clamxy(4),velx(nspdmx,2,6),difx(nspdmx,6),
     *  xdif(nspdmx)
      common/eqstd/ xi(4),ane(10),rho(20),ht(20),pt(20),cp(4),dad(4),
     .  dlt(4),gm1,tprh,trhp,rhxp
      common/cmpstr/ ispxx3, ispcno, icvxx3, icvcno, ibcxx3, ibccno,
     *  isphec, idcomp, iccomp
      common/heavy/ zatmos, zhc, zh(1)
c
c  common giving fundamental constants
c
      common /fconst/ pi, ebase, amu, ame, clight, planck, boltzm, 
     *  cgrav, amsun, echar, ergev, syear
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      if(istdpr.gt.0) write(istdpr,*) 'Enter stdifv. ii3 =',ii3
      noder=.true.
      amms=amsun*am
c
c..      write(6,*) 'n, x, velx, flux:'
c
      icomp=idcomp+iccomp
c
      do 30 n=1,nn
      drad=fdrad(x(n),y(1,n),zh(n),ak,akr,akt,akx)
      grad(1)=min(dad(1),drad)
      t=10.d0**y(3,n)
      xh=y(5,n)
      r=1.d11*10.d0**y(1,n)
      amass=amms*10.d0**x(n)
      itbdf1=0
      call difcff(t,rho,pt,xh,zh(n),r,amass,grad,
     *  idiffus,itbdf1,clamxy,velx,difx,nspdmx,noder)
c..c
c..c  zero velocity of Z for testing
c..c
c..      do 20 i=1,2
c..      do 20 j=1,6
c..   20 velx(2,i,j)=0
c
c  set abundances of diffusing species
c
      xdif(1)=xh
      if(idiffus.eq.2) xdif(2)=y(6+iccomp,n)
c
      do 25 j=1,idcomp
      if(n.eq.1.or.n.eq.nn) velx(j,1,1)=0
      if(j.eq.1) then 
        y(4+icomp+j,n)=(velx(j,1,1)/amms)*xdif(j)
      else
        y(4+icomp+j,n)=
     *    ((velx(j,1,1)/amms)+velx(j,2,1)*y(5+icomp,n))*xdif(j)
      end if
   25 continue
   30 continue
c
      return
      end
