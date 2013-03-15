      subroutine tstcvg(x,y,yp,nn,iy,ea,dt,irscvg,icry)
c
c  Test convergence properties, as defined by ea, and if things 
c  look bad, try to reset solution from previous solution and
c  derivatives stored in common /cyfstr/ in s/r rhs.
c
c  Note: for the time being just test on the first 4 variables
c
c  If irscvg ne 0, do suitable resetting of solution.
c
c  In case of no problems, icry is returned as 0
c
c  Original version: 21/10/02
c
      implicit double precision(a-h, o-z)
c
c  include statement setting storage indices for reaction rate
c  and energy generation parameters
c
      include 'engenr.bz.d.incl'
      parameter (ivrmx4 = 4*ivarmx)
      dimension x(*), y(iy,*), yp(iy,*), ea(iy,*)
      dimension dymax(ivarmx), ndymax(ivarmx)
      common/rhcn/ a1,a2,a3,a4,zdummy,nvar,ifwrt,irhtst,inentr,
     *  ii1,ii2,ii3,icomp
      common/cyfstr/ yzfstr(ivrmx4,nnmax)
      common/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng,idghe3,idgmxc,
     *  idgsdt,idgn75,idgn76,idgn77
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
c  hard-coded limit for maximum change in ea:
c
      ealim=2.d0
      iitest=4
c
      icry=0
c
c  find maximum of maximum change
c
      eamax=0.d0
      do i=1,iitest
	if(ea(i,2).gt.eamax) then
	  imax=i
	  eamax=ea(i,2)
        end if
      end do
c
c  test if below limit and, if so, simply return
c
      if(eamax.le.ealim) then
	return
      end if
c
      nmax=ea(imax,3)
      write(istdou,105) imax, eamax, nmax
      if(istdpr.ne.istdou.and.istdpr.gt.0) 
     *  write(istdpr,105) imax, eamax, nmax
c
c  find actual extreme change in all iitest variables and its location
c  as well as actual maximum change
c
      aamax=0.d0
      do i=1, iitest
	dymax(i)=0
	do n=1,nn
	  ddy=y(i,n)-yzfstr(i,n)
	  if(abs(ddy).gt.abs(dymax(i))) then
	    dymax(i)=ddy
	    ndymax(i)=n
          end if
        end do
	if(abs(dymax(i)).gt.aamax) then
	  aamax=abs(dymax(i))
	  iamax=i
        end if
      end do
c
      write(istdou,110) (i, dymax(i),ndymax(i),x(ndymax(i)),i=1,iitest)
      if(istdpr.ne.istdou.and.istdpr.gt.0) 
     *  write(istdpr,110) (i, dymax(i),ndymax(i),x(ndymax(i)),
     *    i=1,iitest)
c
      if(irscvg.eq.0.or.aamax.lt.ealim) return
c
c  reset solution for a suitable variable, to match derivative,
c  as a rather desparate measure
c
      if(iamax.eq.4) then
        ir=3
        ifr=2*nvar+ir
      else
        write(istdou,115) iamax
        if(istdou.ne.istdpr.and.istdpr.gt.0)
     *    write(istdpr,115) imax
	stop 'tstcvg'
      end if
c
      do i=1,iitest
        do n=1,nn
	  y(i,n)=yzfstr(i,n)
        end do
      end do
c
      ifr=2*nvar+3
      y(3,nn)=yzfstr(3,nn)
      do n=nn-1,1,-1
	y(3,n)=y(3,n+1)
     *         -0.5d0*(x(n+1)-x(n))*(yzfstr(ifr,n+1)+yzfstr(ifr,n))
      end do
c
      ifr=2*nvar+4
      iar=3*nvar+2
      do n=2,nn
        y(2,n)=yp(2,n)+((yzfstr(4,n)-yzfstr(4,n-1))/(x(n)-x(n-1))
     *         -0.5d0*(yzfstr(ifr,n)+yzfstr(ifr,n-1)))*dt/yzfstr(iar,n)
      end do
c
      write(istdou,120) ir
      if(istdou.ne.istdpr.and.istdpr.gt.0)
     *  write(istdpr,120) ir
      icry=1
      return
  105 format(/' ***** Warning in s/r tstcvg: for variable',i3,
     *  ' maximum change =',1pe13.5,' at n =',i5)
  115 format(/' For imax = ', i2, ' no reset defined in tstcvg')
  110 format(/' Variable, extreme change, location of extreme change:'/
     *   (i3, 1pe13.5,i5,e13.5))
  120 format(/' Solution for variable ',i2,' has been reset')
      end
