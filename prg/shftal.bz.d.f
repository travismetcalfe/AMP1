      subroutine shftal(x,y,ynew,nn,iy,init)
c
c  Transforms luminosity, including a shift in the luminosity variable 
c  used for computing when the minimum average energy generation rate
c  is less than the value at the surface by a factor hardcoded to 0.05
c  In this case set alshft such that the value is equal to this minimum
c  Otherwise reset shift to zero
c
c  Resetting only when init = 1.
c
c  Modified 23/6/03, changing hardcoded limit from 0.05 to 0.20.
c
      implicit double precision(a-h, o-z)
      dimension x(1), y(iy,1), ynew(iy,1)
c
      common/clshft/ alshft
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      data eprmin /0.20d0/
c
c  test for setting luminosity shift
c
      alshfp=alshft
      als=10.d0**y(4,1)-alshfp
      if(init.eq.1) then
	do n=1,nn
	  dll=eprmin*10.d0**x(n)-(10.d0**y(4,n)-alshfp)/als
	  if(n.eq.1) then
	    dllmax=dll
	    nmax=1
	  else if(dll.gt.dllmax) then
	    dllmax=dll
	    nmax=n
	  end if
        end do
	if(istdpr.gt.0) write(istdpr,110) dllmax,nmax
	if(dllmax.gt.0) then
	  alshft=dllmax*als
	else
	  alshft=0.d0
        end if
	if(istdpr.gt.0) write(istdpr,120) alshft
      end if
c
c  test for resetting y(4,.)
c
      if(alshft.ne.alshfp) then
	do 20 n=1,nn
        y(4,n)=log10(10.d0**y(4,n)+alshft-alshfp)
        ynew(4,n)=log10(10.d0**ynew(4,n)+alshft-alshfp)
   20   continue
c
      end if
      return
  110 format(/' In shftal, maximum (epsav_lim*q - L/Ls) =',f10.5,
     *   ' at n =',i5)
  120 format(//' In shftal, alshft reset to ',1pe13.5/)
      end
