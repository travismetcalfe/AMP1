      subroutine setrot(x,y,yp,iy,nn,velrot,init)
c
c  set or reset moment of inertia and angular velocity
c
c  If init = -1 store velrot in omgrot and set initial moment of
c  inertia into riner0 and riner
c
c  If init = 1 set initial moment of inertia into riner0 and riner,
c  without changing omgrot
c
c  If init = 0, compute new moment of inertia and set into riner.
c  Also reset omgrot (for now assumed uniform) to keep constant
c  angular momentum.
c
c  Note: angular velocity is given in sec**(-1), and unit of moment 
c  of inertia is 1.e55 g cm**2.
c
c  Original version: 9/4/01
c
      implicit double precision (a-h, o-z)
c
c  include statement setting storage indices for reaction rate
c  and energy generation parameters
c
      include 'engenr.nnz.d.incl'
      parameter(iw=2)
      dimension x(1), y(iy,1),yp(iy,1)
      dimension q(nnmax)
c
c  common giving fundamental constants
c
      common /fconst/ pi, ebase, amu, ame, clight, planck, boltzm, 
     *  cgrav, amsun, echar, ergev, syear
      common/comgrt/ isprot, irotcn, a5, riner0, riner, omgrt0,
     *  omgrot(1)
      common/totmss/ am, rs
      common/work/ w(iw,1)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
c  test for storing initial omega, from velocity or angular velocity 
c  in velrot
c
      if(init.eq.-1) then
	if(isprot.gt.0) then
	  if(velrot.le.1) then
	    omgrts=velrot
          else
	    omgrts=velrot/(1.d11*10.d0**y(1,1))
          end if
	  write(istdou,110) omgrts
	  if(istdpr.ne.istdou.and.istdpr.gt.0) write(istdpr,110) omgrts
        else
	  omgrts=0
        end if
c
	omgrt0=omgrts
        do 20 n=1,nn
   20   omgrot(n)=omgrts
c
      end if
c
c  calculate moment of inertia for model
c
      do 30 n=1,nn
      q(n)=10.d0**x(n)
   30 w(1,n)=10.d0**(2*y(1,n))
c
      call vinta(q,w(1,1),w(2,1),nn,iw,iw)
c
c..      write(6,*) 'Moment of inertia integration'
c..      write(6,'(i5,1p3e13.5)') (n,q(n),w(1,n),w(2,n),n=1,nn,20)
c..      write(6,*) 'amsun, am, velrot, w(2,nn)',amsun,am, velrot,
c..     *  w(2,nn)
      riner=-2.d0*(1.d-33*amsun)*am*w(2,nn)/3.d0
      if(istdpr.gt.0) write(istdpr,120) riner
c
      if(init.eq.1.or.init.eq.-1) then
	riner0=riner
      else
c
c  set new angular velocity
c
	omgrts=omgrt0*riner0/riner
        do 40 n=1,nn
   40   omgrot(n)=omgrts
	if(istdpr.gt.0) write(istdpr,110) omgrot(1)
c
      end if
c
      return
c
  110 format(//' Surface angular velocity set to ',1pe13.5,' sec**(-1)')
  120 format(//' Moment of inertia =',1pe13.5,' 1.e55 g cm**2')
      end
