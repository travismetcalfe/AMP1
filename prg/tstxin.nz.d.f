      subroutine tstxin(x,y,nn,iy,icase,xhint,xhintp,xtime,xtimep,icry,
     *   label)
c
c  Calculate integrated hydrogen abundance and check that
c  its variation is sensible (in a manner to be defined)
c
c  For icase = 1: assume output storage of y
c  For icase = 2: assume internal (tnrkt) storage of y in calculation
c  with diffusion
c
c  xhint returns new integrated value. This is assumed to correspond
c  to time xtime (in seconds) and is compared with the input value
c  xhintp at time xtimep.
c
c  label is a string which might be used for location information
c
c  If a possible inconsistency is detected, icry is set to 1, otherwise
c  icry is returned as 0.
c
c  Original version: 5/12/00
c
      implicit double precision (a-h, o-z)
      character*(*) label
c
      include 'engenr.nnz.d.incl'
      dimension x(1), y(iy,1)
      dimension xh(nnmax),xin(nnmax),q(nnmax)
      common /fconst/ pi, ebase, amu, ame, clight, planck, boltzm, 
     *  cgrav, amsun, echar, ergev, syear
      common/totmss/ am, rs
      common/cmpstr/ ispxx3, ispcno, icvxx3, icvcno, ibcxx3, ibccno,
     *  idcomp, iccomp
      common/cmxcor/ qmxcor,rmxcor,frmxc,nmxcor,nmxcp,qmxcp,rmxcp,
     *  qmxmin,qmxmax
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      icry=0
c
      idiag=1
      if(istdpr.le.0) idiag=0
      epslm1=1.d-5
      epslm2=0.1
c
      if(icase.eq.1) then
	ixh=5
	ilum=4
      else
	ixh=4
	ilum=4+idcomp
      end if
c
      do 20 n=1,nn
      xh(n)=y(ixh,n)
   20 q(n)=10.d0**x(n)
c
      q(nn+1)=0
      xh(nn+1)=xh(nn)
c
      call vinta(q,xh,xin,nn+1,1,1)
c
      xhint=-xin(nn+1)
c
      if(idiag.gt.0) write(istdpr,120) label,xtime,xhint
c
c  test for possible inconsistencies
c
      if(xtime.eq.xtimep) then
	if(abs(xhint-xhintp).gt.epslm1.and.istdpr.gt.0) then
	  write(istdpr,130) label,xtime,xhint,xhintp
	  icry=1
        end if
      else
c
	dxhidt = (xhint - xhintp)/(xtime-xtimep)
c
c  compare with expected change, based on luminosity
c
	dxhadt = -1.66d-19*(10.d0**(33+y(ilum,1)))/(am*amsun)
	if(abs(dxhidt/dxhadt-1).gt.epslm2.and.istdpr.gt.0) then
	  write(istdpr,140) label,xtime,xhint,xtimep,xhintp,
     *      dxhidt,dxhadt
	  icry=1
        end if
      end if
      return
  120 format(/' In tstxin, label ',a,' xtime, xhint =',1p2e13.5)
  130 format(/' ***** Warning in tstxin, at label ',a/
     *  ' At xtime =',1pe13.5,' xhint, xhintp  =',2e13.5,
     *  ' differ')
  140 format(/' ***** Warning in tstxin, at label ',a/
     *  ' At xtime =',1pe13.5,' inconsistent change.'/
     *  ' xhint, xhintp, xtimep  =',3e13.5/
     *  ' dXint/dt = ',e13.5,' dXanl/dt =',e13.5)
      end
