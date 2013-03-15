      subroutine opalkz(rl, tl, xh, zh, ak, akr, akt, akx, ier)
c
c  driver routine for opacity tables from OPAL and Kurucz.
c  Calling G. Houdek routines, version 8
c
c  This is a version of opalku, with zh = Z passed in argument list
c  rather than in common/heavy/, for consistency with version
c  of evolution code allowing varying Z.
c
c  Original version: 19/7/95
c
c  Modified 28/11/95 to use GH package version 8.
c  Type of interpolation is defined by value of iwdopc in call to
c  s/r opingh below:
c  iwdopc = iwdop0 + 10*iwdop1, where
c    iwdop0 = 9: use minimum-norm interpolation (as in previous versions)
c    iwdop0 = 8: use birational splines
c
c    iwdop1 = 0: Do not include electron conduction
c    iwdop1 = 1: Include electron conduction
c
      implicit double precision (a-h, o-z)
      character*80 copcvr
      common/opcver/ copcvr
      common/copcgh/ imodgh
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      data copcvr
     *  /' GH interpolation in OPAL and Kurucz tables, version v8'/
      data epsder /1.e-4/
      data zhp /-1.d0/
c
c  set up to call OPAL routine, depending on value of imodgh
c
      rlg=rl-3*(tl-6)
      imodga=iabs(imodgh)
      if(imodga.eq.0.or.imodga.eq.1) then
        call opintc(xh,zh,tl,rlg,ak,akr,akt,akx,akz,iexp,ier)
      else if(imodga.eq.2) then
        call opints(xh,zh,tl,rlg,ak,akr,akt,akx,akz,iexp,ier)
      else
	write(istdou,100) imodgh
	if(istdpr.gt.0.and.istdpr.ne.istdou) write(istdpr,100) imodgh
	stop 'opalkz'
      end if
c
      return
  100 format(//' ***** Error in s/r opalkz: imodgh =',i5,
     *  ' not allowed.')
      end
