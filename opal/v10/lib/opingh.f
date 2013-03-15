      subroutine opingh(inopc,iwdopc)
c
c  initialize G. Houdek opacity routines
c
c  For consistency, unit inopc is assumed to have been set to
c  file name containing list of opacity tables, with s/r ofiles
c
c  Choice of procedure is set by iwdopc = iwdop0 + 10*iwdop1 + 100*iwdop2
c    iwdop0 = 9: use minimum-norm interpolation (as in previous versions)
c    iwdop0 = 8: use birational splines
c
c    iwdop1 = 0: Do not include electron conduction
c    iwdop1 = 1: Include electron conduction
c
c    iwdop2 = 0: X, Z interpolation with iorder = 4
c    iwdop2 = 1: X, Z interpolation with iorder = 6
c
c  Flag (added 7/6/02) for table files:
c
c  ioptab =  0: Kurucz + OPAL92
c  ioptab =  1: Kurucz + OPAL95
c  ioptab = 11: Alexander + OPAL95
c
c  Modified 7/6/02, changing output to diagnistic string copcvr
c  (note that this was apparently incorrect previously)
c
c  Modified 12/4/06, to control iorder with iwdop2.
c  
      implicit double precision (a-h, o-z)
      character*80 file, filess
      character*256 tabnam, copcvr
      common/opcver/ copcvr
c  Note: the following commons seems never to be used
      common/copcgh/ imodgh,ioptbl
      common/cofile/ nfiles, idsfil(20), file(20), iopen(20), filess(20)
      common/cversn/ ivrbcs,ivrrhs,ivreos,ivropc,ivreng
      common/opatab/ ioptab
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
c     data ioptab /-1/
c
c  Hard-code version number in ivropc
c
      ivropc=9
c
      copcvr = 
     *  ' GH interpolation, version v9'
c
c  decode iwdopc
c
      iwdop0=mod(iwdopc,10)
      iwdop1=mod(iwdopc/10,10)
      iwdop2=mod(iwdopc/100,10)
      if(istdpr.gt.0) 
     *  write(istdpr,*) 'In opingh, iwdopc, iwdop0, iwdop1, iwdop2 =',
     *  iwdopc, iwdop0, iwdop1, iwdop2
c
      if(iwdop0.ne.8) then
        imodgh = 1
      else
        imodgh = 2
      end if
c
c  test for electron conduction
c
      if(iwdop1.eq.0) imodgh=-imodgh
      if(istdpr.gt.0) write(istdpr,*) 'In opingh, imodgh =', imodgh
c
c  test for order of X, Z interpolation
c
      if(iwdop2.eq.0) then
	iorder = 4
      else
	iorder = 6
      end if
      if(istdpr.gt.0) write(istdpr,*) 'In opingh, iorder =', iorder
c
c  set file name for tables
c
      call stfile(inopc,ntab)
      if(istdpr.gt.0) 
     *  write(istdpr,*) 'inopc, ntab', inopc, ntab
      tabnam=file(ntab)
      if(istdpr.gt.0) 
     *  write(istdpr,*) 'Now tabnam set to',tabnam
      close(inopc)
c
      call maceps(epsprc)
      call opinit(epsprc,iorder,tabnam,imodgh)
c
c  set version and table
c
      lcopcv = length(copcvr)
      if(ioptab.eq.0) then
	copcvr=copcvr(1:lcopcv)//'. OPAL92, Kurucz91 tables'
      else if(ioptab.eq.1) then
	copcvr=copcvr(1:lcopcv)//'. OPAL95, Kurucz91 tables'
      else if(ioptab.eq.11) then
	copcvr=copcvr(1:lcopcv)//'. OPAL95, Alexander94 tables'
      else
	copcvr=copcvr(1:lcopcv)//'. Tables undefined'
      end if
      return
      end
