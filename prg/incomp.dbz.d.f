      subroutine incomp(y,nn,iy,aztst,z,iche3,icnocs,iheccs, ihe3bc,
     *  istcno,isthec,icnotr,ihectr,icntsl,itlr,nt)
c
c  Initializes array cvr in common/compvr/, based on trial model 
c  and possibly variables in common/compsz/. Initialization
c  depends on values of parameters iche3 and icnocs,
c  as well as on the possible inclusion of diffusion,
c  flagged by idiffus (passed in common/cdiffu/).
c  Also may initialize storage parameters in common/cmpstr/.
c  May initialize heavy-element and CNO abundances in array y.
c  Initializes xtlcno in common/cnofrc/.
c  Initializes CNO abundances in aztst and in xcnoc.
c
c  Assumes that Z has already been set in common/heavy/ by
c  call of s/r sethvz.
c
c  On input, icnocs and iheccs define treatment of CNO cycle and 4He
c  burning, and icnotr and ihectr define the cases in the trial model
c
c  If nt .lt. 0, no setting of abundances (modification 5/5/95)
c
c  Currently only handles He and heavy-element diffusion.
c
c  Original version: 5/8/91.
c
c  Modified 19/7/95, to allow setting of heavy-element abundance in
c  common/heavy/. Note that for consistency with older versions,
c  the CNO abundances for istcno = 1 are still based on value of
c  z passed in argument list, whereas setting of Y uses Z as 
c  set in common/heavy/
c
c  Modified 5/8/95, to allow setting storage parameters in 
c  case of diffusion.
c
c  Modified 14/7/96, to set also He3 abundance from xzer3
c
c  Modified 13/10/98, to keep unchanged He3 abundance if xzer3 .lt. 0
c
c  Modified 20/7/01, to incorporate icnocs .gt. 10 (with full treatment
c  of CNO elements). So far, only icnocs = 11 has been implemented.
c
c  Modified 13/8/02, preparing for inclusion of 4He burning. This also
c  includes substantial revision of the storage in common/compvr/
c
      implicit double precision(a-h, o-z)
c
c  include statement setting storage indices for reaction rate
c  and energy generation parameters
c  Note: engenr.n.d.incl replaced by engenr.bz.d.incl, 8/8/02
      include 'engenr.bz.d.incl'
c
      dimension y(iy,nn), aztst(1)
      common/rhcn/ a1,a2,a3,a4,zdummy,nvar,ifwrt,irhtst,inentc,
     *  ii1,ii2,ii3,icomp
      common/rcncns/ a(knucst),b(knucst),s(knucst,isnuc),
     *  zz(knucst),zz12(2,knucst),zsmh,acno,awght(knucwg),q(knucq), 
     *  knucpm, kqvals, kawght
      common/compvr/ cvr(icvrmx,1)
      common/compsz/ xzerh, yzer, xzer3, xrz12, xrz13, xrz14, xrz16
      common/cmpstr/ ispxx3, ispcno, icvxx3, icvcno, ibcxx3, ibccno,
     *  isphec, idcomp, iccomp
      common/cstcvr/ icvh1, icvhe4, icvhe3, icvc12, icvc13, icvn14, 
     *  icvo16
      common/cstyst/ iyh1, iyhe3, iyc12, iyc13, iyn14, iyo16, 
     *  iyche4, iycc12
      common/cstazt/ iah1, iahe3, iac12, iac13, ian14, iao16, 
     *  iache4, iacc12
      common/cdiffu/ idiffus, itbdif, rhob, rmxdif, tlfdif, dtldif,
     *  ctdfmx, rctdif, rctdf1
      common/cnofrc/ fcno, xtlcno
      common/heavy/ zatmos, zhc, zh(1)
      common/bccomp/ xhc,xhs,xcnoc(nspcmx),xhecc(2)
      common/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng,idghe3
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      idiag=2
c
      if(idiag.ge.1.and.istdpr.gt.0) 
     *  write(istdpr,105) 
     *  nn,iche3, icnocs, ihe3bc, iheccs, istcno, icnotr, isthec, 
     *  ihectr, icntsl, itlr, idiffus
c
c  test for initializing storage indices and xtlcno
c  (only when not continuing iteration for L and R)
c
      if(itlr.le.1) then
c
c  zero cvr (might need a little care and checking)
c
	call zero(cvr,icvrmx*(nn+1))
c
	iyh1 = 5
c
c  initialize storage indices in y array  (reset below, 
c  depending on icnocs and iheccs)
c
        iyc12 = -1
	iyc13 = -1
        iyn14 = -1
	iyo16 = -1
	iyche4 = -1
	iycc12 = -1
c
c  initialize storage indices in cvr array  (reset below, 
c  depending on icnocs and iheccs)
c
        icvh1  = 1
	icvhe4 = 2
	icvhe3 = 3
c
	icvc12 = -1
	icvc13 = -1
	icvn14 = -1
        icvo16 = -1
c
c  initialize storage indices in aztst array  (reset below, 
c  depending on icnocs and iheccs)
c
	iah1 = 3
	iahe3 = 4
c
        iac12 = -1
	iac13 = -1
        ian14 = -1
	iao16 = -1
	iache4 = -1
	iacc12 = -1
c
	if(iche3.eq.1) then
	  ispxx3 = 2
	  iyhe3 = 6
        else
	  ispxx3 = 1
	  iyhe3 = -1
        end if
c
	icvxx3 = 2
c
	if(icnocs.eq.0) then
	  ispcno = 0
	  icvcno = 0
	else if(icnocs.eq.1) then
	  ispcno = 1
	  icvcno = 2
	  iyn14 = ispxx3 + 5
	  icvn14 = 6
	  icvo16 = 7
	  ian14 = 5
	else if(icnocs.eq.2) then
	  ispcno = 1
	  icvcno = 4
	  iyn14 = ispxx3 + 5
	  icvc12 = 4
	  icvc13 = 5
	  icvn14 = 6
	  icvo16 = 7
	  ian14 = 5
	else if(icnocs.eq.3) then
	  ispcno = 2
	  icvcno = 3
	  iyc13 = ispxx3 + 5
	  iyn14 = ispxx3 + 6
	  icvc12 = 4
	  icvc13 = 5
	  icvn14 = 6
	  iac13 = 5
	  ian14 = 6
	else if(icnocs.eq.4) then
	  ispcno = 3
	  icvcno = 4
	  iyc12 = ispxx3 + 5
	  iyc13 = ispxx3 + 6
	  iyn14 = ispxx3 + 7
	  icvc12 = 4
	  icvc13 = 5
	  icvn14 = 6
	  icvo16 = 7
	  iac12 = 5
	  iac13 = 6
	  ian14 = 7
	else if(icnocs.eq.11) then
	  ispcno = 2
	  icvcno = 2
	  iyn14 = ispxx3 + 5
	  iyo16 = ispxx3 + 6
	  icvn14 = 6
	  icvo16 = 7
	  ian14 = 5
	  iao16 = 6
        else
	  write(istdou,110) icnocs
	  if(istdpr.ne.istdou.and.istdpr.gt.0) write(istdpr,110) icnocs
	  stop 'incomp'
        end if
c
	if(iheccs.ne.0) then
	  isphec = 2
        else
	  isphec = 0
        end if
c
	if(ihe3bc.eq.0) then
	  ibcxx3 = 2
        else
	  ibcxx3 = 1
        end if
c
	ibccno = ispcno
c
c  test for setting storage for diffusion
c
	if(idiffus.eq.1) then
	  idcomp = 1
	  iccomp = ispxx3 + ispcno + isphec - 1
	else if(idiffus.eq.2) then
	  idcomp = 2
	  iccomp = ispxx3 + ispcno + isphec - 1
	else
	  idcomp = 0
	  iccomp = ispxx3 + ispcno + isphec
        end if
c
c
	if(icnocs.eq.1.or.icnocs.eq.11) then
	  xtlcno = fcno*z/acno + xrz16*z/awght(6)
c
c  reset N14 abundance, for consistency with old programme
c
	  xrz14 = fcno*awght(4)/acno
	else if(icnocs.eq.1.or.icnocs.eq.2.or.icnocs.eq.4) then
	  xtlcno = z*(xrz12/awght(2) + xrz13/awght(3) + xrz14/awght(4)
     *              + xrz16/awght(6))
	else if(icnocs.eq.3) then
	  xtlcno = z*(xrz12/awght(2) + xrz13/awght(3) + xrz14/awght(4))
        end if
c
	isphyd = ispxx3 + ispcno
c
c  test for setting 4He burning
c
	if(iheccs.ne.0) then
	  icvc12 = 4
          iyche4 = isphyd + 5
	  iycc12 = isphyd + 6
	  iache4 = ispcno + 5
	  iacc12 = ispcno + 6
        end if
c
        if(idiag.ge.1.and.istdpr.gt.0) then
          write(istdpr,115) ispxx3, ispcno, isphec, icvxx3, icvcno,
     *      ibcxx3, ibccno,xtlcno
          write(istdpr,116) idcomp, iccomp
        end if
c
      end if
c
      icomp = ispxx3 + ispcno + isphec
c
      if(nt.lt.0) then
	if(istdpr.gt.0) write(istdpr,117)
	return
      end if
c
c  determine case for setting of cvr
c
      if(itlr.gt.1) then
	if(nt.gt.0) then
	  isetcn = 1
	  isethc = 1
        else
	  isetcn = 0
	  isethc = 0
        end if
      else if(icntsl.ge.1) then
	if(icnotr.eq.icnocs) then
	  isetcn = 1
        else
	  isetcn = 2
	  if(istdpr.gt.0) write(istdpr,120) icnotr, icnocs
	end if
	if(ihectr.eq.iheccs) then
	  isethc = 1
        else
	  isethc = 2
	  if(istdpr.gt.0) write(istdpr,125) ihectr, iheccs
	end if
      else 
	if(istcno.eq.1) then
	  if(icnotr.eq.icnocs) then
	    isetcn = 1
          else
	    isetcn = 2
	    if(istdpr.gt.0) write(istdpr,130) icnotr, icnocs
	  end if
        else
          isetcn = 2
	end if
	if(isthec.eq.1) then
	  if(ihectr.eq.iheccs) then
	    isethc = 1
          else
	    isethc = 2
	    if(istdpr.gt.0) write(istdpr,135) ihectr, iheccs
	  end if
        else
          isethc = 2
	end if
      end if
c
      if(idiag.ge.1.and.istdpr.gt.0) write(istdpr,140) isetcn, isethc
c
c  set cvr, depending on case defined in isetcn, and on icnocs
c
c#ai#  Note: setting of He3 requires thought, although it is
c#ai#  unlikely to be of any significance.
c
      do 30 n=1,nn
      cvr(1,n) = y(5,n)
      cvr(2,n) = 1 - y(5,n) - zh(n)
      if(iche3.eq.1) then
	cvr(3,n) = y(6,n)
      else
	cvr(3,n) = 0
      end if
   30 continue
c
c  store Z in appropriate place in y, for heavy-element diffusion
c  **** The logics of this needs a little thought
c
c  For subsequent steps in (R, L) iteration, Z is assumed to
c  be read in with trial model and hence must be reset in zh,
c  for ZAMS model
c
      if(idiffus.eq.2) then
	iz=6+iccomp
	if(itlr.le.1) then
	  do 32 n=1,nn
   32     y(iz,n)=zh(n)
	else
	  do 34 n=1,nn
   34     zh(n)=y(iz,n)
	  zatmos=y(iz,1)
	  zhc=y(iz,nn)
        end if
      end if
c
      cvr(1,nn+1)=xhc
      cvr(2,nn+1)=1-zhc-xhc
c
      if(isetcn.eq.1) then
c
c  set variables from values already in y and aztst
c  use of xtlcno requires testing
c
	if(icnocs.eq.1) then
	  do 35 n=1,nn
          cvr(icvn14,n) = y(iyn14,n)
   35     cvr(icvo16,n) = (xtlcno-cvr(icvn14,n)/awght(4))*awght(6)
	  cvr(icvn14,nn+1) = aztst(ian14)
          cvr(icvo16,nn+1) = (xtlcno-cvr(icvn14,nn+1)/awght(4))*awght(6)
	else if(icnocs.eq.2) then
	  do 40 n=1,nn
   40     cvr(icvn14,n) = y(iyn14,n)
	  cvr(icvn14,nn+1) = aztst(ian14)
	else if(icnocs.eq.3) then
	  do 45 n=1,nn
          cvr(icvc13,n) = y(iyc13,n)
   45     cvr(icvn14,n) = y(iyn14,n)
	  cvr(icvc13,nn+1) = aztst(iac13)
	  cvr(icvn14,nn+1) = aztst(ian14)
	else if(icnocs.eq.4) then
	  do 47 n=1,nn
          cvr(icvc12,n) = y(iyc12,n)
          cvr(icvc13,n) = y(iyc13,n)
   47     cvr(icvn14,n) = y(iyn14,n)
	  cvr(icvc12,nn+1) = aztst(iac12)
	  cvr(icvc13,nn+1) = aztst(iac13)
	  cvr(icvn14,nn+1) = aztst(ian14)
	else if(icnocs.eq.11) then
	  do 50 n=1,nn
          cvr(icvn14,n) = y(iyn14,n)
   50     cvr(icvo16,n) = y(iyo16,n)
	  cvr(icvn14,nn+1) = aztst(ian14)
	  cvr(icvo16,nn+1) = aztst(iao16)
        end if
c
	if(icnocs.ge.1) then
	  call store(aztst(5),xcnoc(1),ispcno)
        end if
c
      else
c
c  set variables from values read in, as set in common/compsz/
c
	if(iche3.eq.1.and.icntsl.le.0) then
	  if(xzer3.ge.0) then
	    do 52 n=1,nn
	    y(6,n)=xzer3
   52       cvr(3,n)=xzer3
c
	    aztst(iahe3)=xzer3
	  else
	    do 53 n=1,nn
   53       cvr(3,n)=y(6,n)
          end if
        end if
c
	if(icnocs.eq.1) then
          xz14 = xrz14*z
          xz16 = xrz16*z
	  do 55 n=1,nn
	  y(iyn14,n) = xz14
          cvr(icvn14,n) = xz14
   55     cvr(icvo16,n) = xz16
	  aztst(ian14) = xz14
	  xcnoc(1) = xz14
	  cvr(icvn14,nn+1) = xz14
	  cvr(icvo16,nn+1) = xz16
	else if(icnocs.eq.2) then
          xz12 = xrz12*z
          xz13 = xrz13*z
          xz14 = xrz14*z
          xz16 = xrz16*z
	  do 60 n=1,nn
	  y(iyn14,n) = xz14
          cvr(icvc12,n) = xz12
          cvr(icvc13,n) = xz13
          cvr(icvn14,n) = xz14
   60     cvr(icvo16,n) = xz16
	  aztst(ian14) = xz14
	  xcnoc(1) = xz14
	  cvr(icvc12,nn+1) = xz12
	  cvr(icvc13,nn+1) = xz13
	  cvr(icvn14,nn+1) = xz14
	  cvr(icvo16,nn+1) = xz16
	else if(icnocs.eq.3) then
          xz12 = xrz12*z
          xz13 = xrz13*z
          xz14 = xrz14*z
	  do 65 n=1,nn
	  y(iyc13,n) = xz13
	  y(iyn14,n) = xz14
          cvr(icvc12,n) = xz12
          cvr(icvc13,n) = xz13
   65     cvr(icvn14,n) = xz14
	  aztst(iac13) = xz13
	  aztst(ian14) = xz14
	  xcnoc(1) = xz13
	  xcnoc(2) = xz14
	  cvr(icvc12,nn+1) = xz12
	  cvr(icvc13,nn+1) = xz13
	  cvr(icvn14,nn+1) = xz14
	else if(icnocs.eq.4) then
          xz12 = xrz12*z
          xz13 = xrz13*z
          xz14 = xrz14*z
          xz16 = xrz16*z
	  do 70 n=1,nn
	  y(iyc12,n) = xz12
	  y(iyc13,n) = xz13
	  y(iyn14,n) = xz14
          cvr(icvc12,n) = xz12
          cvr(icvc13,n) = xz13
          cvr(icvn14,n) = xz14
   70     cvr(icvo16,n) = xz16
	  aztst(iac12) = xz12
	  aztst(iac13) = xz13
	  aztst(ian14) = xz14
	  xcnoc(1) = xz12
	  xcnoc(2) = xz13
	  xcnoc(3) = xz14
	  cvr(icvc12,nn+1) = xz12
	  cvr(icvc13,nn+1) = xz13
	  cvr(icvn14,nn+1) = xz14
	  cvr(icvo16,nn+1) = xz16
	else if(icnocs.eq.11) then
          xz14 = xrz14*z
          xz16 = xrz16*z
	  do 72 n=1,nn
	  y(iyn14,n) = xz14
	  y(iyo16,n) = xz16
          cvr(icvn14,n) = xz14
   72     cvr(icvo16,n) = xz16
	  aztst(ian14) = xz14
	  xcnoc(1) = xz14
	  cvr(icvn14,nn+1) = xz14
	  cvr(icvo16,nn+1) = xz16
        end if
      end if
c
c  test for setting variables for 4He burning
c
      if(iheccs.ne.0) then
	if(isethc.eq.1) then
c
c  possibly overwrite variables in cvr
c
	  do 75 n=1,nn
	  cvr(2,n)=y(iyche4,n)
   75     cvr(icvc12,n)=y(iycc12,n)
	  cvr(2,nn+1)=aztst(iache4)
	  cvr(icvc12,nn+1)=aztst(iacc12)
        else
	  do 80 n=1,nn
	  y(iyche4,n)=cvr(2,n)
          y(iycc12,n)=cvr(icvc12,n)
   80     continue
	  aztst(iache4)=cvr(2,nn+1)
          aztst(iacc12)=cvr(icvc12,nn+1)
        end if
	call store(aztst(iache4),xhecc,2)
      end if
c
      if(idiag.ge.2.and.istdpr.gt.0.and.iycc12.gt.0.and.icvc12.gt.0) 
     *  then
	write(istdpr,150) (n,y(iyche4,n),y(iycc12,n),
     *    cvr(icvhe4,n),cvr(icvc12,n),n=1,nn,10)
	write(istdpr,155) aztst(iache4),aztst(iacc12),
     *    xhecc(1),xhecc(2)
      end if
c
      return
  105 format(//' Enter s/r incomp with nn =',i5,' iche3 =',i3,
     *  ' icnocs =',i3,' ihe3bc =',i3,' iheccs =',i3/
     *  ' istcno =',i3,' icnotr =',i3,' isthec =',i3,' ihectr =',i3,
     *  ' icntsl =',i3,' itlr =',i3/
     *  ' idiffus =',i3)
  110 format(//' ***** Error in s/r incomp. icnocs =',i3,
     *  '  not allowed')
  115 format(//' Quantities set in s/r incomp:'/
     *  ' ispxx3 =',i4/
     *  ' ispcno =',i4/
     *  ' isphec =',i4/
     *  ' icvxx3 =',i4/
     *  ' icvcno =',i4/
     *  ' ibcxx3 =',i4/
     *  ' ibccno =',i4/
     *  ' xtlcno =',1pe13.5)
  116 format(/
     *  ' idcomp =',i4/
     *  ' iccomp =',i4)
  117 format(/' Exit incomp without resetting abundances.')
  120 format(//' ***** Warning in s/r incomp. For icntsl = 1,',
     *  ' icnotr = ',i3, ' .ne. icnocs =',i3/
     *  ' CNO abundances set from values read in')
  125 format(//' ***** Warning in s/r incomp. For icntsl = 1,',
     *  ' ihectr = ',i3, ' .ne. iheccs =',i3/
     *  ' 4He and 12C abundances set from values read in')
  130 format(//' ***** Warning in s/r incomp. For istcno = 1,',
     *  ' icnotr = ',i3, ' .ne. icnocs =',i3/
     *  ' CNO abundances set from values read in')
  135 format(//' ***** Warning in s/r incomp. For isthec = 1,',
     *  ' ihectr = ',i3, ' .ne. iheccs =',i3/
     *  ' 4He and 12C abundances set from values read in')
  140 format(/' Setting abundances with isetcn =',i3,' isethc =',i3)
  150 format(//' Abundances set in s/r incomp.'/
     *  ' n,y(iyche4,n),y(iycc12,n),',
     *  ' cvr(icvhe4,n),cvr(icvc12,n):'/(i5,1p4e15.7))
  155 format(/' aztst(iache4), aztst(iacc12), xhecc(1), xhecc(2):'/
     *  1p4e13.5)
      end
